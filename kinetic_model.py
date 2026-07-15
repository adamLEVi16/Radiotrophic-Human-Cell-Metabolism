"""
Radiotrophic Cell ODE Kinetic Model (Phase 3, dose-driven rewrite)
==================================================================
Dynamic model of intracellular ROS during radiation exposure. Unlike the
constraint-based model, this resolves transients: the ROS spike at onset,
approach to steady state, and recovery after the beam stops.

What changed and why
--------------------
The previous version was driven by an abstract "radio_flux" whose ROS source
term (normalized by 120 nmol/min / 60) was ~6 orders of magnitude below the
scavenging capacity of a single SOD pool. As a result ROS never rose above
its initial condition: every dose and every defense knockout produced
identical, flat output. The model could not answer the question it was built
for.

This rewrite drives ROS production from absorbed dose rate (Gy/s) via water-
radiolysis G-values (radiotrophic_common), so production is physically scaled
and the system actually accumulates ROS, saturates enzymes, and discriminates
between defense configurations. Dsup's protection is capped at its measured
~40% efficacy (Hashimoto et al. 2016) rather than being allowed to remove all
hydroxyl radicals.

State variables:  [O2S, H2O2, GSH, GSSG, DNA_damage]
  Hydroxyl radical is treated in quasi-steady state (lifetime << 1 ns) and
  reported from its production/scavenging balance.

Kinetic parameters from literature (see inline citations).

Authors: Adam Labban
Date: March 2026
"""

import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd

import radiotrophic_common as rc

# ============================================================
# KINETIC PARAMETERS (literature-sourced)
# ============================================================
SOD_CONC = 10e-6        # M, Cu/Zn-SOD1 (Crapo et al. 1992)
CAT_CONC = 1e-6         # M, catalase (mostly peroxisomal, low cytosolic)
GPX_CONC = 0.5e-6       # M, GPX1 (selenium-dependent)
GR_CONC = 0.2e-6        # M, glutathione reductase

SOD_KCAT = 2e9          # 1/M/s, diffusion-limited (McCord & Fridovich 1969)
CAT_KCAT = 4e7          # 1/s (Ogura 1955)
CAT_KM = 1.1            # M (very high Km)
GPX_KCAT = 5e2          # 1/s
GPX_KM = 35e-6          # M (Flohe 1971)
GR_KCAT = 200           # 1/s (Carlberg & Mannervik 1975)
GR_KM = 65e-6           # M

# Pseudo-first-order rate constants (1/s) at physiological enzyme levels
K_SOD = SOD_KCAT * SOD_CONC             # ~2e4 /s
K_CAT = CAT_KCAT * CAT_CONC / CAT_KM    # ~36 /s (catalase is slow at low [H2O2])
K_SPONT_O2S = 1e5                       # 1/M/s spontaneous dismutation (2nd order)
VMAX_GPX = GPX_KCAT * GPX_CONC          # M/s
VMAX_GR = GR_KCAT * GR_CONC             # M/s

# Fenton chemistry (Fe2+ + H2O2 -> *OH)
FENTON_K = 76           # 1/M/s (Walling 1975)
FE2_CONC = 1e-6         # M labile iron (Kakhlon & Cabantchik 2002)

# Dsup: measured DNA-damage reduction (Hashimoto et al. 2016)
DSUP_EFFICIENCY = 0.40  # intercepts 40% of OH that would reach chromatin

# Mn-antioxidant (Daly et al. 2004; Barnese et al. 2012; Archibald 1982)
MN_CONC = 0.5e-3
MN_K_O2S = 1e6          # 1/M/s
MN_K_H2O2 = 1e3         # 1/M/s
K_MN_O2S = MN_K_O2S * MN_CONC           # ~500 /s
K_MN_H2O2 = MN_K_H2O2 * MN_CONC         # ~0.5 /s

# Nrf2-enhanced GR (Lewis et al. 2015)
NRF2_GR_BOOST = 2.0

# Glutathione pool (Meister 1988)
GSH_TOTAL = 5e-3

# Partition of hydroxyl radicals
PHI_DNA = 0.02          # fraction of *OH reaching the DNA backbone
FRAC_GSH_SCAV = 0.5     # fraction of the scavenged remainder handled by GSH

# Basal mitochondrial superoxide leak
BASAL_O2S = 0.5e-9      # M/s

# Repair
K_REPAIR = 0.01         # 1/s first-order BER (~100 s half-life)

# Cellular lethal thresholds
LETHAL_O2S = 50e-6      # M
LETHAL_H2O2 = 100e-6    # M

# Dose rate for the onset/ablation/pulse experiments. Set deliberately in the
# FLASH-radiotherapy regime (>40 Gy/s, Favaudon et al. 2014) — far above any
# survivable chronic exposure — so that radiolytic ROS actually stresses the
# defense network and each system's role becomes visible. The energy-budget
# analysis already shows that at survivable dose rates radiolytic ROS is
# negligible and the native network copes trivially.
STRESS_DOSE = 100.0     # Gy/s


def radiotrophic_ode(t, y, dose_rate, use_dsup=True, use_mn=True, use_nrf2=True,
                     use_sod=True, use_cat=True, use_gpx=True, radiation_on=True):
    """ODE system for transient ROS dynamics driven by absorbed dose rate."""
    O2S, H2O2, GSH, GSSG, DNA_dmg = (max(v, 0.0) for v in y)

    # --- ROS production from radiolysis (only while irradiated) ---
    if radiation_on:
        p = rc.ros_production_rates(dose_rate)
        prod_o2s, prod_h2o2, prod_oh = p["o2s"], p["h2o2"], p["oh"]
    else:
        prod_o2s = prod_h2o2 = prod_oh = 0.0

    # --- Superoxide handling ---
    v_sod = K_SOD * O2S if use_sod else 0.0
    v_mn_o2s = K_MN_O2S * O2S if use_mn else 0.0
    v_spont = K_SPONT_O2S * O2S * O2S            # always present

    dO2S = prod_o2s + BASAL_O2S - v_sod - v_mn_o2s - v_spont

    # --- Hydrogen peroxide handling ---
    v_cat = K_CAT * H2O2 if use_cat else 0.0
    v_gpx = (VMAX_GPX * H2O2 / (GPX_KM + H2O2)) * (GSH / (GSH + 1e-4)) if use_gpx else 0.0
    v_mn_h2o2 = K_MN_H2O2 * H2O2 if use_mn else 0.0
    v_fenton = FENTON_K * FE2_CONC * H2O2

    dH2O2 = (0.5 * v_sod + 0.5 * v_spont + prod_h2o2
             - v_cat - v_gpx - v_mn_h2o2 - v_fenton)

    # --- Hydroxyl radicals (quasi-steady) and DNA damage ---
    oh_flux = prod_oh + v_fenton                 # M/s of *OH generated
    dsup_factor = (1.0 - DSUP_EFFICIENCY) if use_dsup else 1.0
    dna_prod = oh_flux * PHI_DNA * dsup_factor
    scav_flux = oh_flux * (1.0 - PHI_DNA)
    v_oh_gsh = scav_flux * FRAC_GSH_SCAV

    dDNA = dna_prod - K_REPAIR * DNA_dmg

    # --- Glutathione redox cycling ---
    gr_boost = NRF2_GR_BOOST if use_nrf2 else 1.0
    v_gr = gr_boost * VMAX_GR * GSSG / (GR_KM + GSSG)
    pool_excess = (GSH + 2 * GSSG) - GSH_TOTAL
    pool_correction = 0.1 * pool_excess

    dGSH = 2 * v_gr - 2 * v_gpx - v_oh_gsh - pool_correction
    dGSSG = v_gpx + 0.5 * v_oh_gsh - v_gr

    return [dO2S, dH2O2, dGSH, dGSSG, dDNA]


def run_simulation(dose_rate, duration=300, pulse_off_time=None, **defenses):
    """Run one kinetic simulation. `defenses` overrides use_* flags."""
    y0 = [1e-9, 1e-7, GSH_TOTAL * 0.95, GSH_TOTAL * 0.05, 0.0]

    def wrapper(t, y):
        rad_on = not (pulse_off_time is not None and t > pulse_off_time)
        return radiotrophic_ode(t, y, dose_rate, radiation_on=rad_on, **defenses)

    t_eval = np.linspace(0, duration, 1000)
    sol = solve_ivp(wrapper, [0, duration], y0, t_eval=t_eval,
                    method="LSODA", rtol=1e-8, atol=1e-14, max_step=1.0)
    if not sol.success:
        print(f"  warning: ODE solver reported: {sol.message}")

    # Reconstruct quasi-steady OH concentration for reporting. Mirror the ODE's
    # radiation-on logic: after pulse_off_time the beam is off, so the radiolytic
    # OH source is zero and only the Fenton contribution remains. (Without this,
    # a pulsed run would report an on-beam OH source forever and mask recovery.)
    radiolytic_oh = rc.ros_production_rates(dose_rate)["oh"]
    oh_series = []
    for i in range(len(sol.t)):
        rad_on = not (pulse_off_time is not None and sol.t[i] > pulse_off_time)
        H2O2 = max(sol.y[1][i], 0)
        oh_flux = (radiolytic_oh if rad_on else 0.0) + FENTON_K * FE2_CONC * H2O2
        # nominal OH scavenging capacity (GSH-dominated) for a steady-state estimate
        k_oh_total = 1e10 * max(sol.y[2][i], 1e-6) + 1e9
        oh_series.append(oh_flux / k_oh_total)
    return pd.DataFrame({
        "time_s": sol.t,
        "superoxide_M": np.maximum(sol.y[0], 0),
        "h2o2_M": np.maximum(sol.y[1], 0),
        "oh_radical_M": np.array(oh_series),
        "gsh_M": np.maximum(sol.y[2], 0),
        "gssg_M": np.maximum(sol.y[3], 0),
        "dna_damage": np.maximum(sol.y[4], 0),
    })


def _run_with_overrides(dose_rate, duration, overrides, **defenses):
    """Run a simulation with selected module-level parameters temporarily
    overridden. The ODE reads these constants from the module namespace at
    call time, so setting them here changes the run and restoring afterward
    keeps the defaults intact for other experiments."""
    import sys
    mod = sys.modules[__name__]
    saved = {k: getattr(mod, k) for k in overrides}
    try:
        for k, v in overrides.items():
            setattr(mod, k, v)
        return run_simulation(dose_rate, duration=duration, **defenses)
    finally:
        for k, v in saved.items():
            setattr(mod, k, v)


def run_sensitivity_analysis(n=150, seed=7):
    """Monte Carlo uncertainty propagation for the kinetic model — the ODE
    analogue of the FBA Monte Carlo, so both arms carry formal uncertainty.

    The most uncertain kinetic parameters are sampled over literature-plausible
    ranges at the stress dose, and the resulting distributions of the harm
    readouts are reported, together with a sensitivity ranking (correlation of
    each parameter with accumulated DNA damage).

    Returns (distributions_df, ranking_df).
    """
    rng = np.random.default_rng(seed)
    sampled = {"K_SOD": [], "K_CAT": [], "VMAX_GPX": [],
               "FE2_CONC": [], "PHI_DNA": [], "K_REPAIR": []}
    out = {"steady_h2o2_uM": [], "dna_damage_120s": [], "gsh_depletion_pct": []}

    for _ in range(n):
        ov = {
            "K_SOD": K_SOD * rng.uniform(0.5, 2.0),      # SOD abundance +/- 2x
            "K_CAT": K_CAT * rng.uniform(0.3, 3.0),      # catalase (cytosolic, uncertain)
            "VMAX_GPX": VMAX_GPX * rng.uniform(0.4, 2.0),
            "FE2_CONC": rng.uniform(2e-7, 5e-6),         # labile iron pool varies widely
            "PHI_DNA": rng.uniform(0.005, 0.05),         # fraction of *OH reaching DNA
            "K_REPAIR": rng.uniform(0.005, 0.02),        # BER rate
        }
        df = _run_with_overrides(STRESS_DOSE, 120, ov)
        for k in sampled:
            sampled[k].append(ov[k])
        out["steady_h2o2_uM"].append(df["h2o2_M"].iloc[-1] * 1e6)
        out["dna_damage_120s"].append(df["dna_damage"].iloc[-1])
        out["gsh_depletion_pct"].append((1 - df["gsh_M"].iloc[-1] / GSH_TOTAL) * 100)

    dist_rows = []
    for name, vals in out.items():
        v = np.array(vals)
        dist_rows.append({
            "output": name,
            "p5": float(np.percentile(v, 5)),
            "median": float(np.median(v)),
            "p95": float(np.percentile(v, 95)),
        })
    dist = pd.DataFrame(dist_rows)
    for c in ("p5", "median", "p95"):
        dist[c] = dist[c].map(lambda x: float(f"{x:.4g}"))

    dna = np.array(out["dna_damage_120s"])
    rank_rows = []
    for k, vals in sampled.items():
        r = np.corrcoef(np.array(vals), dna)[0, 1]
        rank_rows.append({"parameter": k, "corr_with_dna_damage": round(float(r), 3)})
    rank_rows.sort(key=lambda r: -abs(r["corr_with_dna_damage"]))
    return dist, pd.DataFrame(rank_rows)


def run_all_kinetic_experiments():
    """Run all Phase 3 kinetic experiments."""
    results = {}

    # K1: radiation onset transient (all defenses on) --------------------
    print("  K1: Radiation onset transient...")
    results["k1_onset_transient"] = run_simulation(STRESS_DOSE, duration=300)

    # K2: dose-response (peak/steady ROS vs dose rate) -------------------
    print("  K2: Dose-response kinetics...")
    rows = []
    for dose in [0, 1e-3, 1e-2, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100]:
        df = run_simulation(dose, duration=120)
        peak_o2s, peak_h2o2 = df["superoxide_M"].max(), df["h2o2_M"].max()
        rows.append({
            "dose_Gy_s": dose,
            "peak_superoxide_uM": round(peak_o2s * 1e6, 5),
            "peak_h2o2_uM": round(peak_h2o2 * 1e6, 5),
            "peak_oh_nM": round(df["oh_radical_M"].max() * 1e9, 6),
            "steady_superoxide_uM": round(df["superoxide_M"].iloc[-1] * 1e6, 5),
            "steady_h2o2_uM": round(df["h2o2_M"].iloc[-1] * 1e6, 5),
            "dna_damage_120s": round(df["dna_damage"].iloc[-1], 6),
            "gsh_depletion_pct": round((1 - df["gsh_M"].iloc[-1] / GSH_TOTAL) * 100, 2),
            "superoxide_lethal": "YES" if peak_o2s > LETHAL_O2S else "no",
            "h2o2_lethal": "YES" if peak_h2o2 > LETHAL_H2O2 else "no",
        })
    results["k2_dose_response"] = pd.DataFrame(rows)

    # K3: defense ablation kinetics --------------------------------------
    print("  K3: Defense ablation kinetics...")
    configs = [
        ("All defenses",   dict()),
        ("No Dsup",        dict(use_dsup=False)),
        ("No Mn-AOX",      dict(use_mn=False)),
        ("No Nrf2",        dict(use_nrf2=False)),
        ("No SOD",         dict(use_sod=False)),
        ("No catalase",    dict(use_cat=False)),
        ("No GPX",         dict(use_gpx=False)),
        ("Native only",    dict(use_dsup=False, use_mn=False, use_nrf2=False)),
    ]
    # Reported at steady state (end of exposure); peak values are dominated by
    # the resting initial condition and are uninformative here.
    rows = []
    for label, d in configs:
        df = run_simulation(STRESS_DOSE, duration=120, **d)
        rows.append({
            "config": label,
            "steady_superoxide_uM": round(df["superoxide_M"].iloc[-1] * 1e6, 5),
            "steady_h2o2_uM": round(df["h2o2_M"].iloc[-1] * 1e6, 5),
            "dna_damage_120s": round(df["dna_damage"].iloc[-1], 8),
            "gsh_depletion_pct": round((1 - df["gsh_M"].iloc[-1] / GSH_TOTAL) * 100, 3),
        })
    results["k3_ablation_kinetics"] = pd.DataFrame(rows)

    # K4: radiation pulse recovery (beam off at 60 s) --------------------
    print("  K4: Radiation pulse recovery...")
    results["k4_pulse_recovery"] = run_simulation(STRESS_DOSE, duration=300,
                                                  pulse_off_time=60)

    # K5: parameter sensitivity / uncertainty (Monte Carlo) --------------
    print("  K5: Parameter sensitivity (Monte Carlo)...")
    dist, ranking = run_sensitivity_analysis()
    results["k5_sensitivity"] = dist
    results["k5_sensitivity_ranking"] = ranking

    return results


if __name__ == "__main__":
    print("Running Phase 3: ODE Kinetic Model (dose-driven)")
    print("=" * 60)
    results = run_all_kinetic_experiments()

    for name, df in results.items():
        print(f"\n{'='*60}\nEXPERIMENT: {name.upper()}\n{'='*60}")
        if len(df) <= 20:
            print(df.to_string(index=False))
        else:
            print(f"  Time series: {len(df)} points, {df['time_s'].iloc[-1]:.0f}s")
            print(f"  Peak superoxide: {df['superoxide_M'].max()*1e6:.4f} uM")
            print(f"  Peak H2O2:       {df['h2o2_M'].max()*1e6:.4f} uM")
            print(f"  Final DNA damage:{df['dna_damage'].iloc[-1]:.6f}")
            print(f"  Final GSH:       {df['gsh_M'].iloc[-1]*1e3:.2f} mM "
                  f"({df['gsh_M'].iloc[-1]/GSH_TOTAL*100:.1f}%)")

    for name, df in results.items():
        df.to_csv(f"{name}.csv", index=False)
    print("\n\nKinetic results saved to current directory")

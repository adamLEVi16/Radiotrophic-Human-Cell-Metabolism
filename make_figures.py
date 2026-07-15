"""
Figure generation
==================
Renders the analysis figures from the CSV outputs produced by
energy_budget.py, radiotrophic_model.py, and kinetic_model.py.
Run those first (or use run_all.py, which sequences everything).

Outputs PNG files to figures/.
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FIGDIR = "figures"
COL = {"radio": "#c1440e", "normal": "#3b6ea5", "accent": "#4a7c59",
       "warn": "#b5651d", "grid": "#d9d9d9"}


def _style(ax):
    ax.grid(True, color=COL["grid"], linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fig_energy_budget():
    df = pd.read_csv("energy_budget.csv")
    d = df[(df.efficiency == 1.0) & (df.atp_demand == "high_1e9")].copy()
    d = d.sort_values("dose_rate_Gy_s")
    fig, ax = plt.subplots(figsize=(8, 4.6))
    x = np.arange(len(d))
    ax.bar(x, d["radiation_pct_of_demand"], color=COL["radio"], width=0.6)
    ax.set_yscale("log")
    ax.axhline(100, color=COL["normal"], linestyle="--", linewidth=1.2,
               label="100% of ATP demand (break-even)")
    ax.set_xticks(x)
    ax.set_xticklabels([r.replace("_", "\n") for r in d["regime"]], fontsize=8)
    ax.set_ylabel("Radiation energy as % of ATP demand\n(log scale, 100% efficiency)")
    ax.set_title("Radiation cannot power a human cell at survivable dose rates")
    ax.legend(fontsize=8, frameon=False)
    _style(ax)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig1_energy_budget.png", dpi=140)
    plt.close(fig)


def fig_kinetic_dose_response():
    df = pd.read_csv("k2_dose_response.csv")
    d = df[df.dose_Gy_s > 0]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2))
    a1.loglog(d.dose_Gy_s, d.steady_h2o2_uM, "o-", color=COL["radio"])
    a1.set_xlabel("Dose rate (Gy/s)")
    a1.set_ylabel("Steady-state H$_2$O$_2$ (uM)")
    a1.set_title("H$_2$O$_2$ rises with dose (defenses intact)")
    _style(a1)
    a2.loglog(d.dose_Gy_s, d.dna_damage_120s, "s-", color=COL["accent"])
    a2.set_xlabel("Dose rate (Gy/s)")
    a2.set_ylabel("DNA damage at 120 s (a.u.)")
    a2.set_title("DNA damage accumulation vs dose")
    _style(a2)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig2_kinetic_dose_response.png", dpi=140)
    plt.close(fig)


def fig_kinetic_ablation():
    df = pd.read_csv("k3_ablation_kinetics.csv")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
    x = np.arange(len(df))
    a1.bar(x, df.steady_h2o2_uM, color=COL["radio"])
    a1.set_xticks(x); a1.set_xticklabels(df.config, rotation=40, ha="right", fontsize=8)
    a1.set_ylabel("Steady-state H$_2$O$_2$ (uM)")
    a1.set_title("Peroxide control (catalase/GPX dominate)")
    _style(a1)
    a2.bar(x, df.dna_damage_120s, color=COL["accent"])
    a2.set_xticks(x); a2.set_xticklabels(df.config, rotation=40, ha="right", fontsize=8)
    a2.set_ylabel("DNA damage at 120 s (a.u.)")
    a2.set_title("DNA protection (Dsup dominates)")
    _style(a2)
    fig.suptitle("Defense ablation at FLASH-regime dose (100 Gy/s)", y=1.02)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig3_kinetic_ablation.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_kinetic_timecourse():
    k1 = pd.read_csv("k1_onset_transient.csv")
    k4 = pd.read_csv("k4_pulse_recovery.csv")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2), sharey=False)
    a1.plot(k1.time_s, k1.h2o2_M * 1e6, color=COL["radio"], label="H$_2$O$_2$")
    a1.set_xlabel("Time (s)"); a1.set_ylabel("H$_2$O$_2$ (uM)")
    a1.set_title("K1: onset transient (beam on)")
    _style(a1)
    a2.plot(k4.time_s, k4.h2o2_M * 1e6, color=COL["warn"])
    a2.axvline(60, color=COL["normal"], linestyle="--", linewidth=1, label="beam off")
    a2.set_xlabel("Time (s)"); a2.set_ylabel("H$_2$O$_2$ (uM)")
    a2.set_title("K4: recovery after pulse")
    a2.legend(fontsize=8, frameon=False)
    _style(a2)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig4_kinetic_timecourse.png", dpi=140)
    plt.close(fig)


def fig_energy_constrained():
    df = pd.read_csv("energy_constrained_dose.csv")
    d = df[df.efficiency == 1.0].sort_values("dose_Gy_s")
    fig, ax = plt.subplots(figsize=(8, 4.4))
    x = np.arange(len(d))
    ax.bar(x, np.maximum(d.pct_boost, 1e-9), color=COL["radio"])
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([r.replace("_", "\n") for r in d.regime], fontsize=8)
    ax.set_ylabel("ATP boost from radiotrophy (%)\n(energy-constrained, log scale)")
    ax.set_title("With energy conservation enforced, the ATP boost vanishes")
    _style(ax)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig5_energy_constrained.png", dpi=140)
    plt.close(fig)


def fig_monte_carlo():
    # Recompute the raw distribution for a histogram (summary CSV has only stats).
    import radiotrophic_model as R
    rng = np.random.default_rng(42)
    boosts = []
    for _ in range(400):
        ros = rng.uniform(0.1, 1.0); atpo = rng.uniform(0.2, 1.0)
        sodc = rng.uniform(2, 6); melt = rng.uniform(0.01, 0.05)
        m = R.build_model(); radio = m.reactions.get_by_id("RADIO")
        mets = {mt.id: mt for mt in radio.metabolites}
        radio.add_metabolites({mets["o2s_c"]: ros, mets["atp_c"]: -atpo,
                               mets["adp_c"]: atpo, mets["pi_c"]: atpo,
                               mets["melanin_c"]: -melt}, combine=False)
        m.reactions.get_by_id("SODc").upper_bound = sodc
        m.reactions.get_by_id("EX_glc").lower_bound = -5
        ar = m.slim_optimize(); R.disable_engineered(m); an = m.slim_optimize()
        if an and an > 0 and ar is not None:
            boosts.append((ar - an) / an * 100)
    boosts = np.array(boosts)
    fig, ax = plt.subplots(figsize=(8, 4.4))
    ax.hist(boosts, bins=30, color=COL["normal"], edgecolor="white")
    ax.axvline(np.median(boosts), color=COL["radio"], linewidth=1.5,
               label=f"median {np.median(boosts):.1f}%")
    ax.set_xlabel("ATP boost, unconstrained RADIO (%)")
    ax.set_ylabel("Samples")
    ax.set_title("Monte Carlo: abstract-model boost spans 0-32% with parameter uncertainty")
    ax.legend(fontsize=8, frameon=False)
    _style(ax)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig6_monte_carlo.png", dpi=140)
    plt.close(fig)


def fig_kinetic_sensitivity():
    df = pd.read_csv("k5_sensitivity_ranking.csv")
    df = df.reindex(df.corr_with_dna_damage.abs().sort_values().index)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    colors = [COL["radio"] if v >= 0 else COL["normal"]
              for v in df.corr_with_dna_damage]
    ax.barh(df.parameter, df.corr_with_dna_damage, color=colors)
    ax.axvline(0, color="#555", linewidth=0.8)
    ax.set_xlabel("Correlation with accumulated DNA damage")
    ax.set_title("Kinetic sensitivity: DNA damage is driven by ·OH-to-DNA fraction and repair")
    _style(ax)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig7_kinetic_sensitivity.png", dpi=140)
    plt.close(fig)


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    for fn in (fig_energy_budget, fig_kinetic_dose_response, fig_kinetic_ablation,
               fig_kinetic_timecourse, fig_energy_constrained, fig_monte_carlo,
               fig_kinetic_sensitivity):
        print(f"  {fn.__name__} ...")
        fn()
    print(f"Figures written to {FIGDIR}/")


if __name__ == "__main__":
    main()

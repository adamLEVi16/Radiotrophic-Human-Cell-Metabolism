"""
Global Parameter Sensitivity Analysis
======================================
Sweeps ±50% on 6 uncertain parameters and re-runs three headline experiments:
  1. dsup_conditional sweep (threshold location — most important)
  2. Defense ablation (SOD bottleneck robustness)
  3. Dose-response ceiling (max ATP gain)

Key question for each: does the threshold finding hold, and does the threshold
LOCATION shift or hold across the parameter range?

Parameters swept:
  po_ratio_nadh       2.5  ATP per NADH (ETC_N stoichiometry)
  po_ratio_fadh2      1.5  ATP per FADH2 (ETC_F stoichiometry)
  etc_leak_nadh       0.01 superoxide per NADH in ETC_N (electron leak rate)
  melanin_coeff       0.02 melanin consumed per unit RADIO flux
  hashimoto_eff       0.40 Dsup protection ceiling (Hashimoto et al. 2016 40%)
  oh_scav_cap         10   OH_SCAV upper bound (central estimate, range 5-20)

Each parameter is tested at -50%, baseline, and +50% of its central value.
All other parameters are held at their baseline values (one-at-a-time design).

Output:
  sensitivity_summary.csv  — one row per (parameter, level) with headline metrics
  sensitivity_threshold.csv — full dsup_conditional curve at each (parameter, level)

Authors: Adam Labban
"""

import cobra
import pandas as pd
import sys
import os

# Import the build_model function from radiotrophic_model.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from radiotrophic_model import build_model, disable_engineered


# ============================================================
# PARAMETER DEFINITIONS
# ============================================================

BASELINE_PARAMS = {
    'po_ratio_nadh':   2.5,
    'po_ratio_fadh2':  1.5,
    'etc_leak_nadh':   0.01,
    'melanin_coeff':   0.02,
    'hashimoto_eff':   0.40,
    'oh_scav_cap':     10.0,
}

SWEEP_LEVELS = [-0.50, 0.0, +0.50]  # fraction change from baseline


def build_modified_model(params):
    """
    Build the radiotrophic model with modified parameter values.
    Applies changes directly to reaction stoichiometry / bounds.
    """
    model = build_model()

    # --- P/O ratio: NADH ---
    # ETC_N: nadh_m:-1, o2_m:-0.5, adp_m:-2.5, pi_m:-2.5,
    #        nad_m:1, h2o_m:0.5, atp_m:2.5, o2s_m:0.01
    etc_n = model.reactions.get_by_id('ETC_N')
    po_n = params['po_ratio_nadh']
    leak_n = params['etc_leak_nadh']
    atp_m = model.metabolites.get_by_id('atp_m')
    adp_m = model.metabolites.get_by_id('adp_m')
    pi_m  = model.metabolites.get_by_id('pi_m')
    o2s_m = model.metabolites.get_by_id('o2s_m')
    # Set new stoichiometry relative to baseline (delta from default)
    etc_n.add_metabolites({
        atp_m: po_n - 2.5,    # delta from default 2.5
        adp_m: -(po_n - 2.5),
        pi_m:  -(po_n - 2.5),
        o2s_m: leak_n - 0.01, # delta from default 0.01
    })

    # --- P/O ratio: FADH2 ---
    # ETC_F: fadh2_m:-1, o2_m:-0.5, adp_m:-1.5, pi_m:-1.5,
    #        fad_m:1, h2o_m:0.5, atp_m:1.5, o2s_m:0.005
    etc_f = model.reactions.get_by_id('ETC_F')
    po_f = params['po_ratio_fadh2']
    atp_m2 = model.metabolites.get_by_id('atp_m')
    adp_m2 = model.metabolites.get_by_id('adp_m')
    pi_m2  = model.metabolites.get_by_id('pi_m')
    etc_f.add_metabolites({
        atp_m2: po_f - 1.5,
        adp_m2: -(po_f - 1.5),
        pi_m2:  -(po_f - 1.5),
    })

    # --- Melanin consumption coefficient in RADIO ---
    radio = model.reactions.get_by_id('RADIO')
    mel_c = model.metabolites.get_by_id('melanin_c')
    mc = params['melanin_coeff']
    # Current stoichiometry: melanin_c: -0.02
    # New stoichiometry: melanin_c: -mc
    radio.add_metabolites({mel_c: -(mc - 0.02)})  # delta (more negative = more consumed)

    # --- Hashimoto efficiency: sets DSUP upper bound ---
    # DSUP cap = hashimoto_eff × (max_RADIO_flux × OH_stoichiometry)
    # = hashimoto_eff × (50 × 0.24) = hashimoto_eff × 12
    dsup_ub = params['hashimoto_eff'] * 50 * 0.24
    model.reactions.get_by_id('DSUP').upper_bound = dsup_ub

    # --- OH_SCAV upper bound ---
    model.reactions.get_by_id('OH_SCAV').upper_bound = params['oh_scav_cap']

    return model


def run_dsup_conditional(model, oh_scav_caps=None):
    """
    Run the dsup_conditional sweep on a given model.
    Returns DataFrame with threshold metrics.
    """
    if oh_scav_caps is None:
        oh_scav_caps = [0, 1, 2, 3, 5, 7, 10, 15, 20, 50, 100]

    dsup_ub = model.reactions.get_by_id('DSUP').upper_bound
    rows = []
    for cap in oh_scav_caps:
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -5
            model.reactions.get_by_id('RADIO').upper_bound = 50
            model.reactions.get_by_id('DSUP').upper_bound = dsup_ub
            model.reactions.get_by_id('OH_SCAV').upper_bound = cap
            sol = model.optimize()
            atp_with = sol.objective_value if sol.status == 'optimal' else 0

        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -5
            model.reactions.get_by_id('RADIO').upper_bound = 50
            model.reactions.get_by_id('DSUP').upper_bound = 0
            model.reactions.get_by_id('OH_SCAV').upper_bound = cap
            sol = model.optimize()
            atp_no = sol.objective_value if sol.status == 'optimal' else 0

        rows.append({
            'oh_scav_cap': cap,
            'atp_with_dsup': round(atp_with, 2),
            'atp_no_dsup': round(atp_no, 2),
            'atp_delta': round(atp_with - atp_no, 2),
            'dsup_essential': 'YES' if (atp_with - atp_no) > 1.0 else 'no'
        })
    return pd.DataFrame(rows)


def find_threshold(df_cond):
    """
    Find the OH_SCAV cap at which Dsup transitions from essential to dispensable.
    Returns the cap value just before delta drops to zero, or None if always essential/dispensable.
    """
    essential_caps = df_cond[df_cond['dsup_essential'] == 'YES']['oh_scav_cap'].tolist()
    dispensable_caps = df_cond[df_cond['dsup_essential'] == 'no']['oh_scav_cap'].tolist()
    if not essential_caps:
        return None, 'always_dispensable'
    if not dispensable_caps:
        return None, 'always_essential'
    return max(essential_caps), 'threshold_exists'


def run_ablation_atp(model):
    """Run the SOD-ablation case and return ATP with/without SOD."""
    with model:
        model.reactions.get_by_id('EX_glc').lower_bound = -5
        model.reactions.get_by_id('RADIO').upper_bound = 50
        sol = model.optimize()
        atp_all = sol.objective_value if sol.status == 'optimal' else 0

    with model:
        model.reactions.get_by_id('EX_glc').lower_bound = -5
        model.reactions.get_by_id('RADIO').upper_bound = 50
        model.reactions.get_by_id('SODc').upper_bound = 0
        model.reactions.get_by_id('SODm').upper_bound = 0
        sol = model.optimize()
        atp_no_sod = sol.objective_value if sol.status == 'optimal' else 0

    return round(atp_all, 2), round(atp_no_sod, 2)


def run_max_atp_gain(model):
    """Return the ATP boost at max RADIO flux vs normal cell."""
    with model:
        model.reactions.get_by_id('EX_glc').lower_bound = -5
        sol = model.optimize()
        atp_radio = sol.objective_value if sol.status == 'optimal' else 0

    with model:
        model.reactions.get_by_id('EX_glc').lower_bound = -5
        disable_engineered(model)
        sol = model.optimize()
        atp_norm = sol.objective_value if sol.status == 'optimal' else 0

    gain = atp_radio - atp_norm
    pct = gain / atp_norm * 100 if atp_norm > 0 else 0
    return round(atp_radio, 2), round(atp_norm, 2), round(gain, 2), round(pct, 1)


# ============================================================
# MAIN SWEEP
# ============================================================

if __name__ == '__main__':
    print("Global Parameter Sensitivity Analysis")
    print("=" * 60)
    print("Sweeping 6 parameters at -50%, baseline, +50%")
    print("Most important metric: Dsup threshold location\n")

    summary_rows = []
    threshold_rows = []

    for param_name, baseline_val in BASELINE_PARAMS.items():
        for level_frac in SWEEP_LEVELS:
            level_label = {-0.5: 'low (-50%)', 0.0: 'baseline', 0.5: 'high (+50%)'}[level_frac]
            param_val = baseline_val * (1 + level_frac)

            # Build params dict with this one parameter varied
            params = dict(BASELINE_PARAMS)
            params[param_name] = param_val

            print(f"  {param_name} = {param_val:.4f} ({level_label})...")

            try:
                model = build_modified_model(params)

                # Run dsup_conditional
                df_cond = run_dsup_conditional(model)
                threshold_cap, threshold_status = find_threshold(df_cond)

                # Run ablation (SOD bottleneck)
                atp_all, atp_no_sod = run_ablation_atp(model)
                sod_essential = atp_all > atp_no_sod + 1.0

                # Run max ATP gain
                atp_radio, atp_norm, gain, pct = run_max_atp_gain(model)

                summary_rows.append({
                    'parameter': param_name,
                    'baseline_value': baseline_val,
                    'swept_value': round(param_val, 4),
                    'level': level_label,
                    'dsup_threshold_oh_scav': threshold_cap,
                    'dsup_threshold_status': threshold_status,
                    'atp_all_defenses': atp_all,
                    'atp_no_sod': atp_no_sod,
                    'sod_essential': 'YES' if sod_essential else 'no',
                    'atp_radio': atp_radio,
                    'atp_normal': atp_norm,
                    'atp_gain': gain,
                    'pct_boost': pct,
                })

                # Store full threshold curve for this param/level
                df_cond['parameter'] = param_name
                df_cond['swept_value'] = round(param_val, 4)
                df_cond['level'] = level_label
                threshold_rows.append(df_cond)

            except Exception as e:
                print(f"    ERROR: {e}")
                summary_rows.append({
                    'parameter': param_name,
                    'baseline_value': baseline_val,
                    'swept_value': round(param_val, 4),
                    'level': level_label,
                    'dsup_threshold_oh_scav': None,
                    'dsup_threshold_status': f'ERROR: {e}',
                    'atp_all_defenses': None,
                    'atp_no_sod': None,
                    'sod_essential': None,
                    'atp_radio': None,
                    'atp_normal': None,
                    'atp_gain': None,
                    'pct_boost': None,
                })

    df_summary = pd.DataFrame(summary_rows)
    df_threshold = pd.concat(threshold_rows, ignore_index=True) if threshold_rows else pd.DataFrame()

    # ============================================================
    # PRINT RESULTS
    # ============================================================
    print("\n" + "=" * 60)
    print("SENSITIVITY SUMMARY")
    print("=" * 60)
    print(df_summary.to_string(index=False))

    print("\n" + "=" * 60)
    print("THRESHOLD ROBUSTNESS SUMMARY")
    print("=" * 60)
    threshold_summary = df_summary[['parameter', 'swept_value', 'level',
                                     'dsup_threshold_oh_scav', 'dsup_threshold_status',
                                     'sod_essential', 'pct_boost']].copy()
    print(threshold_summary.to_string(index=False))

    # ============================================================
    # KEY FINDINGS CHECK
    # ============================================================
    print("\n" + "=" * 60)
    print("KEY FINDINGS ROBUSTNESS")
    print("=" * 60)

    all_threshold_exist = all(
        r == 'threshold_exists' for r in df_summary['dsup_threshold_status']
        if r is not None
    )
    all_sod_essential = all(
        r == 'YES' for r in df_summary['sod_essential']
        if r is not None
    )
    threshold_locations = df_summary[df_summary['dsup_threshold_status'] == 'threshold_exists']['dsup_threshold_oh_scav']
    pct_boosts = df_summary['pct_boost'].dropna()

    print(f"  Dsup threshold exists across all parameter sweeps: {'YES' if all_threshold_exist else 'NO — see table'}")
    if len(threshold_locations) > 0:
        print(f"  Threshold location range: OH_SCAV {threshold_locations.min()}–{threshold_locations.max()} (baseline: 5)")
    print(f"  SOD essential across all sweeps: {'YES' if all_sod_essential else 'NO — see table'}")
    print(f"  ATP boost range: {pct_boosts.min():.1f}%–{pct_boosts.max():.1f}% (baseline: 5.5%)")

    # Save
    df_summary.to_csv('sensitivity_summary.csv', index=False)
    df_threshold.to_csv('sensitivity_threshold.csv', index=False)
    print("\nResults saved: sensitivity_summary.csv, sensitivity_threshold.csv")

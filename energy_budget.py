"""
Energy budget analysis (Phase 3)
================================
The load-bearing question the metabolic model cannot answer on its own:
does ionizing radiation carry enough energy to meaningfully power a human
cell? This module compares absorbed radiation power to ATP-demand power
across realistic dose regimes, using only conservation of energy.

The metabolic (FBA) model treats radiotrophic NADH generation as a source
bounded by an arbitrary flux cap. Here that source is bounded instead by the
energy actually deposited, which is the physically correct constraint. The
result determines whether the whole premise is viable.

Output: energy_budget.csv
"""

import pandas as pd
import radiotrophic_common as rc


def build_energy_budget():
    """Radiation power vs ATP demand across dose regimes and efficiencies."""
    rows = []
    for regime, dose in rc.DOSE_REGIMES.items():
        p_rad = rc.radiation_power_per_cell(dose)
        for atp_label, atp in (("low_1e7", rc.ATP_DEMAND_LOW),
                               ("high_1e9", rc.ATP_DEMAND_HIGH)):
            p_atp = rc.atp_demand_power(atp)
            for eff in (1.0, 0.01):   # thermodynamic ceiling and 1% (optimistic melanin)
                frac = rc.energy_fraction(dose, atp, eff)
                rows.append({
                    "regime": regime,
                    "dose_rate_Gy_s": dose,
                    "atp_demand": atp_label,
                    "efficiency": eff,
                    "radiation_power_W": p_rad,
                    "atp_demand_power_W": p_atp,
                    "radiation_fraction_of_demand": frac,
                    "radiation_pct_of_demand": frac * 100.0,
                })
    return pd.DataFrame(rows)


def build_breakeven_table():
    """Dose rates required for radiation to matter, vs lethality context."""
    LETHAL_TOTAL_DOSE = 5.0  # Gy, approximate whole-body LD50 (acute)
    rows = []
    for atp_label, atp in (("low_1e7", rc.ATP_DEMAND_LOW),
                           ("high_1e9", rc.ATP_DEMAND_HIGH)):
        for eff in (1.0, 0.01):
            d_be = rc.breakeven_dose_rate(atp, eff)
            rows.append({
                "atp_demand": atp_label,
                "efficiency": eff,
                "breakeven_dose_Gy_s": d_be,
                "seconds_to_lethal_5Gy": LETHAL_TOTAL_DOSE / d_be,
                "verdict": ("survivable" if LETHAL_TOTAL_DOSE / d_be > 3600
                            else "rapidly lethal"),
            })
    return rows


def build_flux_reconciliation():
    """Dose rate implied by the original FBA RADIO fluxes (energy-limited)."""
    rows = []
    for flux in (2, 5, 10, 25, 28.57, 50):
        for eff in (1.0, 0.01):
            rows.append({
                "radio_flux_mmol_gDW_h": flux,
                "efficiency": eff,
                "required_dose_Gy_s": rc.dose_for_radio_flux(flux, eff),
                "x_lethal_dose_rate": rc.dose_for_radio_flux(flux, eff) / 5.0,
            })
    return rows


def main():
    df = build_energy_budget()
    df.to_csv("energy_budget.csv", index=False)

    print("ENERGY BUDGET: radiation power vs ATP demand")
    print("=" * 70)
    view = df[(df["efficiency"] == 1.0)]
    for _, r in view.iterrows():
        print(f"  {r['regime']:22s} {r['atp_demand']:9s} "
              f"D={r['dose_rate_Gy_s']:.1e} Gy/s  "
              f"radiation = {r['radiation_pct_of_demand']:.2e}% of demand")

    print("\nBREAK-EVEN (radiation power = ATP demand)")
    print("-" * 70)
    for r in build_breakeven_table():
        print(f"  demand {r['atp_demand']:9s} eff {r['efficiency']:.2f}: "
              f"{r['breakeven_dose_Gy_s']:.2f} Gy/s -> lethal in "
              f"{r['seconds_to_lethal_5Gy']:.1f} s  [{r['verdict']}]")

    print("\nRECONCILIATION: dose rate implied by original FBA RADIO fluxes")
    print("-" * 70)
    for r in build_flux_reconciliation():
        if r["efficiency"] == 1.0:
            print(f"  RADIO flux {r['radio_flux_mmol_gDW_h']:6.2f} needs "
                  f"{r['required_dose_Gy_s']:8.1f} Gy/s "
                  f"({r['x_lethal_dose_rate']:.0f}x the whole-body-lethal dose rate)")

    print("\nCONCLUSION: at every survivable dose rate radiation supplies a")
    print("vanishing fraction of ATP demand; break-even requires promptly")
    print("lethal dose rates even at 100% conversion efficiency.")
    print("\nSaved: energy_budget.csv")


if __name__ == "__main__":
    main()

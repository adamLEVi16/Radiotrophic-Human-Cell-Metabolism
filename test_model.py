"""
Tests for the radiotrophic cell models
=======================================
Runs under pytest (`pytest test_model.py`) or standalone (`python test_model.py`).

Covers three things:
  1. Structural sanity of the FBA model (feasibility, conserved-moiety balance).
  2. That the key scientific results hold (Dsup cap enforced, energy constraint
     kills the boost at survivable doses, kinetic model discriminates defenses).
  3. A regression guard on the exact bug fixed during development: perturbing
     the radiotrophic reaction must keep ATP/ADP/Pi balanced.
"""

import numpy as np
import radiotrophic_common as rc
import radiotrophic_model as R
import kinetic_model as K


# --- 1. Structural sanity -------------------------------------------------
def test_model_builds_and_is_feasible():
    m = R.build_model()
    m.reactions.get_by_id("EX_glc").lower_bound = -5
    sol = m.optimize()
    assert sol.status == "optimal"
    assert sol.objective_value > 0


def test_radio_reaction_adenylate_balanced():
    # ATP consumed must equal ADP produced and Pi produced (no phantom energy).
    m = R.build_model()
    radio = m.reactions.get_by_id("RADIO")
    c = {mt.id: radio.metabolites[mt] for mt in radio.metabolites}
    assert abs(c["atp_c"] + c["adp_c"]) < 1e-9, "ATP/ADP not balanced in RADIO"
    assert abs(c["adp_c"] - c["pi_c"]) < 1e-9, "ADP/Pi not balanced in RADIO"


def test_atp_hydrolysis_reactions_balanced():
    m = R.build_model()
    for rid in ("ATPM", "BER"):
        c = {mt.id: m.reactions.get_by_id(rid).metabolites[mt]
             for mt in m.reactions.get_by_id(rid).metabolites}
        assert abs(c.get("atp_c", 0) + c.get("adp_c", 0)) < 1e-9, f"{rid} adenylate"


# --- 2. Key scientific results -------------------------------------------
def test_dsup_cap_enforced():
    # Dsup flux must not exceed 40% of the hydroxyl radicals generated.
    m = R.build_model()
    m.reactions.get_by_id("EX_glc").lower_bound = -5
    sol = m.optimize()
    assert sol.fluxes["DSUP"] <= R.DSUP_MAX_FRACTION * 0.24 * sol.fluxes["RADIO"] + 1e-6


def test_energy_budget_negligible_at_survivable_dose():
    # At ISS dose rate radiation supplies a vanishing fraction of ATP demand.
    frac = rc.energy_fraction(rc.DOSE_REGIMES["iss_leo"], rc.ATP_DEMAND_HIGH)
    assert frac < 1e-6


def test_energy_constrained_radio_flux_tiny_at_iss():
    cap = rc.max_radio_flux(rc.DOSE_REGIMES["iss_leo"], efficiency=1.0)
    assert cap < 1e-6, "energy-limited RADIO flux should be negligible at ISS dose"


def test_original_flux_requires_lethal_dose():
    # The original headline flux (~28.6) implies a promptly lethal dose rate.
    dose = rc.dose_for_radio_flux(28.57, efficiency=1.0)
    assert dose > 100, "original flux should require >100 Gy/s"


# --- 3. Kinetic model discriminates --------------------------------------
def test_kinetic_no_sod_raises_superoxide():
    base = K.run_simulation(K.STRESS_DOSE, duration=60)
    nosod = K.run_simulation(K.STRESS_DOSE, duration=60, use_sod=False)
    assert nosod.superoxide_M.iloc[-1] > 10 * base.superoxide_M.iloc[-1]


def test_kinetic_no_dsup_raises_dna_damage():
    base = K.run_simulation(K.STRESS_DOSE, duration=120)
    nodsup = K.run_simulation(K.STRESS_DOSE, duration=120, use_dsup=False)
    ratio = nodsup.dna_damage.iloc[-1] / base.dna_damage.iloc[-1]
    assert 1.5 < ratio < 1.8, f"expected ~1/0.6=1.67x, got {ratio:.2f}"


def test_kinetic_not_degenerate():
    # The former bug: identical output across doses. Ensure H2O2 responds to dose.
    lo = K.run_simulation(1.0, duration=60).h2o2_M.iloc[-1]
    hi = K.run_simulation(100.0, duration=60).h2o2_M.iloc[-1]
    assert hi > 5 * lo


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL  {t.__name__}: {e}")
        except Exception as e:
            print(f"  ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{passed}/{len(tests)} tests passed")
    raise SystemExit(0 if passed == len(tests) else 1)

"""
Invariant tests for the radiotrophic cell model.
================================================
Guards the core scientific claims against silent regressions. Run with:
    python3 test_model.py        (plain asserts, no pytest needed)
    pytest test_model.py         (if pytest is installed)

These are not unit tests of arithmetic; they assert that the *findings*
still hold: the model builds, the baseline is stable, Dsup protection is
capped at 40%, glutathione scavenging is finite, radiation is a usable
resource at low dose, and the defenses can be overwhelmed (lethality).
"""

from radiotrophic_model import (
    build_model, disable_engineered,
    DSUP_PROTECTION_FRACTION, OH_PER_RADIO, GSH_SCAV_CAP,
)

GLUCOSE = 5.0
TOL = 1e-6


def _solve_forced(model, dose, glucose=GLUCOSE):
    with model:
        model.reactions.get_by_id('EX_glc').lower_bound = -glucose
        model.reactions.get_by_id('RADIO').bounds = (dose, dose)
        return model.optimize()


def test_model_builds():
    m = build_model()
    assert len(m.reactions) >= 54
    assert len(m.metabolites) == 52
    assert any(c.name == 'dsup_protection_cap' for c in m.constraints)


def test_baseline_atp_stable():
    """Normal cell (no engineered pathways) at glucose=5 makes ~111.4 ATP."""
    m = build_model()
    with m:
        m.reactions.get_by_id('EX_glc').lower_bound = -GLUCOSE
        disable_engineered(m)
        sol = m.optimize()
    assert sol.status == 'optimal'
    assert abs(sol.objective_value - 111.39) < 0.5, sol.objective_value


def test_dsup_capped_at_40_percent():
    """Dsup may intercept at most 40% of generated OH radicals."""
    m = build_model()
    sol = _solve_forced(m, 20)
    assert sol.status == 'optimal'
    oh_generated = 20 * OH_PER_RADIO
    assert sol.fluxes['DSUP'] <= DSUP_PROTECTION_FRACTION * oh_generated + TOL


def test_gsh_scavenging_is_finite():
    """OH scavenging never exceeds the finite glutathione pool cap."""
    m = build_model()
    sol = _solve_forced(m, 28)
    assert sol.status == 'optimal'
    assert sol.fluxes['OH_SCAV'] <= GSH_SCAV_CAP + TOL


def test_low_dose_is_a_resource():
    """At low forced dose, radiation yields net ATP with zero DNA lesions."""
    m = build_model()
    base = _solve_forced(m, 0).objective_value
    sol = _solve_forced(m, 10)
    assert sol.status == 'optimal'
    assert sol.objective_value > base + TOL          # net energy gain
    assert sol.fluxes['FENTON'] < TOL                # no DNA damage


def test_high_dose_strains_then_kills():
    """Mid dose produces repairable lesions; high dose has no viable state."""
    m = build_model()
    strained = _solve_forced(m, 28)
    assert strained.status == 'optimal'
    assert strained.fluxes['FENTON'] > TOL           # lesions occur
    lethal = _solve_forced(m, 35)
    assert lethal.status != 'optimal'                # defenses overwhelmed


def test_sod_is_single_point_of_failure():
    """Knocking out SOD collapses radiotrophic flux (known bottleneck)."""
    m = build_model()
    with m:
        m.reactions.get_by_id('EX_glc').lower_bound = -GLUCOSE
        m.reactions.get_by_id('SODc').upper_bound = 0
        sol = m.optimize()
    assert sol.status == 'optimal'
    assert sol.fluxes['RADIO'] < TOL


def test_deterministic():
    """Same model solved twice gives identical objective (reproducibility)."""
    a = _solve_forced(build_model(), 15).objective_value
    b = _solve_forced(build_model(), 15).objective_value
    assert a == b


if __name__ == '__main__':
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    passed = 0
    for t in tests:
        try:
            t()
            print(f'PASS  {t.__name__}')
            passed += 1
        except AssertionError as e:
            print(f'FAIL  {t.__name__}: {e}')
        except Exception as e:
            print(f'ERROR {t.__name__}: {type(e).__name__}: {e}')
    print(f'\n{passed}/{len(tests)} tests passed')
    raise SystemExit(0 if passed == len(tests) else 1)

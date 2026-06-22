"""
Radiation dose calibration (illustrative / anchored).
=====================================================
The constraint-based model's fluxes are RELATIVE, not physically
dimensioned. A first-principles check confirms this: mapping the model's
hydroxyl-radical production (0.24 mmol per unit RADIO flux) to absorbed
energy via radiolysis G-values (Buxton et al. 1988, G(.OH) ~ 0.28 umol/J)
implies dose rates of ~1e5-1e6 Gy/h at the lethal threshold -- six orders
of magnitude above any real lethal dose (~single-digit Gy). The ROS
coefficients were chosen for their *ratio* (superoxide:OH), not absolute
scaling, so the radiation axis cannot be read as Gy from first principles.

We therefore provide an ANCHORED, ILLUSTRATIVE conversion: the model's
lethal threshold (the dose at which no feasible steady state exists) is
matched to a representative acute lethal dose for mammalian cells. This is
a calibration *target* for future work, NOT a validated physical mapping.
Any absolute Gy values below are order-of-magnitude illustrations only.

Reproducible: uses the pinned-solver model in radiotrophic_model.py.
"""

import numpy as np
from radiotrophic_model import build_model

# --- Anchor assumptions (explicit and adjustable) ---------------------
# Representative acute lethal dose for mammalian cells. Clonogenic survival
# falls steeply over ~2-10 Gy; the melanin-NP mouse study (Nat. Commun.
# 2025) used 6 Gy. We anchor the model's lethal threshold to 10 Gy as a
# round, defensible acute-lethal reference.
REFERENCE_LETHAL_GY = 10.0
GLUCOSE = 5.0


def find_lethal_flux(step=0.5, hi=60):
    """Lowest forced RADIO flux with no feasible steady state (= lethal)."""
    model = build_model()
    for dose in np.arange(0, hi + step, step):
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -GLUCOSE
            model.reactions.get_by_id('RADIO').bounds = (dose, dose)
            if model.optimize().status != 'optimal':
                return float(dose)
    return None


# Conversion factor: Gy per unit model flux, anchored at the lethal point.
_LETHAL_FLUX = find_lethal_flux()
FLUX_TO_GRAY = REFERENCE_LETHAL_GY / _LETHAL_FLUX if _LETHAL_FLUX else None


def flux_to_gray(flux):
    """Illustrative absolute dose (Gy) for a model flux. Anchored, approximate."""
    if FLUX_TO_GRAY is None:
        return float('nan')
    return flux * FLUX_TO_GRAY


if __name__ == '__main__':
    print('Radiation calibration (ILLUSTRATIVE / ANCHORED -- not validated)')
    print('-' * 64)
    print(f'Model lethal threshold flux : {_LETHAL_FLUX:g}')
    print(f'Anchored to                 : {REFERENCE_LETHAL_GY:g} Gy (acute, mammalian)')
    print(f'=> illustrative factor      : {FLUX_TO_GRAY:.3f} Gy per flux unit')
    print()
    for f in [21, 25, 29]:
        print(f'  flux {f:>3}  ~  {flux_to_gray(f):>4.1f} Gy (illustrative)')
    print()
    print('NOTE: first-principles radiolysis calibration overshoots by ~1e6,')
    print('confirming fluxes are relative. Absolute Gy values are illustrative.')

"""
Plot radiation dose vs. ATP production and lethality.
========================================================
Sweeps a FORCED radiation dose (the cell cannot opt out of a radiation
field) and records ATP output until the ROS defenses are overwhelmed and
no feasible steady state exists (lethality). Reproducible: GLPK solver.

Output: radiation_atp_lethality.png
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from radiotrophic_model import build_model, OH_PER_RADIO
from calibration import flux_to_gray, REFERENCE_LETHAL_GY

GLUCOSE = 5.0

model = build_model()

# Baseline ATP with no radiation (forced dose 0).
with model:
    model.reactions.get_by_id('EX_glc').lower_bound = -GLUCOSE
    model.reactions.get_by_id('RADIO').bounds = (0, 0)
    atp_baseline = model.optimize().objective_value

doses, atps, lesions = [], [], []
lethal_dose = None

for dose in np.arange(0, 40.01, 0.5):
    with model:
        model.reactions.get_by_id('EX_glc').lower_bound = -GLUCOSE
        model.reactions.get_by_id('RADIO').bounds = (dose, dose)
        s = model.optimize()
        doses.append(dose)
        if s.status == 'optimal':
            atps.append(s.objective_value)
            lesions.append(s.fluxes['FENTON'])
        else:
            atps.append(np.nan)        # no viable state -> cell dies
            lesions.append(np.nan)
            if lethal_dose is None:
                lethal_dose = dose

doses = np.array(doses)
atps = np.array(atps)
lesions = np.array(lesions)

# First dose at which DNA lesions begin (transition RESOURCE -> STRAINED).
strain_idx = np.where(np.nan_to_num(lesions) > 1e-6)[0]
strain_dose = doses[strain_idx[0]] if len(strain_idx) else None

fig, ax1 = plt.subplots(figsize=(10, 6))

# --- ATP production curve ---
ax1.plot(doses, atps, color='#1f77b4', lw=2.5, marker='o', ms=3,
         label='ATP production', zorder=5)
ax1.axhline(atp_baseline, color='#1f77b4', ls=':', lw=1, alpha=0.6)
ax1.annotate('no-radiation baseline', xy=(1, atp_baseline),
             xytext=(1, atp_baseline + 1.2), fontsize=8, color='#1f77b4')
ax1.set_xlabel('Radiation dose (forced radiotrophic flux)', fontsize=11)
ax1.set_ylabel('ATP production (flux units)', color='#1f77b4', fontsize=11)
ax1.tick_params(axis='y', labelcolor='#1f77b4')

# --- DNA lesions on secondary axis ---
ax2 = ax1.twinx()
ax2.plot(doses, lesions, color='#d62728', lw=2, ls='--', marker='s', ms=3,
         label='DNA lesions (unrepaired-rate)', zorder=4)
ax2.set_ylabel('DNA lesions requiring repair', color='#d62728', fontsize=11)
ax2.tick_params(axis='y', labelcolor='#d62728')
ax2.set_ylim(bottom=0)

# --- Regime shading ---
xmax = doses.max()
if strain_dose is not None:
    ax1.axvspan(0, strain_dose, color='#2ca02c', alpha=0.08)
    ax1.axvspan(strain_dose, lethal_dose if lethal_dose else xmax,
                color='#ff7f0e', alpha=0.10)
if lethal_dose is not None:
    ax1.axvspan(lethal_dose, xmax, color='#d62728', alpha=0.12)
    ax1.axvline(lethal_dose, color='#d62728', lw=2, ls='-')
    ax1.annotate(f'LETHAL threshold\n(dose = {lethal_dose:g})',
                 xy=(lethal_dose, atp_baseline),
                 xytext=(lethal_dose + 0.5, atp_baseline - 4),
                 fontsize=9, color='#d62728', fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color='#d62728'))

legend_elems = [
    Patch(facecolor='#2ca02c', alpha=0.25, label='RESOURCE (radicals fully neutralized)'),
    Patch(facecolor='#ff7f0e', alpha=0.30, label='STRAINED (net gain, DNA lesions occur)'),
    Patch(facecolor='#d62728', alpha=0.30, label='LETHAL (defenses overwhelmed)'),
]
lines1, labs1 = ax1.get_legend_handles_labels()
lines2, labs2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2 + legend_elems,
           labs1 + labs2 + [e.get_label() for e in legend_elems],
           loc='lower left', fontsize=8, framealpha=0.9)

# --- Illustrative absolute-dose axis (anchored, see calibration.py) ---
ax_top = ax1.secondary_xaxis('top', functions=(flux_to_gray,
                                               lambda g: g / flux_to_gray(1)))
ax_top.set_xlabel('Illustrative absolute dose (Gy, anchored -- approximate)',
                  fontsize=9, color='dimgray')
ax_top.tick_params(colors='dimgray', labelsize=8)

ax1.set_title('Radiotrophic human cell: ATP production and lethality vs. radiation dose',
              fontsize=12, fontweight='bold', pad=34)
fig.tight_layout()
fig.savefig('radiation_atp_lethality.png', dpi=150)
print(f'baseline ATP (no radiation): {atp_baseline:.2f}')
print(f'STRAINED onset (first DNA lesions) at dose: {strain_dose}')
print(f'LETHAL threshold at dose: {lethal_dose}')
print('Saved radiation_atp_lethality.png')

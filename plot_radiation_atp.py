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

from radiotrophic_model import build_model
from calibration import flux_to_gray

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

xmax = doses.max()
# Peak ATP (the RESOURCE ceiling) among the feasible doses.
feasible = ~np.isnan(atps)
peak_i = np.nanargmax(atps)
peak_dose, peak_atp = doses[peak_i], atps[peak_i]

fig, ax1 = plt.subplots(figsize=(12, 7))

# --- Regime shading (drawn first, behind the curves) ---
resource_end = strain_dose if strain_dose is not None else (
    lethal_dose if lethal_dose is not None else xmax)
ax1.axvspan(0, resource_end, color='#2ca02c', alpha=0.10, zorder=0)
if strain_dose is not None:
    ax1.axvspan(strain_dose, lethal_dose if lethal_dose else xmax,
                color='#ff7f0e', alpha=0.13, zorder=0)
if lethal_dose is not None:
    ax1.axvspan(lethal_dose, xmax, color='#d62728', alpha=0.13, zorder=0)

# --- ATP production curve (left axis) ---
ax1.plot(doses, atps, color='#1f77b4', lw=2.6, marker='o', ms=4,
         label='ATP production', zorder=5)
ax1.axhline(atp_baseline, color='#1f77b4', ls=':', lw=1.2, alpha=0.7, zorder=1)
ax1.set_xlabel('Radiation dose  (forced radiotrophic flux, model units)', fontsize=12)
ax1.set_ylabel('ATP production (flux units)', color='#1f77b4', fontsize=12)
ax1.tick_params(axis='y', labelcolor='#1f77b4')
ax1.set_xlim(-1, xmax + 1)
# Headroom so annotations sit above the curve without touching it.
ax1.set_ylim(atp_baseline - 1.2, peak_atp + 2.2)

# --- DNA lesions (right axis) ---
ax2 = ax1.twinx()
ax2.plot(doses, lesions, color='#d62728', lw=2.2, ls='--', marker='s', ms=4,
         label='DNA lesions (repair load)', zorder=4)
ax2.set_ylabel('DNA lesions requiring repair (flux units)', color='#d62728', fontsize=12)
ax2.tick_params(axis='y', labelcolor='#d62728')
ax2.set_ylim(0, float(np.nanmax(lesions)) * 1.5 + 1e-6)

# --- In-plot annotations (no overlaps: each sits in its own band) ---
# baseline label, parked at the far left just above the dotted line
ax1.text(0.2, atp_baseline + 0.12, 'no-radiation baseline', fontsize=9,
         color='#1f77b4', va='bottom')
# RESOURCE ceiling = ATP peak; if lesions start at the same dose, say so once.
ceiling_note = f'RESOURCE ceiling (dose {peak_dose:g})\npeak ATP'
if strain_dose is not None and abs(strain_dose - peak_dose) < 1.0:
    ceiling_note += ' — DNA lesions begin here'
ax1.annotate(ceiling_note,
             xy=(peak_dose, peak_atp), xytext=(peak_dose - 11, peak_atp + 1.1),
             fontsize=9.5, color='#2e7d32', fontweight='bold', ha='left',
             arrowprops=dict(arrowstyle='->', color='#2e7d32'))
# If lesions begin well before the peak, call that out separately.
if strain_dose is not None and abs(strain_dose - peak_dose) >= 1.0:
    ax1.annotate('first DNA lesions\n(STRAINED onset)',
                 xy=(strain_dose, atp_baseline + 0.1),
                 xytext=(strain_dose - 1.5, peak_atp - 0.3),
                 fontsize=9, color='#b9770e', ha='center',
                 arrowprops=dict(arrowstyle='->', color='#b9770e'))
# lethal threshold: vertical line + boxed label INSIDE the plot
if lethal_dose is not None:
    ax1.axvline(lethal_dose, color='#d62728', lw=2.2, zorder=3)
    ax1.text(lethal_dose + 0.6, peak_atp + 0.9,
             f'LETHAL threshold\ndose = {lethal_dose:g}', fontsize=10,
             color='#d62728', fontweight='bold', va='top', ha='left',
             bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#d62728', lw=1.2))

# --- Illustrative absolute-dose axis on top ---
ax_top = ax1.secondary_xaxis('top', functions=(flux_to_gray,
                                               lambda g: g / flux_to_gray(1)))
ax_top.set_xlabel('Illustrative absolute dose  (Gy, anchored — approximate, not calibrated)',
                  fontsize=10, color='dimgray')
ax_top.tick_params(colors='dimgray', labelsize=9)

# --- Combined legend, OUTSIDE the axes at the bottom (no data overlap) ---
legend_elems = [
    Patch(facecolor='#2ca02c', alpha=0.30, label='RESOURCE — radicals fully neutralized'),
    Patch(facecolor='#ff7f0e', alpha=0.35, label='STRAINED — net gain, lesions repaired at ATP cost'),
    Patch(facecolor='#d62728', alpha=0.35, label='LETHAL — defenses overwhelmed (cell dies)'),
]
lines1, labs1 = ax1.get_legend_handles_labels()
lines2, labs2 = ax2.get_legend_handles_labels()
fig.legend(lines1 + lines2 + legend_elems, labs1 + labs2 + [e.get_label() for e in legend_elems],
           loc='lower center', ncol=3, fontsize=9, framealpha=0.95,
           bbox_to_anchor=(0.5, -0.02))

ax1.set_title('Radiotrophic human cell: ATP production and lethality vs. radiation dose',
              fontsize=13, fontweight='bold', pad=30)
fig.tight_layout(rect=(0, 0.08, 1, 1))
fig.savefig('radiation_atp_lethality.png', dpi=150, bbox_inches='tight')
print(f'baseline ATP (no radiation): {atp_baseline:.2f}')
print(f'STRAINED onset (first DNA lesions) at dose: {strain_dose}')
print(f'LETHAL threshold at dose: {lethal_dose}')
print('Saved radiation_atp_lethality.png')

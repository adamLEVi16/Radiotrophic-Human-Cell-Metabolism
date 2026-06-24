"""
Mechanistic schematic: radiation impacting an engineered radiotrophic cell.
===========================================================================
Conceptual figure (not data) showing how the melanin + engineered-defense
"treatment" turns incoming radiation into energy while neutralising the
radical byproducts to protect DNA.

Output: fig_mechanism.png
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import (FancyBboxPatch, Circle, Ellipse, FancyArrowPatch,
                                FancyArrow, Patch, Wedge)
import numpy as np

# Colour scheme
ENERGY = '#2ca02c'      # green - energy/good outcome
DAMAGE = '#d62728'      # red   - ROS/damage
ENGINEER = '#8e44ad'    # purple- engineered / cross-species parts
NATIVE = '#1f77b4'      # blue  - native human defenses
RAD = '#f1c40f'         # yellow- radiation
MELANIN = '#2c2c2c'     # near-black melanin

fig, ax = plt.subplots(figsize=(14, 9.6))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')


def box(x, y, w, h, text, fc, ec='black', fs=10, tc='black', style='round', lw=1.5):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f'{style},pad=0.3',
                       fc=fc, ec=ec, lw=lw, zorder=3)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center',
            fontsize=fs, color=tc, zorder=4, wrap=True)


def arrow(x1, y1, x2, y2, color='black', lw=2, style='-|>', ls='-', mut=18):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                 mutation_scale=mut, color=color, lw=lw, linestyle=ls, zorder=2))


# ---- Cell body ----
cell = FancyBboxPatch((8, 8), 84, 76, boxstyle='round,pad=0.6',
                      fc='#eaf6ff', ec=NATIVE, lw=2.5, zorder=0)
ax.add_patch(cell)
ax.text(50, 97, 'Engineered radiotrophic human cell', ha='center',
        fontsize=15, fontweight='bold')

# ---- Incoming radiation (top-left), jagged arrows toward melanin ----
ax.text(13, 84, 'Ionizing radiation\n(γ / X-ray)', ha='center', fontsize=11,
        color='#b8860b', fontweight='bold')
for dy in (-2, 0, 2):
    arrow(4, 80 + dy, 26, 70 + dy, color=RAD, lw=3, mut=22)

# ---- Melanin granules (absorber) ----
for (mx, my) in [(30, 70), (34, 73), (33, 67), (37, 70), (29, 74)]:
    ax.add_patch(Circle((mx, my), 2.4, fc=MELANIN, ec='none', zorder=3))
ax.text(33, 60.5, 'Melanin\n(radiation absorber)', ha='center', fontsize=10,
        fontweight='bold')

# ---- Radiotrophic energy path: melanin -> NADH -> mitochondrion -> ATP ----
box(46, 64, 20, 8, 'Radiotrophic\nNADH generation', '#d5f5e3', ENERGY, 10, lw=2)
arrow(40, 69, 46, 68, color=ENERGY, lw=2.5)

# Mitochondrion (label sits inside on a white pad so cristae don't clash)
mito = Ellipse((79, 60), 22, 14, fc='#fdebd0', ec='#e67e22', lw=2, zorder=2)
ax.add_patch(mito)
for i in range(4):  # cristae
    ax.add_patch(Ellipse((73 + i * 3.3, 60), 2.0, 7.5, fc='none', ec='#e67e22',
                 lw=1, zorder=2))
ax.text(79, 60, 'Mitochondrion\n(ETC → ATP)', ha='center', va='center',
        fontsize=9, color='#a04000', fontweight='bold', zorder=5,
        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.78))
arrow(66, 68, 73, 63, color=ENERGY, lw=2.5)
box(85, 69, 11, 7, 'ATP ⚡', '#abebc6', ENERGY, 12, lw=2)
arrow(86, 60, 90, 69, color=ENERGY, lw=2.5)

# ---- ROS byproducts (red) branching down from the radiotrophic step ----
ax.text(42, 58, 'byproduct: reactive oxygen species', ha='center', fontsize=9,
        color=DAMAGE, style='italic')
box(20, 46, 16, 7, 'Superoxide\nO₂•⁻', '#fadbd8', DAMAGE, 10, lw=1.5)
box(58, 46, 16, 7, 'Hydroxyl radical\n•OH', '#fadbd8', DAMAGE, 10, lw=1.5)
arrow(50, 64, 30, 53, color=DAMAGE, lw=2, ls=':')
arrow(54, 64, 64, 53, color=DAMAGE, lw=2, ls=':')

# ---- Superoxide defense chain: SOD -> H2O2 -> catalase/GPX -> water ----
box(16, 33, 13, 7, 'SOD\n+ MnSOD2*', '#e8daef', ENGINEER, 9, lw=1.5)
arrow(27, 46, 24, 40, color=DAMAGE, lw=2)
box(33, 33, 11, 7, 'H₂O₂', '#fdebd0', '#e67e22', 9, lw=1.5)
arrow(29, 36.5, 33, 36.5, color='black', lw=1.8)
box(15, 22, 15, 7, 'Catalase / GPX', '#d6eaf8', NATIVE, 9, lw=1.5)
arrow(36, 33, 26, 29, color='black', lw=1.8)
box(33, 22, 12, 7, 'H₂O\n(safe)', '#d5f5e3', ENERGY, 9.5, lw=1.5)
arrow(30, 25.5, 33, 25.5, color=ENERGY, lw=2)

# Failure-point callout: superoxide is the species that overwhelms SOD first.
ax.annotate('⚠ FAILS FIRST\nO₂•⁻ overwhelms SOD\npast the dose ceiling',
            xy=(17, 33.5), xytext=(2.5, 13.5), fontsize=8.5, color=DAMAGE,
            fontweight='bold', ha='left', va='center', zorder=6,
            arrowprops=dict(arrowstyle='->', color=DAMAGE, lw=1.6))

# ---- OH defense: GSH scavenging + Dsup shield -> protect DNA ----
box(49, 33, 16, 7, 'GSH scavenging\n(finite pool)', '#d6eaf8', NATIVE, 9, lw=1.5)
arrow(64, 46, 60, 40, color=DAMAGE, lw=2)

# Nucleus with Dsup shield
ax.add_patch(Wedge((78, 21), 10, 0, 360, width=2.0, fc='none', ec=ENGINEER,
             lw=2.5, ls='--', zorder=3))
ax.add_patch(Circle((78, 21), 8, fc='#d6eaf8', ec=NATIVE, lw=1.5, zorder=3))
t = np.linspace(0, 4 * np.pi, 100)
ax.plot(78 + 2.0 * np.sin(t), 16.5 + t * 0.5, color=NATIVE, lw=1.2, zorder=4)
ax.plot(78 - 2.0 * np.sin(t), 16.5 + t * 0.5, color=NATIVE, lw=1.2, zorder=4)
ax.text(78, 21, 'DNA', ha='center', va='center', fontsize=9, zorder=5,
        fontweight='bold', color='#154360')
ax.text(78, 32.5, 'Dsup* DNA shield', ha='center', fontsize=9,
        color=ENGINEER, fontweight='bold')
arrow(60, 33, 71, 26, color=DAMAGE, lw=2, ls=':')      # residual OH toward nucleus
arrow(75, 29, 76, 26, color=ENGINEER, lw=2)             # Dsup intercept

# ---- Outcome banner (bottom-left, clear of the nucleus on the right) ----
box(20, 9.5, 42, 5.5,
    'OUTCOME (up to the dose ceiling):   ⚡ ATP gained    ✓ DNA protected',
    '#eafaf1', ENERGY, 10, ENERGY, lw=2)

# ---- Legend ----
legend = [
    Patch(fc=RAD, label='ionizing radiation'),
    Patch(fc=MELANIN, label='melanin (absorber)'),
    Patch(fc='#d5f5e3', ec=ENERGY, label='energy pathway / safe products'),
    Patch(fc='#fadbd8', ec=DAMAGE, label='reactive oxygen species (damage)'),
    Patch(fc='#d6eaf8', ec=NATIVE, label='native human defenses'),
    Patch(fc='#e8daef', ec=ENGINEER, label='engineered / cross-species (*)'),
]
ax.legend(handles=legend, loc='center', fontsize=9.5, framealpha=0.95,
          ncol=3, bbox_to_anchor=(0.5, 0.905))

ax.text(50, 4.5,
        '* engineered components: tardigrade Dsup, overexpressed MnSOD2, '
        'Deinococcus Mn-antioxidant, enhanced Nrf2.  Beyond the dose ceiling, '
        'ROS overwhelms the defenses → cell death.',
        ha='center', fontsize=8.5, color='#555555', style='italic')

fig.tight_layout()
fig.savefig('fig_mechanism.png', dpi=150, bbox_inches='tight')
print('wrote fig_mechanism.png')

"""
Generate the supporting paper figures from experiment CSVs.
===========================================================
  Fig 2  fig_ablation.png          - defense knockouts (single points of failure)
  Fig 3  fig_bottleneck_relief.png - lifting the radiotrophic ceiling
  Fig 4  fig_ros_sensitivity.png   - robustness to the ROS:NADH coefficient

Run after radiotrophic_model.py. Reproducible (reads committed CSVs).
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.patches import Patch

BLUE, RED, GREEN, GRAY = '#1f77b4', '#d62728', '#2ca02c', '#888888'
AMBER = '#e69500'


# ----------------------------------------------------------------------
# Fig 2: Defense ablation
# ----------------------------------------------------------------------
def fig_ablation():
    df = pd.read_csv('ablation.csv')
    full = df['radio_flux'].max()
    fig, ax = plt.subplots(figsize=(11.5, 6.5))
    y = np.arange(len(df))

    # Three tiers: dead (0), impaired (reduced), intact (full ceiling).
    def tier(rf):
        if rf < 1e-6:
            return RED
        if rf < full - 1e-6:
            return AMBER
        return BLUE
    colors = [tier(rf) for rf in df['radio_flux']]
    ax.barh(y, df['radio_flux'], color=colors, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(df['config'], fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(0, full * 1.30)
    ax.set_xlabel('Radiotrophic flux achieved  (energy capture, model units)', fontsize=11)
    ax.set_title('Defense ablation — which knockouts kill radiotrophy?',
                 fontsize=13, fontweight='bold')
    # Reference line at the intact ceiling.
    ax.axvline(full, color=BLUE, ls=':', lw=1.2, alpha=0.7, zorder=1)
    ax.text(full, -0.7, f'intact ceiling = {full:g}', color=BLUE, fontsize=8,
            ha='center', va='bottom')
    # Per-bar labels: flux + the resulting ATP, no overlaps (text past bar end).
    for yi, rf, atp, c in zip(y, df['radio_flux'], df['atp'], colors):
        if c == RED:
            ax.text(0.3, yi, f'✗ collapses to 0   (ATP {atp:g})', va='center',
                    fontsize=9, color=RED, fontweight='bold')
        else:
            tag = '  (full)' if c == BLUE else '  (impaired)'
            ax.text(rf + 0.3, yi, f'{rf:g}{tag}   ATP {atp:g}', va='center',
                    fontsize=8.5, color='black')
    ax.axvline(0, color='k', lw=0.8)
    ax.legend(handles=[
        Patch(color=BLUE, label='intact — full radiotrophic ceiling'),
        Patch(color=AMBER, label='impaired — capture reduced'),
        Patch(color=RED, label='collapsed — single point of failure'),
    ], loc='lower right', fontsize=9, framealpha=0.95)
    fig.tight_layout()
    fig.savefig('fig_ablation.png', dpi=150)
    plt.close(fig)
    print('wrote fig_ablation.png')


# ----------------------------------------------------------------------
# Fig 3: Bottleneck relief
# ----------------------------------------------------------------------
def fig_bottleneck_relief():
    df = pd.read_csv('bottleneck_relief.csv')
    x = np.arange(len(df))
    fig, ax1 = plt.subplots(figsize=(11.5, 6.8))

    bars = ax1.bar(x - 0.2, df['radio_flux'], width=0.4, color=BLUE, zorder=3)
    ax1.set_ylabel('Radiotrophic ceiling — flux achieved', color=BLUE, fontsize=11)
    ax1.tick_params(axis='y', labelcolor=BLUE)
    ax1.set_ylim(0, df['radio_flux'].max() * 1.22)   # headroom for labels + note
    for xi, v in zip(x, df['radio_flux']):
        ax1.text(xi - 0.2, v + 0.8, f'{v:g}', ha='center', fontsize=9,
                 color=BLUE, fontweight='bold')

    ax2 = ax1.twinx()
    atp_bars = ax2.bar(x + 0.2, df['atp'], width=0.4, color=GREEN, alpha=0.85, zorder=3)
    ax2.set_ylabel('ATP production', color=GREEN, fontsize=11)
    ax2.tick_params(axis='y', labelcolor=GREEN)
    atp_lo = df['atp'].min() - 4
    ax2.set_ylim(atp_lo, df['atp'].max() + (df['atp'].max() - atp_lo) * 0.22)
    for xi, v in zip(x, df['atp']):
        ax2.text(xi + 0.2, v + 0.4, f'{v:g}', ha='center', fontsize=8.5,
                 color='#1e7d32')

    ax1.set_xticks(x)
    ax1.set_xticklabels(df['strategy'], rotation=12, ha='right', fontsize=9.5)
    ax1.set_title('Bottleneck relief — lifting the radiotrophic ceiling',
                  fontsize=13, fontweight='bold')
    # Explanatory note parked in clear space over the short middle bars.
    ax1.text(0.5, 0.86,
             'MnSOD2 alone (C) does nothing — SOD is not the binding limit.\n'
             'Relieving OH-neutralisation (GSH, B/D) is what lifts the ceiling.',
             transform=ax1.transAxes, ha='center', va='top', fontsize=9,
             color=GRAY, style='italic',
             bbox=dict(boxstyle='round,pad=0.4', fc='white', ec=GRAY, alpha=0.85))
    ax1.legend([bars, atp_bars], ['radiotrophic ceiling (flux)', 'ATP production'],
               loc='upper left', fontsize=9, framealpha=0.95)
    fig.tight_layout()
    fig.savefig('fig_bottleneck_relief.png', dpi=150)
    plt.close(fig)
    print('wrote fig_bottleneck_relief.png')


# ----------------------------------------------------------------------
# Fig 4: ROS-coefficient sensitivity
# ----------------------------------------------------------------------
def fig_ros_sensitivity():
    df = pd.read_csv('ros_sensitivity.csv')
    fig, ax1 = plt.subplots(figsize=(11.5, 6.8))

    ax1.plot(df['ros_coefficient'], df['net_gain'], color=GREEN, lw=2.6,
             marker='o', ms=6, label='net ATP gain from radiation', zorder=5)
    ax1.plot(df['ros_coefficient'], df['radio_flux'], color=BLUE, lw=2.2,
             marker='s', ms=5, ls='--', label='radiotrophic flux', zorder=4)
    ax1.axhline(0, color='k', lw=0.8, alpha=0.5)
    ax1.set_xlabel('ROS production coefficient  (superoxide per NADH) — the most uncertain parameter',
                   fontsize=10.5)
    ax1.set_ylabel('Net ATP gain  /  radiotrophic flux', fontsize=11)
    ax1.set_ylim(0, df['net_gain'].max() * 1.15)

    # Literature range band + label (kept clear of the curves, low on the axis).
    ax1.axvspan(0.28, 0.32, color=GRAY, alpha=0.18, zorder=0)
    ax1.text(0.30, df['net_gain'].max() * 0.30, 'literature\nrange',
             ha='center', va='center', fontsize=8.5, color='#555555')
    # "always positive" takeaway, parked on the right where the curves are low.
    ax1.text(0.97, 0.30, 'net gain stays > 0\nacross the whole range',
             transform=ax1.transAxes, ha='right', va='center', fontsize=9,
             color=GREEN, style='italic',
             bbox=dict(boxstyle='round,pad=0.35', fc='white', ec=GREEN, alpha=0.9))

    ax2 = ax1.twinx()
    ax2.plot(df['ros_coefficient'], df['dna_repair'], color=RED, lw=2,
             marker='^', ms=6, ls=':', label='DNA lesions (repair load)', zorder=4)
    ax2.set_ylabel('DNA lesions requiring repair', color=RED, fontsize=11)
    ax2.tick_params(axis='y', labelcolor=RED)
    ax2.set_ylim(0, df['dna_repair'].max() * 1.5 + 1e-6)
    # Mark where lesions appear (low ROS cost -> high flux -> damage).
    ax2.annotate('lesions only at low ROS cost\n(high flux drives damage)',
                 xy=(0.2, df['dna_repair'].max() * 0.66),
                 xytext=(0.6, df['dna_repair'].max() * 0.9),
                 fontsize=8.5, color=RED,
                 arrowprops=dict(arrowstyle='->', color=RED))

    ax1.set_title('Robustness — net benefit persists across the ROS-cost range',
                  fontsize=13, fontweight='bold')
    l1, lab1 = ax1.get_legend_handles_labels()
    l2, lab2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lab1 + lab2, loc='upper center', fontsize=9.5,
               framealpha=0.95)
    fig.tight_layout()
    fig.savefig('fig_ros_sensitivity.png', dpi=150)
    plt.close(fig)
    print('wrote fig_ros_sensitivity.png')


if __name__ == '__main__':
    fig_ablation()
    fig_bottleneck_relief()
    fig_ros_sensitivity()
    print('All supporting figures generated.')

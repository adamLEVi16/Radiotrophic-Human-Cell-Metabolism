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

BLUE, RED, GREEN, GRAY = '#1f77b4', '#d62728', '#2ca02c', '#888888'


# ----------------------------------------------------------------------
# Fig 2: Defense ablation
# ----------------------------------------------------------------------
def fig_ablation():
    df = pd.read_csv('ablation.csv')
    fig, ax = plt.subplots(figsize=(10, 6))
    y = np.arange(len(df))
    collapsed = df['radio_flux'] < 1e-6
    colors = [RED if c else BLUE for c in collapsed]
    ax.barh(y, df['radio_flux'], color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(df['config'])
    ax.invert_yaxis()
    ax.set_xlim(0, df['radio_flux'].max() * 1.18)
    ax.set_xlabel('Radiotrophic flux achieved (energy capture)', fontsize=11)
    ax.set_title('Defense ablation: which knockouts kill radiotrophy?',
                 fontsize=12, fontweight='bold')
    # Label every bar; mark the collapses explicitly in red at x=0.
    for yi, rf, c in zip(y, df['radio_flux'], collapsed):
        if c:
            ax.text(0.3, yi, '✗ collapses to 0', va='center', fontsize=9,
                    color=RED, fontweight='bold')
        else:
            ax.text(rf + 0.3, yi, f'{rf:g}', va='center', fontsize=8)
    ax.axvline(0, color='k', lw=0.8)
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color=BLUE, label='radiotrophy viable'),
        Patch(color=RED, label='radiotrophy collapses (single point of failure)'),
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
    fig, ax1 = plt.subplots(figsize=(10, 6))

    bars = ax1.bar(x - 0.2, df['radio_flux'], width=0.4, color=BLUE,
                   label='radiotrophic ceiling (flux)')
    ax1.set_ylabel('Radiotrophic flux achieved', color=BLUE, fontsize=11)
    ax1.tick_params(axis='y', labelcolor=BLUE)
    for xi, v in zip(x, df['radio_flux']):
        ax1.text(xi - 0.2, v + 0.6, f'{v:g}', ha='center', fontsize=8, color=BLUE)

    ax2 = ax1.twinx()
    atp_bars = ax2.bar(x + 0.2, df['atp'], width=0.4, color=GREEN, alpha=0.8,
                       label='ATP production')
    ax2.set_ylabel('ATP production', color=GREEN, fontsize=11)
    ax2.tick_params(axis='y', labelcolor=GREEN)
    ax2.set_ylim(110, df['atp'].max() * 1.05)

    ax1.set_xticks(x)
    ax1.set_xticklabels(df['strategy'], rotation=15, ha='right', fontsize=9)
    ax1.set_title('Bottleneck relief: lifting the radiotrophic ceiling',
                  fontsize=12, fontweight='bold')
    ax1.text(0.5, 0.95,
             'MnSOD2 alone (C) does nothing — SOD is not the binding limit.\n'
             'Relieving OH-neutralisation (GSH) is what lifts the ceiling.',
             transform=ax1.transAxes, ha='center', va='top', fontsize=8,
             color=GRAY, style='italic')
    ax1.legend([bars, atp_bars], ['radiotrophic ceiling (flux)', 'ATP production'],
               loc='upper left', fontsize=8)
    fig.tight_layout()
    fig.savefig('fig_bottleneck_relief.png', dpi=150)
    plt.close(fig)
    print('wrote fig_bottleneck_relief.png')


# ----------------------------------------------------------------------
# Fig 4: ROS-coefficient sensitivity
# ----------------------------------------------------------------------
def fig_ros_sensitivity():
    df = pd.read_csv('ros_sensitivity.csv')
    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(df['ros_coefficient'], df['net_gain'], color=GREEN, lw=2.5,
             marker='o', label='net ATP gain from radiation')
    ax1.plot(df['ros_coefficient'], df['radio_flux'], color=BLUE, lw=2,
             marker='s', ls='--', label='radiotrophic flux')
    ax1.set_xlabel('ROS production coefficient (superoxide per NADH) — most uncertain parameter',
                   fontsize=10)
    ax1.set_ylabel('Net ATP gain / radiotrophic flux', fontsize=11)
    ax1.axvspan(0.28, 0.32, color=GRAY, alpha=0.15)
    ax1.text(0.30, ax1.get_ylim()[1] * 0.92, 'literature\nrange', ha='center',
             fontsize=8, color=GRAY)

    ax2 = ax1.twinx()
    ax2.plot(df['ros_coefficient'], df['dna_repair'], color=RED, lw=2,
             marker='^', ls=':', label='DNA repair (lesions)')
    ax2.set_ylabel('DNA repair flux (lesions)', color=RED, fontsize=11)
    ax2.tick_params(axis='y', labelcolor=RED)
    ax2.set_ylim(bottom=0)

    ax1.set_title('Robustness: net benefit persists across the ROS-cost range',
                  fontsize=12, fontweight='bold')
    l1, lab1 = ax1.get_legend_handles_labels()
    l2, lab2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lab1 + lab2, loc='upper right', fontsize=9)
    fig.tight_layout()
    fig.savefig('fig_ros_sensitivity.png', dpi=150)
    plt.close(fig)
    print('wrote fig_ros_sensitivity.png')


if __name__ == '__main__':
    fig_ablation()
    fig_bottleneck_relief()
    fig_ros_sensitivity()
    print('All supporting figures generated.')

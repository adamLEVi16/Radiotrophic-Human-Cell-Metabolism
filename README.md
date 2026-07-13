# Radiotrophic Human Cell Metabolism

A computational feasibility study of engineering **radiation-powered
("radiotrophic") metabolism** into human cells, using melanin for energy
capture plus three cross-species radioprotection systems: tardigrade **Dsup**,
*Deinococcus radiodurans* **Mn-antioxidant** complexes, and enhanced **Nrf2**
antioxidant signalling.

The project asks two questions and answers them with three complementary
models:

1. **Can radiation meaningfully power a human cell?** (energy budget)
2. **If a radiotrophic pathway ran, could the cell survive the ROS?** (metabolic
   + kinetic ROS models)

## Headline finding

**Radiotrophy cannot meaningfully power human cells at any survivable dose
rate.** By conservation of energy, absorbed ionizing radiation supplies between
**10⁻⁶ % and 10⁻¹⁰ %** of a cell's ATP demand at realistic exposures (natural
background, low-Earth-orbit/ISS, Chernobyl hotspots). Break-even — radiation
supplying as much power as metabolism — requires **~0.3–28 Gy/s**, a dose rate
that delivers a whole-body-lethal dose in seconds and would kill the cell long
before any energetic benefit accrued (`fig1_energy_budget.png`).

The radioprotection genes (Dsup, Mn-AOX, Nrf2, SOD2) remain independently
useful for radiation *tolerance* — but that is shielding and repair, not energy
metabolism. The two ideas should not be conflated.

> An earlier version of this model reported a "+18 % ATP boost." That number was
> an artifact of a radiotrophic reaction that generated reducing equivalents
> without any link to the energy actually deposited by radiation. Once energy
> conservation is enforced (`energy_constrained_dose`), the boost at every
> survivable dose rate is ~0 (`fig5_energy_constrained.png`).

## Repository layout

| File | Purpose |
|---|---|
| `radiotrophic_common.py` | Shared physical basis: water-radiolysis G-values, dose→ROS and dose→energy conversions, reference dose regimes. Single source of truth for both models. |
| `energy_budget.py` | Phase 3 — radiation power vs ATP demand across dose regimes → `energy_budget.csv` |
| `radiotrophic_model.py` | Constraint-based (FBA) metabolic model, COBRApy. 11 experiments incl. energy-constrained dose, FVA, Monte Carlo. |
| `kinetic_model.py` | Dose-driven ODE model of transient ROS dynamics (scipy). |
| `dsup_analysis.py` | Literature-synthesized structural feasibility assessment for Dsup. |
| `make_figures.py` | Renders all figures to `figures/`. |
| `test_model.py` | Test suite (structural, scientific, regression). |
| `run_all.py` | Reproduces every CSV, JSON, and figure. |

## Reproduce

```bash
pip install -r requirements.txt
python run_all.py      # regenerates all CSVs, the model JSON, and figures/
python test_model.py   # or: pytest test_model.py
```

## The three models

### 1. Energy budget (the decisive test)
Uses only conservation of energy and water-radiolysis yields. Absorbed power
per cell is `dose_rate × cell_mass`; ATP-demand power is
`turnover × ΔG_ATP`. The ratio is tiny at every survivable dose. This is the
result the metabolic model cannot produce on its own, because a lumped
"radiotrophic" reaction hides the energy accounting.

### 2. Constraint-based metabolic model (FBA)
A 56-reaction COBRApy model (glycolysis, TCA, ETC, melanin synthesis, four ROS
defense systems, DNA-damage/repair). Key improvements over the original:

- **Energy-constrained RADIO flux** (`energy_constrained_dose`): the
  radiotrophic reaction's upper bound is set by the radiation energy available
  at each dose rate, not an arbitrary cap. This is the physically correct
  version of the dose-response experiment.
- **Dsup capped at its measured 40 %** of hydroxyl-radical interception
  (Hashimoto et al. 2016), enforced as a flux-ratio constraint. Previously FBA
  routed 100 % of ·OH through Dsup, overstating its role.
- **Flux Variability Analysis** (`fva`): results reported as ranges at 95 % of
  optimal, not single arbitrary optima.
- **Monte Carlo** (`monte_carlo`): the abstract-model ATP boost is characterized
  over literature-justified parameter ranges (median ~11 %, 5th–95th percentile
  0–26 %) — its *sensitivity*, not an endorsement. The energy budget overrides
  it regardless.

### 3. Kinetic ROS model (ODE)
Transient intracellular ROS during exposure, **driven by absorbed dose rate**
via the shared G-values. The previous version was driven by an abstract flux
~6 orders of magnitude below scavenging capacity, so ROS never rose above its
initial condition and every dose and knockout produced identical output. The
rewrite accumulates ROS, saturates enzymes, and discriminates defenses
(run at a deliberately extreme FLASH-radiotherapy dose, 100 Gy/s, to expose
each system's role):

- Remove **SOD** → superoxide rises ~40× (`fig3`)
- Remove **catalase** → H₂O₂ rises ~6×
- Remove **Dsup** → DNA damage rises exactly 1/0.6 = 67 % (the 40 % cap)
- Remove **Nrf2** → glutathione depletes ~2× faster

## What is trustworthy vs. what is illustrative

- **Trustworthy:** the energy budget (first-principles), the dose→ROS scaling
  (standard radiolysis), and the *relative* roles of the defense systems.
- **Illustrative / order-of-magnitude:** absolute DNA-damage units, the melanin
  energy-transduction efficiency (assumed generous), and the FBA flux
  magnitudes. None of these change the headline, which holds across the full
  Monte Carlo and efficiency range.

## Limitations

- Lumped reactions (glycolysis, ETC, melanin synthesis) trade mechanistic
  detail for tractability.
- The metabolic model optimizes ATP maintenance, not growth/biomass.
- The Dsup structural assessment is a literature synthesis, not a folding or
  docking calculation.
- ODE parameters are single literature point estimates; only the FBA arm has
  formal uncertainty propagation.

## Falsifiable predictions (for wet-lab validation)

1. Melanized human cells under survivable chronic radiation will show **no
   measurable ATP or growth advantage** attributable to energy capture.
2. Any survival benefit from melanin + Dsup/Mn-AOX/Nrf2 will track **radiation
   tolerance markers** (γH2AX foci, ROS levels), not bioenergetic markers
   (ATP/ADP ratio, oxygen consumption).
3. Dsup will reduce DNA damage by ~40 %, independent of melanin.

## Key references

Buxton et al. 1988 (radiolysis G-values); Dadachova et al. 2007 (fungal
radiotrophy); Hashimoto et al. 2016 (Dsup); Daly et al. 2004 (Mn-antioxidant);
Lewis et al. 2015 (Nrf2); Favaudon et al. 2014 (FLASH regime); Kolesnikova et
al. 2023 (SOD2). Full parameter provenance is in the source docstrings.

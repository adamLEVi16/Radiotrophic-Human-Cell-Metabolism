# Revision notes: how the conclusion changed

`plan.pdf` is the original project proposal. Its framing — engineer melanin-based
radiotrophic metabolism into human cells and gain usable energy — is **not**
supported by the analysis in this repository. This note records what changed and
why, so the repo is not self-contradictory. The proposal is kept for provenance;
the README reflects the current, evidence-based conclusion.

## What the original proposal expected

- Melanin captures ionizing-radiation energy and drives NADH production, giving
  engineered human cells a metabolic **energy** advantage.
- Cross-species radioprotection genes (Dsup, Mn-antioxidant, Nrf2) manage the
  resulting ROS, making the system viable.
- An earlier model reported a **+18 % ATP boost** as supporting evidence.

## What the rigorous analysis found

1. **The energy is negligible (decisive).** By conservation of energy, absorbed
   radiation supplies ~10⁻⁶–10⁻¹⁰ % of a cell's ATP demand at every survivable
   dose rate (background, ISS, Chernobyl hotspots). Break-even requires
   ~0.3–28 Gy/s — a promptly lethal dose rate. See `energy_budget.py`,
   `fig1_energy_budget.png`. This is an upper bound: it holds even at 100 %
   energy-conversion efficiency and regardless of the biological mechanism.

2. **The "+18 % boost" was an artifact.** The radiotrophic reaction generated
   reducing equivalents with no link to the energy actually deposited by
   radiation; its flux was bounded only by an arbitrary cap. Once the flux is
   constrained by available radiation energy (`energy_constrained_dose`), the
   boost at any real dose rate is ~0 (`fig5_energy_constrained.png`).

3. **Two modeling defects were corrected.**
   - The kinetic ROS model's source term was ~6 orders of magnitude below
     scavenging capacity, so every dose and knockout produced identical output.
     It was rewritten to be dose-driven and now discriminates defenses.
   - Perturbing the radiotrophic reaction had broken ATP/ADP/Pi conservation
     (undetected because the metabolites carry no elemental formulas). This is
     now balanced and guarded by a test.

## What still stands from the original idea

- The **radioprotection genes remain independently useful** for radiation
  *tolerance* (shielding and repair) — just not for energy metabolism. Dsup
  reduces DNA damage ~40 % (now enforced as a cap); SOD is the ROS bottleneck.
- The distinction that matters: radiotrophy conflated **radiation tolerance**
  with **radiation-powered metabolism**. The first is real and useful; the
  second is energetically impossible for human cells at survivable doses.

## Bottom line

The project's value is now a **quantitative negative result**: a first-principles
demonstration of why radiotrophic metabolism cannot meaningfully power human
cells, with the radioprotection toolkit validated separately. That is a stronger
and more defensible contribution than the original optimistic framing.

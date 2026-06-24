# Paper Plan — Radiotrophic Human Cell Metabolism

Working planning document for turning this project into a written paper.

## 1. Honest readiness assessment

**Is it ready to submit as a peer-reviewed research paper today? No.**
**Is it a real, defensible foundation for one? Yes.**

What we have is a self-consistent *computational hypothesis study* with a clear,
novel framing (radiation as a survivable energy input for human cells, not just a
hazard) and a clean central result (the RESOURCE → STRAINED → LETHAL arc, with SOD
as the limiting defense). That is genuinely publishable **as a hypothesis /
methods / modelling paper in a preprint or specialist venue** — provided we are
upfront about the limitations rather than overselling.

What blocks it from a high-impact experimental journal:
1. **No physical dimensioning of the radiation axis** (flux ≠ Gy). This is the
   biggest single credibility gap and we should fix or fully own it.
2. **No wet-lab validation** — everything is in silico; comparisons are to
   third-party literature, not our own experiments.
3. **No external/independent model validation** — the network is bespoke and
   uses lumped reactions with several hand-set coefficients.

**Recommended target:** a methods/modelling preprint (bioRxiv) and/or submission
to a venue that accepts computational feasibility / hypothesis work, e.g.
*PLOS ONE*, *Frontiers in Microbiology / Bioengineering*, *Scientific Reports*,
or *Journal of Theoretical Biology*. Framed as **"a constraint-based feasibility
analysis,"** not "we showed human cells are radiotrophic."

## 2. The honest one-sentence claim

> Using a constraint-based metabolic model, we show that melanin-based
> radiotrophy *could* make a human cell treat ionizing radiation as a net energy
> source up to a defense-limited dose, beyond which engineered antioxidant
> capacity — chiefly SOD — sets a hard survival ceiling.

Note the modal verb. The paper is about *feasibility and limits*, not a
demonstrated phenotype.

## 3. Gap analysis (done → needed)

| Area | Have | Needed for paper |
|---|---|---|
| Core model | ✅ 56-rxn FBA, reproducible | Sensitivity/robustness sweep over uncertain params (partly done: ROS coeff) |
| Central result | ✅ 3-regime arc + figure | Same result under varied assumptions to show it's not an artifact |
| Dsup biology | ✅ 40% cap, feasibility score | Fine as supporting analysis |
| Dose units | ⚠️ illustrative anchor only | Either rigorous re-dimensioning OR explicit "relative units" framing |
| Validation | ⚠️ literature comparison table | Frame honestly; ideally add 1+ independent model cross-check |
| Wet-lab | ❌ none | Out of scope for v1; propose as future work |
| Kinetic model | ✅ exists | Integrate its dynamics into the narrative (transients, recovery) |
| Tests/repro | ✅ pinned solver, 8 invariant tests | Mention reproducibility in Methods |

## 4. Proposed structure

1. **Abstract** — feasibility framing, the RESOURCE/STRAINED/LETHAL result, SOD
   ceiling, explicit "in silico, relative units" caveat.
2. **Introduction** — radiotrophy in fungi; why human cells only lose to
   radiation; the engineering question; H₀ vs H₁ as stated in the README.
3. **Methods**
   - Constraint-based model construction (compartments, reactions, objective).
   - ROS production and the four defense systems; literature-grounded parameters.
   - Key constraints: Dsup 40% coupling cap; finite glutathione pool.
   - Forced-dose protocol; reproducibility (pinned GLPK, versions, tests).
   - Calibration: state plainly that fluxes are relative and the Gy axis is an
     anchored illustration (this *strengthens* credibility).
4. **Results**
   - R1: Dose–response and the three-regime arc (main figure).
   - R2: Defense ablation → SOD as single point of failure.
   - R3: Dsup capped vs uncapped → why the constraint matters.
   - R4: Bottleneck relief → the limit is OH-neutralisation (Dsup 40% + finite
     GSH), NOT SOD; MnSOD2 alone does nothing, relieving GSH lifts the ceiling.
   - R5: Combined stress (glucose/O₂ starvation) → radiotrophy as a survival
     supplement; the zero-glucose "lives on radiation alone" case.
   - R6: ROS-coefficient sensitivity → result is robust to the most uncertain param.
5. **Discussion** — what feasibility means here; comparison to fungal data; the
   engineering shopping list (genes, delivery); why SOD must be co-engineered.
6. **Limitations** — the four from the README, stated without flinching.
7. **Conclusion / future work** — re-dimensioning, kinetic survival modelling,
   and the wet-lab validation path.

## 5. Figures & tables (mostly already generated)

- **Fig 1** ✅ `radiation_atp_lethality.png` — ATP + lesions vs dose, three regimes.
- **Fig 2** ✅ `fig_ablation.png` — defense knockouts; SOD & OH-scavenging are
  single points of failure.
- **Fig 3** ✅ `fig_bottleneck_relief.png` — lifting the ceiling; MnSOD2 alone
  has no effect, GSH relief does.
- **Fig 4** ✅ `fig_ros_sensitivity.png` — net benefit robust across the ROS-cost
  range; lesions appear only at low ROS cost / high flux.
- **Fig 5** ⬜ Kinetic transients/recovery (from `kinetic_model.py` outputs).
- **Table 1** ✅ Gene-transfer feasibility (`gene_comparison.csv` / Dsup scores).
- **Table 2** ✅ Literature validation (`experimental_validation.csv`, reframed).

## 6. Roadmap to a submittable draft

- **Milestone A — Solidify the science (highest priority)**
  - ✅ Ablation / bottleneck / sensitivity plots (Figs 2–4) — `plot_experiments.py`.
  - ✅ ROS-coefficient sensitivity (Fig 4) shows the net benefit is robust across
    the most uncertain parameter.
  - ⬜ Decide the dose-units stance: commit to "relative units + anchored
    illustration," or attempt a defensible re-dimensioning. Recommend the former
    for v1.
  - ⬜ One more robustness check: confirm the 3-regime arc survives variation in
    GSH_SCAV_CAP and the Dsup fraction (not just the ROS coefficient).
- **Milestone B — Write Methods + Results** around the existing figures.
  ✅ Full draft in `MANUSCRIPT.md` (Abstract → Conclusion, all 5 figures + 2
  tables cited, honest Limitations section).
- **Milestone C — Write Intro/Discussion/Limitations**, get the framing right
  (feasibility, not phenotype). ✅ Included in `MANUSCRIPT.md`.
- **Milestone D — Internal review pass**, then bioRxiv preprint.
- **Milestone E (optional, major)** — design a wet-lab validation proposal
  (melanized HEK293 + Dsup under controlled dose) as a follow-up grant/paper.
  ✅ Drafted in `WETLAB_VALIDATION.md` — the decisive design isolates energy
  *capture* (the novel claim) from radio*protection* (already known) via a
  melanin × radiation × glucose factorial.

## 7. Top risks / objections to pre-empt

- *"The dose axis is meaningless."* → Own it; report relative units; the arc and
  the SOD ceiling are qualitative claims that don't depend on absolute Gy.
- *"FBA can't model damage/survival."* → That's why we added the finite-GSH cap
  and the forced-dose protocol; pair with the kinetic model for dynamics.
- *"Lumped reactions hide infeasibility."* → Report the sensitivity analysis and
  the invariant tests; release code for reproduction.
- *"This is speculative bioengineering."* → Correct, and we frame it as
  feasibility, with an explicit gene/delivery shopping list and known risks
  (e.g. Nrf2 oncogenicity, Dsup neurotoxicity).

# Engineering Radiotrophic Metabolism in Human Cells: A Constraint-Based Feasibility Analysis

**Adam Labban**

*Preprint draft — computational feasibility study. In silico results; no experimental validation. See Limitations.*

---

## Abstract

Melanized fungi colonizing high-radiation environments such as the Chernobyl
reactor appear to use ionizing radiation as a metabolic resource — a phenomenon
termed radiotrophy. Whether the same capacity could be engineered into human
cells is unknown. Here we ask a deliberately narrow question: not whether
radiotrophy could out-compete glucose metabolism, but whether a human cell could
be engineered to treat ionizing radiation as a **usable, survivable energy input
rather than purely as damage.** We built a constraint-based (flux-balance)
metabolic model of a human cell (56 reactions, 52 metabolites) augmented with a
melanin-driven radiotrophic NADH-generating pathway and a layered set of native
and cross-species engineered reactive-oxygen-species (ROS) defenses (superoxide
dismutase, catalase, glutathione peroxidase, tardigrade Dsup, *Deinococcus*
Mn-antioxidant, Nrf2, and overexpressed MnSOD2). Under a forced-dose protocol,
in which the cell cannot decline the radiation it absorbs, the model produces a
clear three-regime response: a **RESOURCE** regime in which radiation yields net
ATP with the radical load fully neutralized and zero modeled DNA lesions; a
**STRAINED** regime in which defenses saturate, DNA lesions accrue, and net ATP
declines while remaining positive; and a **LETHAL** regime in which no feasible
metabolic steady state exists (interpreted cautiously as loss of viability).
Defense-ablation analysis identifies
superoxide dismutase and glutathione-mediated radical scavenging as the two
single points of failure, and a targeted analysis shows that **superoxide, not
hydroxyl radical, is the species that overwhelms the cell first**, because
hydroxyl damage can be paid down through ATP-dependent repair whereas superoxide
has only the capacity-limited SOD exit. The qualitative result is robust to the
most uncertain model parameter. These findings frame a concrete, testable
hypothesis and a wet-lab experiment designed to isolate energy capture from
radioprotection. We emphasize the boundaries of the claim: the central
radiotrophic reaction is assumed rather than demonstrated, the model's radiation
axis is not physically dimensioned, and no result here has been experimentally
validated. The model tests the consequences of the radiotrophy hypothesis; it
does not establish that biological radiotrophy is achievable in human cells.

---

## 1. Introduction

Ionizing radiation is, for essentially all human cells, a pure liability: it
generates reactive oxygen species and DNA lesions, drives mutagenesis, and at
sufficient dose causes death. Yet certain melanized fungi not only tolerate
intense radiation but appear to grow *better* under it. Isolates recovered from
the walls of the damaged Chernobyl reactor, and later flown on the International
Space Station, show enhanced growth under elevated radiation, an effect linked to
their melanin content (Dadachova et al. 2007; Shunk et al. 2022). The proposed
mechanism — "radiotrophy" — is that melanin acts as a broadband absorber that
transduces radiation energy into a usable form, increasing the reducing capacity
(NADH) available to metabolism (Dadachova et al. 2007; Turick et al. 2011).

The remarkable claim of radiotrophy is not that it is an *efficient* energy
source — it is not, and melanized cells do not outgrow well-fed controls on rich
media. The remarkable claim is qualitative: that radiation, normally a hazard,
can serve as a metabolic *input at all*. For human cells the corresponding
question is sharper still, because a baseline human cell has no positive use for
radiation whatsoever. Could one engineer a human cell in which radiation becomes
net-tolerable — even net-beneficial — instead of strictly destructive?

This question sits at the intersection of two literatures that have so far
remained largely separate. The first is **radioprotection**: melanin shielding,
the tardigrade damage-suppressor protein Dsup (Hashimoto et al. 2016; Chavez et
al. 2019), the *Deinococcus radiodurans* manganese-antioxidant system (Daly et
al. 2004), and antioxidant-pathway activation (Lewis et al. 2015) all reduce
radiation harm and are, individually, well established. The second is
**radiotrophy**: energy *capture* from radiation, demonstrated phenomenologically
in fungi but mechanistically contested. Combining them — asking whether a human
cell could both survive radiation *and* harvest energy from it — has not, to our
knowledge, been modeled.

We formalize the question with two hypotheses:

- **Null (H₀):** Radiation is a strict liability for an engineered human cell;
  any energy captured by melanin is outweighed by the ROS and DNA damage it
  generates, so the cell is no better off, or dies.
- **Alternative (H₁):** Melanin plus engineered ROS defenses allow a human cell
  to route radiation energy into ATP while keeping the resulting damage
  survivable, up to some dose ceiling.

We test these with a constraint-based metabolic model. Constraint-based modeling
is well suited to a feasibility question: it asks what a metabolic network *can*
do at steady state under stoichiometric and capacity constraints, without
requiring kinetic parameters that are unavailable for an organism that does not
yet exist. We stress at the outset that this is an **in silico feasibility study**:
it can show that a hypothesis is internally consistent and worth testing, not
that it is true of real cells.

---

## 2. Methods

### 2.1 Model construction

The model is a constraint-based reconstruction built with COBRApy, comprising 56
reactions and 52 metabolites across three compartments (cytosol, mitochondrion,
extracellular). Core energy metabolism is represented by lumped glycolysis, a
full tricarboxylic-acid cycle, and a lumped electron transport chain with a
physiological P/O ratio and a small (~1%) electron leak to superoxide. Melanin
biosynthesis is modeled from tyrosine via tyrosinase and a lumped melanin-
synthesis reaction. The objective is ATP maintenance (a non-growth ATP demand),
maximized by flux-balance analysis.

This is a small, purpose-built reconstruction rather than a genome-scale human
model such as Recon3D. Reactions were selected to cover only the pathways relevant
to the hypothesis (central energy metabolism, ROS production and defense, and
melanin synthesis); the large majority of human metabolism (biosynthesis, nutrient
signaling, alternative substrates, compartment detail) is deliberately excluded.
The motivation is transparency: every reaction and capacity constraint can be
inspected and justified individually, which a genome-scale model does not permit.
The cost is metabolic coverage, and this restricted scope is a central limitation
(Section 5); a genome-scale re-grounding is proposed as future work.

### 2.2 The radiotrophic pathway and ROS production

The novel reaction, RADIO, represents melanin-mediated transduction of radiation
into reducing power: it consumes melanin and NAD⁺ and produces NADH, with reactive
oxygen species as obligate byproducts. Per unit RADIO flux it generates 0.28
superoxide and 0.24 hydroxyl radical; these coefficients derive from water-
radiolysis G-values (Buxton et al. 1988) adjusted for partial melanin quenching
(Schweitzer et al. 2009), and a small ATP overhead per NADH represents the cost
of energy transduction (consistent with the ATP decline reported by Bryan et al.
2011). The ROS-per-flux constants are enforced to match the reaction
stoichiometry by a build-time assertion, so they cannot silently drift.

RADIO represents a *hypothesized* transduction mechanism, not an experimentally
confirmed biochemical reaction. The change in melanin's electronic properties
under irradiation is supported (Dadachova et al. 2007; Turick et al. 2011), but
the direct conversion of that energy into cellular reducing equivalents (NADH) is
far less established. This reaction is therefore the central assumption of the
study rather than one of its conclusions: the model tests what follows *if* such a
transduction step exists, and to that extent assumes the phenomenon it
investigates. All downstream results are conditional on this assumption.

### 2.3 ROS defenses

Native defenses (superoxide dismutase, catalase, glutathione peroxidase and
reductase, glutathione-mediated hydroxyl scavenging) are capacity-capped to
finite physiological pools. Cross-species engineered defenses are added as
modular reactions: the tardigrade **Dsup** protein (hydroxyl interception at
chromatin), the *Deinococcus* **Mn-antioxidant** complex, an enhanced **Nrf2**
glutathione-recycling pathway, and overexpressed **MnSOD2**.

Two constraints encode the key biology. First, Dsup is capped at intercepting at
most **40% of generated hydroxyl radicals** — matching the ~40% DNA-damage
reduction measured by Hashimoto et al. (2016) — via a coupling constraint linking
Dsup flux to RADIO flux; without it, flux-balance optimization routes 100% of
radicals through the near-free Dsup reaction and pays no DNA-repair cost,
overstating its benefit. Second, glutathione scavenging is given a **finite
capacity**, so that beyond it hydroxyl radicals must either be intercepted by
Dsup or become DNA lesions repaired by base-excision repair at a cost of ~4 ATP
per lesion (Lindahl & Barnes 2000). This finite cap is what allows DNA damage to
occur at high dose, making the survival question testable rather than assumed
away. The overall mechanism is summarized in **Figure 1**.

### 2.4 Forced-dose protocol

A cell in a radiation field cannot decline the dose it absorbs. To represent this,
the central experiment fixes RADIO flux to a series of values (rather than letting
the optimizer choose it) and asks, at each: how much net ATP is produced relative
to zero dose; how the radical load partitions among interception, scavenging, and
DNA lesions; and whether any feasible steady state exists at all. Infeasibility is
interpreted as the cell being unable to balance the radical load, that is, loss of
viability (with the caveat in Section 5 that this is not equivalent to biological
cell death).

### 2.5 Reproducibility

The linear-program solver is pinned to GLPK (deterministic), all model parameters
are named constants, and dependency versions are pinned. Outputs are byte-
identical across runs, and a suite of eight invariant tests guards the central
findings (baseline ATP, the Dsup 40% cap, finite scavenging, the three-regime
arc, the SOD single-point-of-failure, and determinism). All code and data are
maintained in a version-controlled project repository (see Data and code
availability).

### 2.6 Dose calibration (illustrative only)

We attempted to convert the model's flux axis to absorbed dose (Gray) using
radiolysis G-values. This first-principles calibration overshoots realistic
lethal doses by ~10⁶-fold (about six orders of magnitude), confirming that the
model's fluxes are *relative*, not
physically dimensioned. We therefore report an explicitly **illustrative,
anchored** scale, matching the model's lethal threshold to a representative acute
mammalian lethal dose (~10 Gy); absolute Gy values in figures are approximate
analogies, not measurements.

---

## 3. Results

### 3.1 Radiation as a survivable energy source: the three-regime response

Under the forced-dose protocol, the model produces a characteristic three-regime
response (**Figure 2**):

- **RESOURCE** (low dose): ATP production rises monotonically with dose, the
  radical load is **100% neutralized**, and **no modeled DNA lesions** occur. Net ATP
  gain relative to the unirradiated baseline (111.4) reaches +4.3 at the ceiling
  (peak ATP ≈ 115.8). Here radiation is unambiguously a usable, harm-free energy
  input.
- **STRAINED** (intermediate dose): defenses saturate. The fraction of radicals
  neutralized falls (to ~85% at the upper end), DNA lesions appear and are
  repaired at ATP cost, and net ATP declines but remains positive. Radiation is
  still net-beneficial, but increasingly costly.
- **LETHAL** (high dose): no feasible metabolic steady state exists; the ROS load
  cannot be balanced. We interpret this cautiously as loss of viability, noting
  that flux-balance infeasibility is a steady-state statement and is not equivalent
  to biological cell death, which involves dynamics and adaptive responses the
  model does not capture (Section 5).

This arc is the central result: it shows that, within the model, radiation can be
a net-positive, survivable input **up to a defense-limited ceiling**, beyond which
it becomes first costly and then fatal. On the illustrative anchored scale, the
ceiling and lethal threshold correspond to roughly 7 and 10 Gy respectively (with
the strong caveat of Section 2.6).

### 3.2 Radiotrophy as a survival supplement under nutrient stress

The relative benefit of radiotrophy grows under energy limitation. At standard
glucose and oxygen, the radiotrophic cell gains ~4% ATP over a non-radiotrophic
control; under combined glucose/oxygen restriction this rises to ~12%. Most
strikingly, **under zero glucose the non-radiotrophic cell produces no ATP while
the radiotrophic cell continues to produce a substantial amount** (≈29 flux units
at moderate oxygen) — i.e., it subsists on radiation alone. This positions
radiotrophy as a survival supplement that matters most precisely when
conventional metabolism fails, consistent with the environments where natural
radiotrophs are found.

### 3.3 Superoxide dismutase and glutathione scavenging are single points of failure

Systematically knocking out each defense (**Figure 3**) shows that radiotrophy is
robust to loss of most individual defenses — removing the Mn-antioxidant complex,
Nrf2, or catalase has no effect on the achievable ceiling, and removing Dsup
merely reduces it (from ~21 to ~12.5 flux). Two knockouts, however, **collapse
radiotrophy entirely**: superoxide dismutase and glutathione-mediated hydroxyl
scavenging. Each is a single point of failure. This identifies the load-bearing
components of the engineered system and, by extension, the highest-priority
targets for any engineering effort.

### 3.4 Superoxide, not hydroxyl radical, overwhelms the cell first

Given that both superoxide and hydroxyl are produced, which one sets the lethal
threshold? A targeted analysis is decisive: augmenting superoxide handling
(adding MnSOD2) raises the lethal threshold from a flux of 29 to 59 — nearly
doubling the survivable dose — whereas a twentyfold increase in glutathione
scavenging capacity leaves it unchanged. **Superoxide is the species that
overwhelms the cell first.** The reason is structural rather than quantitative:
hydroxyl radicals possess an "escape valve" — unneutralized hydroxyl becomes DNA
damage that is repaired at ATP cost, a route that is always available — whereas
superoxide's only exit is the capacity-limited SOD reaction. Once superoxide
production exceeds total SOD capacity, no steady state exists. This makes SOD
augmentation the single most valuable intervention for raising the survivable
dose, and is annotated in the mechanism schematic (Figure 1).

### 3.5 Relieving the binding bottleneck lifts the ceiling

Because the radiotrophic ceiling is set by radical-neutralization capacity, the
correct engineering levers are those that relieve it (**Figure 4**). Adding
MnSOD2 *alone* produces no improvement — confirming that SOD is not the binding
constraint at the baseline ceiling — whereas relieving glutathione (OH-
neutralization) capacity lifts the ceiling from ~21 to ~29 flux. Relieving both
ROS bottlenecks together raises it to ~42, and a combined intervention (expanded
scavenging, MnSOD2, and synthetic melanin loading to bypass the biosynthetic
cost) reaches ~69 flux with ATP rising to ~151. The lesson is that interventions
must target the *binding* constraint: the same MnSOD2 that is useless at baseline
becomes essential once scavenging is relieved and superoxide becomes limiting.

### 3.6 The result is robust to the most uncertain parameter

The superoxide-per-NADH coefficient is the least-constrained parameter in the
model. Varying it across a 40-fold range (**Figure 5**) leaves the qualitative
conclusion intact: the net ATP gain from radiation remains **positive across the
entire range**, including the literature-derived band. DNA lesions appear only at
the low-cost end, where high radiotrophic flux drives damage. The central claim —
that radiation can be a net-positive, survivable input — does not depend on the
precise value of the most uncertain parameter.

---

## 4. Discussion

Taken together, the results support the alternative hypothesis **within the
model**: under its stated assumptions, the model predicts that an engineered human
cell could treat radiation as a usable and survivable energy input rather than as
pure damage, up to a defense-limited dose. It is worth stating plainly what this
does and does not mean. The model demonstrates internal, mathematical feasibility
given the radiotrophic transduction assumption; it does not demonstrate that
biological radiotrophy is achievable in human cells, and a skeptical reading is
that it shows the hypothesis is self-consistent rather than correct. With that
boundary fixed, the contribution is threefold. First, it reframes radiation for engineered human
cells from a strictly-harmful agent to a conditional resource, and makes that
reframing quantitative through the RESOURCE/STRAINED/LETHAL arc. Second, it
identifies the load-bearing biology: radical-neutralization capacity, with
superoxide dismutase as the component that sets the survival ceiling and
glutathione scavenging as a co-essential second bottleneck. Third, it yields a
concrete engineering priority list — SOD/MnSOD2 augmentation first, then
glutathione capacity — that differs from what an intuition built on
radioprotection alone (e.g., a focus on Dsup) would suggest.

The model's relationship to existing experimental data is consistent under the
correct framing (Table 2). The ISS *Cladosporium* growth advantage (~21%; Shunk
et al. 2022) is of the same order of magnitude as the model's stress-condition ATP
benefit; we stress that growth rate and steady-state ATP are distinct quantities
and this comparison is only qualitative, not a like-for-like validation. The
ATP decline reported in melanized cells at high dose (Bryan et al. 2011)
corresponds to the model's STRAINED-to-LETHAL transition rather than contradicting
it. The Dsup 40% protection (Hashimoto et al. 2016) is imposed directly. Crucially,
the model does *not* reproduce the large fungal CFU boosts (Dadachova et al.
2007), and should not — those reflect growth and repair processes outside a
steady-state ATP model, and the present thesis is about tolerance and capture, not
boost magnitude.

Among the engineered components, Dsup is the most feasible single gene transfer:
it has been expressed in human cells with confirmed nuclear localization and DNA
protection, targets a universally conserved nucleosome surface, and requires no
cofactors (Table 1). Notably, however, the model shows Dsup is *protective but not
energy-determining* — it raises DNA integrity but does not set the survival
ceiling, which superoxide does. The components that most extend survivable dose
(SOD systems) are different from the one with the cleanest prior validation
(Dsup), a tension worth keeping in view for experimental design.

---

## 5. Limitations

These results must be read with several hard limitations.

1. **The radiation axis is not physically dimensioned.** Model fluxes are
   relative; the first-principles calibration fails by ~10⁶-fold. Absolute Gy values
   are illustrative anchors, not measurements, and no quantitative dose claim
   should be drawn from them.
2. **No experimental validation.** Every result is in silico. Comparisons are to
   third-party literature, not to experiments performed here. The model is a
   hypothesis generator, not evidence about real cells.
3. **The model is bespoke and simplified.** Glycolysis, the ETC, and melanin
   synthesis are lumped; several ROS and melanin stoichiometries are literature-
   guided estimates rather than measured values. It is not a validated genome-
   scale reconstruction.
4. **"Lethal" means "no feasible steady state."** This is a reasonable proxy for
   death but not a calibrated survival curve, and flux-balance analysis omits
   dynamics, transient damage spikes, protein and membrane oxidation, and chronic
   effects.
5. **The underlying premise is itself debated.** Whether melanin radiotrophy
   yields *usable* metabolic energy (versus radioprotection alone) is not settled
   even in fungi; the engineering hypothesis inherits that uncertainty.

---

## 6. Conclusion and future work

Under the assumptions of a constraint-based model, and conditional on a melanin
radiotrophic transduction step that is itself hypothetical, melanin-based
radiotrophy combined with engineered antioxidant defenses *could* allow a human
cell to treat ionizing radiation as a usable, survivable energy input up to a
defense-limited dose, beyond which superoxide overwhelms its dismutase capacity
and sets a hard survival ceiling. The result is internally consistent, robust to
the most uncertain parameter, and yields a specific engineering priority
(SOD/MnSOD2 first). It does not establish that biological radiotrophy is achievable
in human cells; it establishes that the hypothesis is self-consistent and
identifies what would have to be true, and experimentally tested, for it to hold.

The decisive next step is experimental. We outline a proof-of-concept design
(Supplementary: Wet-Lab Validation) built around isolating energy *capture* — the
novel claim — from radio*protection*, which is already established. The core is a
melanin × radiation × glucose factorial in which only the radiotrophic hypothesis
predicts a radiation- and melanin-dependent ATP rescue under energy limitation; a
shielded-sham control and a synthetic-melanin arm isolate the effect, and a clean
negative result would itself meaningfully constrain the hypothesis. Beyond that,
re-grounding the network in a validated genome-scale human reconstruction and a
physical dose calibration would convert the present qualitative claims into
quantitative ones.

---

## Figures

- **Figure 1.** Mechanism schematic: radiation → melanin → radiotrophic NADH → ATP,
  with ROS byproducts routed through native and engineered defenses to protect
  DNA; superoxide/SOD marked as the failure-first branch. (`fig_mechanism.png`)
- **Figure 2.** ATP production and DNA lesions vs. forced radiation dose, showing
  the RESOURCE / STRAINED / LETHAL regimes. (`radiation_atp_lethality.png`)
- **Figure 3.** Defense ablation: superoxide dismutase and glutathione scavenging
  are single points of failure. (`fig_ablation.png`)
- **Figure 4.** Bottleneck relief: relieving radical-neutralization capacity
  (not SOD alone) lifts the radiotrophic ceiling. (`fig_bottleneck_relief.png`)
- **Figure 5.** Robustness of the net benefit to the superoxide-per-NADH
  coefficient. (`fig_ros_sensitivity.png`)

## Tables

- **Table 1.** Cross-species gene-transfer feasibility (`gene_comparison.csv`,
  Dsup assessment).
- **Table 2.** Model predictions vs. published experimental observations
  (`experimental_validation.csv`).

---

## Data and code availability

All model code, experiment scripts, generated data, figures, the dose-calibration
module, the test suite, and the wet-lab validation proposal are maintained in a
version-controlled project repository and will be made publicly available upon
publication.

## References

1. Bryan R. et al. (2011) *Fungal Biology* 115:945.
2. Buxton G.V. et al. (1988) *J. Phys. Chem. Ref. Data* 17:513.
3. Chavez C. et al. (2019) *eLife* 8:e47682.
4. Dadachova E. et al. (2007) *PLoS ONE* 2:e457.
5. Daly M.J. et al. (2004) *Science* 306:1025.
6. Hashimoto T. et al. (2016) *Nature Communications* 7:12808.
7. Lewis K.N. et al. (2015) *PNAS* 112:3722.
8. Lindahl T. & Barnes D.E. (2000) *Cold Spring Harb. Symp. Quant. Biol.* 65:127.
9. Schweitzer A.D. et al. (2009) *Int. J. Radiat. Oncol. Biol. Phys.* 73:1494.
10. Shunk G.K. et al. (2022) *Frontiers in Microbiology* (ISS *Cladosporium*).
11. Turick C.E. et al. (2011) *Bioelectrochemistry*.
12. Kolesnikova et al. (2023) PMC10744337.
13. Engineered melanin nanoparticle radioprotection (2025) *Nature Communications*.

*(Reference list to be completed with full citation details at submission.)*

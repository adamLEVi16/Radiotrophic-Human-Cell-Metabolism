# Wet-Lab Validation Proposal

A proof-of-concept experimental plan to test the central prediction of the
computational model: that melanin + engineered antioxidant defenses can let a
human cell **capture usable energy from radiation while tolerating it**, not
merely survive it.

This is the experiment that would move the project from "internally consistent
hypothesis" to "experimentally supported finding." It is written to be dropped
into the paper's Discussion / Future Work and to seed a grant or a collaboration
with a lab that has irradiation facilities.

---

## 0. The one design principle that matters

Radio**protection** (melanin shielding, Dsup DNA protection, antioxidant
overexpression) is **already established** in the literature. If we only show
that melanized/Dsup cells survive radiation better, we have shown nothing new.

The novel, model-specific claim is **radiotrophy = net energy *capture* from
radiation.** Every experiment below is designed to *isolate energy capture from
protection*. The decisive signature is:

> Under energy limitation, ATP / proliferation **increases with radiation dose**
> (up to a ceiling) in melanized cells — a *positive* dose-dependent benefit that
> pure protection can never produce (protection can only reduce harm toward zero,
> never above the unirradiated baseline).

If we see that positive, radiation-dependent, melanin-dependent energy signal,
the radiotrophic hypothesis is supported. If melanized cells only ever do "less
badly" than controls, the result collapses to ordinary radioprotection.

---

## 1. Hypotheses (mapped to the model)

| # | Model prediction | Testable hypothesis |
|---|---|---|
| H1 | RESOURCE regime: net energy + no net DNA damage at low/moderate dose | Melanized cells under energy stress gain ATP/growth *with* dose; γH2AX not elevated vs unirradiated |
| H2 | Energy benefit requires radiation AND melanin | Benefit disappears if radiation is shielded OR melanin is absent |
| H3 | STRAINED→LETHAL arc at high dose | Above a threshold, DNA damage rises and viability falls steeply |
| H4 | SOD / OH-scavenging are single points of failure | Inhibiting SOD (or GSH synthesis) abolishes the benefit and sharply lowers the lethal threshold |
| H5 | Dsup is protective but marginal to *energy* capture | Dsup improves survival/DNA integrity but does **not** create the energy signal by itself |

---

## 2. Cell system & constructs

**Base line:** HEK293T — the line in which Dsup nuclear localization and DNA
protection are already validated (Hashimoto 2016), so it removes one variable.
Secondary line: an immortalized fibroblast or RPE-1 for a non-transformed check.

**Melanin — two independent routes (to avoid a construct-specific artifact):**
- *Biosynthetic:* lentiviral **tyrosinase (TYR)** ± melanosomal genes →
  endogenous melanin; tune with tyrosine in media.
- *Synthetic:* load **melanin nanoparticles** (Nat. Commun. 2025 approach) into
  wild-type cells. This bypasses biosynthesis and decouples "melanin present"
  from "cell engineered" — a clean orthogonal control.

**Defense constructs (modular, so each can be ablated):**
- Tardigrade **Dsup** (DNA protection).
- **SOD2** overexpression (the predicted bottleneck).
- Optional: Nrf2 pathway activator.

**Arms (minimum viable factorial):**
1. WT (no melanin, no defenses) — baseline
2. + melanin only
3. + melanin + Dsup
4. + melanin + Dsup + SOD2 (full engineered model)
5. + Dsup only (no melanin) — protection-without-capture control
6. Arm 4 + SOD inhibitor / GSH-synthesis inhibitor (BSO) — bottleneck test

---

## 3. Radiation

- **Source:** Cs-137 / X-ray irradiator, or a low-dose-rate chronic source.
- **Two modes** (the model distinguishes them):
  - **Chronic low-dose-rate** (the radiotrophy-relevant regime; mirrors the
    Chernobyl/ISS context) — primary.
  - **Acute dose-response** (0–10+ Gy) to map the STRAINED→LETHAL arc — secondary.
- **Critical control: shielded sham.** Identical handling, radiation blocked.
  This is what isolates a radiation-dependent effect from construct effects.

---

## 4. Readouts

| Quantity | Assay | Tests |
|---|---|---|
| ATP / energy charge | CellTiter-Glo; LC-MS for ATP/ADP/AMP | H1, H2 (the core energy signal) |
| Proliferation / survival | clonogenic assay, growth curves | H1, H3 |
| DNA damage | γH2AX foci, alkaline comet assay | H1, H3, H5 |
| ROS | DCFDA (general), MitoSOX (superoxide) | H4 |
| Metabolic flux | Seahorse (OCR/ECAR); ¹³C-glucose tracing | H1, H2 — does radiation offset glucose demand? |
| Melanin content | spectrophotometry (405 nm), Fontana-Masson | dose normalization |

**The killer plot** (direct analog of `radiation_atp_lethality.png`):
ATP (and proliferation) vs. radiation dose, under **glucose-limited** medium,
for melanized vs non-melanized cells. A radiotrophic signal = the melanized
curve rises above its own unirradiated baseline and above the non-melanized
curve, in a melanin- and radiation-dependent way.

---

## 5. The decisive 2×2×2 experiment

Factors: **melanin** (±) × **radiation** (± shielded) × **glucose** (normal /
limited). Readout: ATP and proliferation.

Predicted radiotrophic signature (and *only* radiotrophy predicts the starred cell):

| melanin | radiation | glucose | predicted ATP/growth |
|---|---|---|---|
| – | – | limited | low (starving) |
| – | + | limited | low (radiation = damage only) |
| + | – | limited | low (melanin alone ≠ energy) |
| **+** | **+** | **limited** | **rescued / elevated ★** |
| + | + | normal | small or no benefit (glucose already sufficient) |

The star cell being elevated — and *nothing else* — is the experiment. If the
benefit shows up without radiation, or without melanin, the hypothesis fails.

---

## 6. Controls & rigor

- Shielded-sham radiation control (isolates radiation effect).
- Synthetic-melanin arm (isolates melanin from the engineering).
- Melanin dose-normalization (benefit should scale with melanin content).
- Blinded image analysis for γH2AX/comet.
- Power analysis up front; biological triplicate minimum; pre-registered
  primary endpoint (ATP under limited glucose, star cell vs. others).
- Metabolic-inhibitor controls (rotenone/oligomycin) to confirm ATP source.

---

## 7. Phasing

- **Phase 1 (proof of signal, ~3–6 mo):** Arms 1–4, chronic low-dose, the 2×2×2.
  Go/no-go: does the star cell show a radiation- and melanin-dependent ATP rescue?
- **Phase 2 (mechanism, ~6–9 mo):** dose-response arc (H3), bottleneck ablation
  (H4, SOD/BSO), Dsup decoupling (H5), Seahorse/¹³C flux.
- **Phase 3 (robustness):** second cell line; chronic long-term viability;
  karyotype/transformation safety (Nrf2 oncogenicity, Dsup chromatin effects).

---

## 8. Risks & honest caveats

- **The effect may be too small to detect.** The model's net energy gain is
  modest (single-digit %); chronic low-dose energy capture may be below assay
  noise. Mitigation: energy-limited conditions amplify the relative signal
  (per the model's combined-stress result), and ¹³C tracing is more sensitive
  than bulk ATP.
- **Radioprotection confound.** The single biggest threat to interpretation;
  the 2×2×2 design and the "Dsup-only" arm exist specifically to control it.
- **Melanin phototoxicity / ROS burden** may dominate at high melanin loads.
- **It may simply not work** — the underlying biological premise (that melanin
  radiotrophy yields usable energy) is still debated even in fungi. A clean
  negative result is itself publishable and would constrain the hypothesis.

---

## 9. What a positive result would mean

A radiation- and melanin-dependent energy rescue under limitation would be the
first evidence that a human cell can be engineered to treat ionizing radiation
as a usable energy input rather than purely as damage — turning the
computational hypothesis into an experimental finding and justifying the
genome-scale modeling and translational follow-up.

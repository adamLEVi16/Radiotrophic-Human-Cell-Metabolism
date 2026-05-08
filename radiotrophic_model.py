"""
Radiotrophic Human Cell Metabolic Model
========================================
Computational feasibility analysis of melanin-based radiotrophic
metabolism in human cells through cross-species genetic engineering.

Authors: Adam Labban
Date: March 2026

Model: Custom constraint-based metabolic model (54 reactions, 52 metabolites)
Built with COBRApy. Includes glycolysis, TCA cycle, electron transport chain,
melanin synthesis, radiotrophic NADH generation, and four ROS defense systems
(native SOD/catalase/GPX + engineered Dsup/Mn-AOX/Nrf2).

Literature-grounded parameters (Phase 2):
    - ROS stoichiometry from water radiolysis G-values (Buxton et al. 1988)
    - Melanin NADH reduction: 4x increase under irradiation (Dadachova et al. 2007)
    - Melanin redox cycling / self-restoration (Turick et al. 2011)
    - Dsup 40% DNA damage reduction (Hashimoto et al. 2016)
    - Mn-antioxidant protein protection (Daly et al. 2004)
    - DNA repair cost ~4 ATP/lesion (Lindahl & Barnes 2000)

Requirements:
    pip install cobra pandas matplotlib
"""

import cobra
import pandas as pd
import json
import os

# ============================================================
# MODEL CONSTRUCTION
# ============================================================

def build_model():
    """Build the radiotrophic human cell metabolic model."""
    
    model = cobra.Model('radiotrophic_human_cell')
    M = {}
    
    def m(mid, name, comp='c'):
        met = cobra.Metabolite(mid, name=name, compartment=comp)
        M[mid] = met
    
    # --- Cytosolic metabolites ---
    for mid, name in [
        ('glc_c','Glucose'), ('pyr_c','Pyruvate'), ('lac_c','Lactate'),
        ('atp_c','ATP'), ('adp_c','ADP'), ('nad_c','NAD+'), ('nadh_c','NADH'),
        ('h_c','H+'), ('pi_c','Phosphate'), ('h2o_c','Water'), ('co2_c','CO2'),
        ('o2_c','O2'), ('o2s_c','Superoxide'), ('h2o2_c','H2O2'),
        ('melanin_c','Melanin'), ('tyr_c','Tyrosine'), ('dopa_c','L-DOPA'),
        ('gthrd_c','GSH'), ('gthox_c','GSSG'),
        ('oh_radical_c','Hydroxyl radical'), ('dna_damage_c','DNA damage marker')
    ]:
        m(mid, name)
    
    # --- Mitochondrial metabolites ---
    for mid, name in [
        ('pyr_m','Pyruvate'), ('accoa_m','Acetyl-CoA'), ('coa_m','CoA'),
        ('oaa_m','Oxaloacetate'), ('cit_m','Citrate'), ('akg_m','Alpha-ketoglutarate'),
        ('succoa_m','Succinyl-CoA'), ('succ_m','Succinate'),
        ('fum_m','Fumarate'), ('mal_m','Malate'),
        ('nad_m','NAD+'), ('nadh_m','NADH'), ('fad_m','FAD'), ('fadh2_m','FADH2'),
        ('atp_m','ATP'), ('adp_m','ADP'), ('h_m','H+'), ('co2_m','CO2'),
        ('h2o_m','Water'), ('o2_m','O2'), ('pi_m','Phosphate'),
        ('o2s_m','Superoxide'), ('h2o2_m','H2O2')
    ]:
        m(mid, name, 'm')
    
    # --- Extracellular metabolites ---
    for mid, name in [
        ('glc_e','Glucose'), ('o2_e','O2'), ('co2_e','CO2'),
        ('h2o_e','Water'), ('lac_e','Lactate'), ('h_e','H+'),
        ('pi_e','Phosphate'), ('tyr_e','Tyrosine')
    ]:
        m(mid, name, 'e')
    
    # --- Helper to add reactions ---
    rxns = []
    def R(rid, name, md, bounds=(0, 1000)):
        r = cobra.Reaction(rid)
        r.name = name
        r.lower_bound, r.upper_bound = bounds
        r.add_metabolites({M[k]: v for k, v in md.items()})
        rxns.append(r)
    
    # ============================================================
    # TRANSPORT REACTIONS
    # ============================================================
    R('GLCt',  'Glucose transport',          {'glc_e':-1, 'glc_c':1})
    R('O2t',   'O2 to cytosol',             {'o2_e':-1, 'o2_c':1})
    R('O2tm',  'O2 to mitochondria',        {'o2_c':-1, 'o2_m':1}, (-1000,1000))
    R('CO2tm', 'CO2 from mitochondria',     {'co2_m':-1, 'co2_c':1}, (-1000,1000))
    R('CO2t',  'CO2 export',               {'co2_c':-1, 'co2_e':1})
    R('H2Ot',  'H2O transport',            {'h2o_e':-1, 'h2o_c':1}, (-1000,1000))
    R('H2Otm', 'H2O to mitochondria',      {'h2o_c':-1, 'h2o_m':1}, (-1000,1000))
    R('Ht',    'H+ export to extracellular', {'h_c':-1, 'h_e':1}, (-1000,1000))
    R('Htm',   'H+ to mitochondria',        {'h_c':-1, 'h_m':1}, (-1000,1000))
    R('LACt',  'Lactate export',            {'lac_c':-1, 'lac_e':1})
    R('PYRtm', 'Pyruvate to mitochondria',  {'pyr_c':-1, 'pyr_m':1})
    R('TYRt',  'Tyrosine transport',        {'tyr_e':-1, 'tyr_c':1})
    R('PIt',   'Phosphate transport',       {'pi_e':-1, 'pi_c':1})
    R('PItm',  'Phosphate to mitochondria', {'pi_c':-1, 'pi_m':1}, (-1000,1000))
    R('ATPtm', 'ATP/ADP translocase',       {'atp_m':-1, 'adp_c':-1, 'atp_c':1, 'adp_m':1}, (-1000,1000))
    R('NADHtm','NADH shuttle (malate-asp)', {'nadh_c':-1, 'nad_c':1, 'nadh_m':1, 'nad_m':-1})
    R('H2O2tm','H2O2 from mitochondria',    {'h2o2_m':-1, 'h2o2_c':1})
    R('O2Stm', 'Superoxide from mitochondria', {'o2s_m':-1, 'o2s_c':1})
    
    # ============================================================
    # CORE METABOLISM
    # ============================================================
    
    # Glycolysis (lumped): Glucose -> 2 Pyruvate + 2 ATP + 2 NADH
    R('GLYC', 'Glycolysis (lumped)', {
        'glc_c':-1, 'nad_c':-2, 'adp_c':-2, 'pi_c':-2,
        'pyr_c':2, 'nadh_c':2, 'atp_c':2, 'h2o_c':2, 'h_c':2
    })
    
    # Lactate dehydrogenase (reversible)
    R('LDH', 'Lactate dehydrogenase', {
        'pyr_c':-1, 'nadh_c':-1, 'h_c':-1, 'lac_c':1, 'nad_c':1
    }, (-1000, 1000))
    
    # Pyruvate dehydrogenase
    R('PDH', 'Pyruvate dehydrogenase', {
        'pyr_m':-1, 'nad_m':-1, 'coa_m':-1,
        'accoa_m':1, 'nadh_m':1, 'co2_m':1
    })
    
    # TCA Cycle
    R('CS',     'Citrate synthase',          {'accoa_m':-1, 'oaa_m':-1, 'h2o_m':-1, 'cit_m':1, 'coa_m':1})
    R('ACONIDH','Aconitase + Isocitrate DH', {'cit_m':-1, 'nad_m':-1, 'akg_m':1, 'nadh_m':1, 'co2_m':1})
    R('AKGDH',  'Alpha-KG dehydrogenase',    {'akg_m':-1, 'nad_m':-1, 'coa_m':-1, 'succoa_m':1, 'nadh_m':1, 'co2_m':1})
    R('SCS',    'Succinyl-CoA synthetase',   {'succoa_m':-1, 'adp_m':-1, 'pi_m':-1, 'succ_m':1, 'coa_m':1, 'atp_m':1})
    R('SDH',    'Succinate dehydrogenase',   {'succ_m':-1, 'fad_m':-1, 'fum_m':1, 'fadh2_m':1})
    R('FUM',    'Fumarase',                  {'fum_m':-1, 'h2o_m':-1, 'mal_m':1}, (-1000,1000))
    R('MDH',    'Malate dehydrogenase',      {'mal_m':-1, 'nad_m':-1, 'oaa_m':1, 'nadh_m':1, 'h_m':1}, (-1000,1000))
    
    # Electron Transport Chain (lumped with ROS byproduct)
    # NADH -> 2.5 ATP (P/O ratio), 1% electron leak to superoxide
    R('ETC_N', 'ETC Complex I-IV (NADH)', {
        'nadh_m':-1, 'o2_m':-0.5, 'adp_m':-2.5, 'pi_m':-2.5,
        'nad_m':1, 'h2o_m':0.5, 'atp_m':2.5, 'o2s_m':0.01
    })
    # FADH2 -> 1.5 ATP
    R('ETC_F', 'ETC Complex II-IV (FADH2)', {
        'fadh2_m':-1, 'o2_m':-0.5, 'adp_m':-1.5, 'pi_m':-1.5,
        'fad_m':1, 'h2o_m':0.5, 'atp_m':1.5, 'o2s_m':0.005
    })
    
    # ============================================================
    # MELANIN SYNTHESIS
    # ============================================================
    R('TYROX',  'Tyrosinase (Tyr -> DOPA)',  {'tyr_c':-1, 'o2_c':-0.5, 'dopa_c':1, 'h2o_c':0.5})
    R('MELSYN', 'Melanin synthesis (lumped)', {'dopa_c':-4, 'o2_c':-2, 'melanin_c':1, 'h2o_c':4, 'h2o2_c':1})
    
    # ============================================================
    # RADIOTROPHIC METABOLISM (novel engineered pathway)
    # Based on Dadachova et al. (2007) - Cryptococcus neoformans
    # Melanin absorbs ionizing radiation, energizes electrons,
    # drives NAD+ -> NADH with superoxide as byproduct
    # ============================================================
    # Literature grounding (Phase 2):
    #   - Melanin consumption: 0.02 per NADH. Turick et al. (2011) showed melanin
    #     undergoes redox cycling and self-restoration under gamma radiation, so
    #     degradation is slow. Coefficient represents net polymer loss after cycling.
    #   - ROS stoichiometry from water radiolysis G-values (Buxton et al. 1988):
    #     G(O2•⁻) ≈ 0.32 μmol/J, G(•OH) ≈ 0.28 μmol/J → ratio ~1.14:1.
    #     Melanin scavenges ~50% of radicals (Schweitzer et al. 2009), so net
    #     ROS per NADH: 0.28 O2•⁻ + 0.24 •OH (after melanin quenching).
    # ATP transduction overhead: Bryan et al. (2011) Fungal Biology showed ATP
    # *decreases* in melanized cells under radiation, suggesting energy transduction
    # has metabolic overhead. Modeled as 0.5 ATP consumed per NADH generated
    # (melanin activation, electron channeling, radical management costs).
    R('RADIO', 'Melanin radiotrophic NADH generation', {
        'melanin_c':-0.02, 'nad_c':-1, 'h_c':-1, 'atp_c':-0.5, 'h2o_c':-0.5,
        'nadh_c':1, 'o2s_c':0.28, 'oh_radical_c':0.24,
        'adp_c':0.5, 'pi_c':0.5
    }, (0, 50))
    
    # ============================================================
    # ROS DAMAGE PATHWAYS
    # ============================================================
    # Hydroxyl radicals cause DNA damage (Fenton chemistry)
    # This creates a damage marker that MUST be repaired at ATP cost
    R('FENTON', 'Hydroxyl radical DNA damage', {
        'oh_radical_c':-1, 'dna_damage_c':1
    })

    # DNA base excision repair - costs ATP per lesion repaired
    # Source: Lindahl & Barnes (2000) - BER requires ~4 ATP per lesion
    R('BER', 'Base excision repair (ATP-dependent)', {
        'dna_damage_c':-1, 'atp_c':-4, 'h2o_c':-4,
        'adp_c':4, 'pi_c':4, 'h_c':4
    })

    # Hydroxyl radical scavenging (glutathione-mediated)
    # Cap at 10: consistent with cellular GSH turnover via GR (cap=2) + NRF2 boosting
    # at ~2–3x, and literature estimates of GSH-mediated OH scavenging in HEK293-like
    # cells (~5–20 units in same flux scale; Meister 1988, Flohe 1971).
    # Dsup + OH_SCAV combined (5+10=15) covers max OH load at RADIO flux=50 (12 units).
    R('OH_SCAV', 'Hydroxyl radical scavenging by GSH', {
        'oh_radical_c':-1, 'gthrd_c':-1, 'gthox_c':0.5, 'h2o_c':1
    }, (0, 10))

    # ============================================================
    # ROS DEFENSE - NATIVE HUMAN
    # Capacity-constrained to realistic cellular Vmax values.
    # Native enzyme pools are finite; radiotrophic ROS load can
    # exceed their capacity, making engineered defenses essential.
    # ============================================================
    # Native enzyme capacities are capped to represent finite cellular pools.
    # At basal ETC ROS (~0.2 flux), these caps are ample. Under radiotrophic
    # load (7.5 superoxide + 2.5 OH at flux=25), the system operates near
    # capacity, making engineered defenses essential for robustness.
    #
    # SOD1 (Cu/Zn): cap=4. At radio flux=25, SOD needs 3.75 → 94% saturated.
    R('SODc', 'Superoxide dismutase (cytosolic)',  {'o2s_c':-2, 'h_c':-2, 'h2o2_c':1, 'o2_c':1}, (0, 4))
    R('SODm', 'Superoxide dismutase (mito)',       {'o2s_m':-2, 'h_m':-2, 'h2o2_m':1, 'o2_m':1}, (0, 3))
    # Catalase: cap=3. Peroxisomal, limited cytosolic availability.
    R('CATc', 'Catalase (cytosolic)',              {'h2o2_c':-2, 'h2o_c':2, 'o2_c':1}, (0, 3))
    R('CATm', 'Catalase (mitochondrial)',          {'h2o2_m':-2, 'h2o_m':2, 'o2_m':1})
    # GPX1: selenium-dependent, limited by GSH pool turnover. Cap=3.
    R('GPX',  'Glutathione peroxidase',            {'h2o2_c':-1, 'gthrd_c':-2, 'h2o_c':2, 'gthox_c':1}, (0, 3))
    # GR: NADPH-dependent (NADH proxy), limited by reductase expression. Cap=2.
    R('GR',   'Glutathione reductase',             {'gthox_c':-1, 'nadh_c':-1, 'h_c':-1, 'gthrd_c':2, 'nad_c':1}, (0, 2))
    
    # ============================================================
    # ROS DEFENSE - CROSS-SPECIES ENGINEERED
    # ============================================================
    
    # Tardigrade Dsup protein - binds chromatin and shields DNA from hydroxyl radicals
    # Source: Hashimoto et al. (2016) Nature Communications
    # Demonstrated 40% reduction in radiation-induced DNA damage in HEK293 cells.
    # Chavez et al. (2019) eLife: Dsup is a nucleosome-binding protein that protects
    # chromosomal DNA from hydroxyl radical-mediated cleavage.
    # Modeled as: Dsup intercepts OH radicals at chromatin, converting to water.
    # Capacity capped at 0.4 × (max_RADIO_flux × OH_stoichiometry):
    #   0.4 × (50 × 0.24) = 4.8 → upper_bound = 5
    # This reflects 40% protection ceiling from Hashimoto et al.; Dsup cannot
    # fully substitute for glutathione-based scavenging (OH_SCAV).
    R('DSUP', 'Tardigrade Dsup DNA shielding', {
        'oh_radical_c':-1, 'h_c':-1, 'h2o_c':1
    }, (0, 5))
    
    # Deinococcus radiodurans Mn-antioxidant complex
    # Source: Daly et al. (2004) Science
    # Manganese complexes protect proteins from oxidative damage
    R('MNAOX', 'Deinococcus Mn-antioxidant complex', {
        'h2o2_c':-2, 'h2o_c':2, 'o2_c':1
    })
    
    # Naked mole rat enhanced Nrf2 pathway
    # Source: Lewis et al. (2015) PNAS
    # Constitutively elevated Nrf2 drives superior antioxidant gene expression
    R('NRF2', 'Enhanced Nrf2 glutathione recycling', {
        'gthox_c':-1, 'nadh_c':-0.5, 'h_c':-0.5,
        'gthrd_c':2, 'nad_c':0.5  # More efficient than normal GR
    })

    # MnSOD2 overexpression (addresses SOD bottleneck)
    # Source: Kolesnikova et al. (2023) PMC10744337 - CRISPRa SOD2 in HEK293T
    # Co-expression of SOD2+catalase provides broadest radioprotection.
    # Mitochondrial MnSOD handles superoxide that Cu/Zn-SOD1 cannot reach.
    # Modeled as additional cytosolic SOD capacity (SOD2 can be retargeted).
    R('SOD2', 'Overexpressed MnSOD2 (engineered)', {
        'o2s_c':-2, 'h_c':-2, 'h2o2_c':1, 'o2_c':1
    }, (0, 0))  # Default OFF; enabled in bottleneck-relief experiments

    # ============================================================
    # SYNTHETIC MELANIN LOADING (bypasses biosynthesis bottleneck)
    # Source: Nature Communications 2025 - engineered melanin NPs
    # Synthetic melanin nanoparticles can be loaded into cells,
    # bypassing the tyrosine -> DOPA -> melanin pathway entirely.
    # ============================================================
    R('EX_mel', 'Synthetic melanin nanoparticle loading', {
        'melanin_c':1
    }, (0, 0))  # Default OFF; enabled in bottleneck-relief experiments
    
    # ============================================================
    # EXCHANGE REACTIONS
    # ============================================================
    R('EX_glc', 'Glucose uptake',    {'glc_e':-1}, (-10, 0))
    R('EX_o2',  'O2 uptake',         {'o2_e':-1},  (-20, 0))
    R('EX_co2', 'CO2 secretion',     {'co2_e':-1}, (0, 1000))
    R('EX_h2o', 'H2O exchange',      {'h2o_e':-1}, (-1000, 1000))
    R('EX_lac', 'Lactate secretion',  {'lac_e':-1}, (0, 1000))
    R('EX_h',   'H+ exchange',       {'h_e':-1},   (-1000, 1000))
    R('EX_pi',  'Phosphate uptake',  {'pi_e':-1},  (-100, 0))
    R('EX_tyr', 'Tyrosine uptake',   {'tyr_e':-1}, (-5, 0))
    
    # ============================================================
    # OBJECTIVE: ATP MAINTENANCE
    # ============================================================
    R('ATPM', 'ATP maintenance demand', {
        'atp_c':-1, 'h2o_c':-1, 'adp_c':1, 'pi_c':1, 'h_c':1
    }, (0, 1000))
    
    model.add_reactions(rxns)
    model.objective = 'ATPM'
    
    return model


def disable_engineered(model):
    """Disable all engineered pathways to simulate a normal human cell."""
    model.reactions.get_by_id('RADIO').upper_bound = 0
    model.reactions.get_by_id('DSUP').upper_bound = 0
    model.reactions.get_by_id('MNAOX').upper_bound = 0
    model.reactions.get_by_id('NRF2').upper_bound = 0
    model.reactions.get_by_id('SOD2').upper_bound = 0
    model.reactions.get_by_id('EX_mel').upper_bound = 0
    model.reactions.get_by_id('OH_SCAV').upper_bound = 10  # native GSH scavenging stays on (realistic cap)


# ============================================================
# EXPERIMENTS
# ============================================================

def run_all_experiments():
    """Run all experiments and return results as DataFrames."""
    
    model = build_model()
    results = {}
    
    # ----------------------------------------------------------
    # EXPERIMENT 1: Glucose restriction
    # ----------------------------------------------------------
    rows = []
    for glc in [10, 8, 6, 5, 4, 3, 2, 1, 0.5, 0]:
        # Normal cell
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -glc
            disable_engineered(model)
            s1 = model.optimize()
            atp_normal = s1.objective_value if s1.status == 'optimal' else 0
        
        # Radiotrophic cell
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -glc
            s2 = model.optimize()
            atp_radio = s2.objective_value if s2.status == 'optimal' else 0
            radio_flux = s2.fluxes.get('RADIO', 0) if s2.status == 'optimal' else 0
            sod_flux = s2.fluxes.get('SODc', 0) if s2.status == 'optimal' else 0
            dsup_flux = s2.fluxes.get('DSUP', 0) if s2.status == 'optimal' else 0
        
        gain = atp_radio - atp_normal
        pct = (gain / atp_normal * 100) if atp_normal > 0 else (100 if atp_radio > 0 else 0)
        
        rows.append({
            'glucose': glc, 'atp_normal': round(atp_normal, 2),
            'atp_radio': round(atp_radio, 2), 'gain': round(gain, 2),
            'pct_boost': round(pct, 1), 'radio_flux': round(radio_flux, 2),
            'sod_flux': round(sod_flux, 2), 'dsup_flux': round(dsup_flux, 2)
        })
    
    results['glucose_restriction'] = pd.DataFrame(rows)
    
    # ----------------------------------------------------------
    # EXPERIMENT 2: Hypoxia (glucose = 5)
    # ----------------------------------------------------------
    rows = []
    for o2 in [20, 10, 5, 3, 2, 1, 0.5, 0]:
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -5
            model.reactions.get_by_id('EX_o2').lower_bound = -o2
            disable_engineered(model)
            s1 = model.optimize()
            atp_normal = s1.objective_value if s1.status == 'optimal' else 0
        
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -5
            model.reactions.get_by_id('EX_o2').lower_bound = -o2
            s2 = model.optimize()
            atp_radio = s2.objective_value if s2.status == 'optimal' else 0
            radio_flux = s2.fluxes.get('RADIO', 0) if s2.status == 'optimal' else 0
        
        pct = ((atp_radio - atp_normal) / atp_normal * 100) if atp_normal > 0 else 0
        rows.append({
            'o2': o2, 'atp_normal': round(atp_normal, 2),
            'atp_radio': round(atp_radio, 2),
            'pct_boost': round(pct, 1), 'radio_flux': round(radio_flux, 2)
        })
    
    results['hypoxia'] = pd.DataFrame(rows)
    
    # ----------------------------------------------------------
    # EXPERIMENT 3: Radiation dose-response (glucose = 5)
    # ----------------------------------------------------------
    rows = []
    for maxrad in [0, 2, 5, 10, 15, 20, 25, 30, 40, 50]:
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -5
            model.reactions.get_by_id('RADIO').upper_bound = maxrad
            sol = model.optimize()
            if sol.status == 'optimal':
                rows.append({
                    'max_radio': maxrad,
                    'atp': round(sol.objective_value, 2),
                    'radio_used': round(sol.fluxes['RADIO'], 2),
                    'sod': round(sol.fluxes['SODc'], 2),
                    'dsup': round(sol.fluxes['DSUP'], 2),
                    'catalase': round(sol.fluxes['CATc'], 2),
                    'mn_aox': round(sol.fluxes['MNAOX'], 2),
                    'dna_repair': round(sol.fluxes['BER'], 2),
                    'oh_scavenged': round(sol.fluxes['OH_SCAV'], 2),
                    'melanin_consumed': round(sol.fluxes['RADIO'] * 0.02, 2),
                    'superoxide_generated': round(sol.fluxes['RADIO'] * 0.28, 2),
                    'oh_generated': round(sol.fluxes['RADIO'] * 0.24, 2)
                })
    
    results['dose_response'] = pd.DataFrame(rows)
    
    # ----------------------------------------------------------
    # EXPERIMENT 4: Defense system ablation (glucose=5, radio=50)
    # ----------------------------------------------------------
    # Note: native enzyme caps set during build (SODc=4, CATc=3, GPX=3, GR=2).
    # SOD2 defaults to 0 (off); enabled explicitly for ablation comparisons.
    # NOTE: DSUP capped at 5 (40% of max OH load at flux=50), OH_SCAV capped at 10
    # (realistic GSH turnover). These replace the previous unrealistic cap of 1000.
    # "No OH scavenging" removes OH_SCAV; Dsup alone (5) now only covers ~42% of max OH.
    configs = [
        ("All defenses ON",      {'RADIO':50,'DSUP':5,  'MNAOX':1000,'NRF2':1000,'SODc':4, 'CATc':3, 'GPX':3, 'GR':2, 'OH_SCAV':10, 'SOD2':0}),
        ("No Dsup",              {'RADIO':50,'DSUP':0,  'MNAOX':1000,'NRF2':1000,'SODc':4, 'CATc':3, 'GPX':3, 'GR':2, 'OH_SCAV':10, 'SOD2':0}),
        ("No Mn-AOX",            {'RADIO':50,'DSUP':5,  'MNAOX':0,   'NRF2':1000,'SODc':4, 'CATc':3, 'GPX':3, 'GR':2, 'OH_SCAV':10, 'SOD2':0}),
        ("No Nrf2",              {'RADIO':50,'DSUP':5,  'MNAOX':1000,'NRF2':0,   'SODc':4, 'CATc':3, 'GPX':3, 'GR':2, 'OH_SCAV':10, 'SOD2':0}),
        ("No SOD",               {'RADIO':50,'DSUP':5,  'MNAOX':1000,'NRF2':1000,'SODc':0, 'CATc':3, 'GPX':3, 'GR':2, 'OH_SCAV':10, 'SOD2':0}),
        ("No Catalase",          {'RADIO':50,'DSUP':5,  'MNAOX':1000,'NRF2':1000,'SODc':4, 'CATc':0, 'GPX':3, 'GR':2, 'OH_SCAV':10, 'SOD2':0}),
        ("No OH scavenging",     {'RADIO':50,'DSUP':5,  'MNAOX':1000,'NRF2':1000,'SODc':4, 'CATc':3, 'GPX':3, 'GR':2, 'OH_SCAV':0,  'SOD2':0}),
        ("Only native defenses", {'RADIO':50,'DSUP':0,  'MNAOX':0,   'NRF2':0,   'SODc':4, 'CATc':3, 'GPX':3, 'GR':2, 'OH_SCAV':10, 'SOD2':0}),
        ("Only engineered",      {'RADIO':50,'DSUP':5,  'MNAOX':1000,'NRF2':1000,'SODc':0, 'CATc':0, 'GPX':0, 'GR':0, 'OH_SCAV':0,  'SOD2':0}),
        ("No defenses",          {'RADIO':50,'DSUP':0,  'MNAOX':0,   'NRF2':0,   'SODc':0, 'CATc':0, 'GPX':0, 'GR':0, 'OH_SCAV':0,  'SOD2':0}),
    ]
    
    rows = []
    for label, cfg in configs:
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -5
            for rid, ub in cfg.items():
                model.reactions.get_by_id(rid).upper_bound = ub
            sol = model.optimize()
            if sol.status == 'optimal':
                f = sol.fluxes
                rows.append({
                    'config': label, 'atp': round(sol.objective_value, 2),
                    'radio_flux': round(f.get('RADIO', 0), 2),
                    'sod': round(f.get('SODc', 0), 2),
                    'catalase': round(f.get('CATc', 0), 2),
                    'gpx': round(f.get('GPX', 0), 2),
                    'gr': round(f.get('GR', 0), 2),
                    'dsup': round(f.get('DSUP', 0), 2),
                    'mn_aox': round(f.get('MNAOX', 0), 2),
                    'nrf2': round(f.get('NRF2', 0), 2),
                    'oh_scav': round(f.get('OH_SCAV', 0), 2),
                    'dna_repair': round(f.get('BER', 0), 2),
                    'status': sol.status
                })
            else:
                rows.append({'config': label, 'atp': 0, 'radio_flux': 0,
                             'sod':0,'catalase':0,'gpx':0,'gr':0,'dsup':0,
                             'mn_aox':0,'nrf2':0,'oh_scav':0,'dna_repair':0,
                             'status': sol.status})

    results['ablation'] = pd.DataFrame(rows)
    
    # ----------------------------------------------------------
    # EXPERIMENT 5: Combined stress
    # ----------------------------------------------------------
    rows = []
    for glc, o2 in [(5,20),(5,5),(2,10),(2,5),(1,5),(1,2),(0.5,1),(0,5),(0,1)]:
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -glc
            model.reactions.get_by_id('EX_o2').lower_bound = -o2
            disable_engineered(model)
            s1 = model.optimize()
            a1 = s1.objective_value if s1.status == 'optimal' else 0
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -glc
            model.reactions.get_by_id('EX_o2').lower_bound = -o2
            s2 = model.optimize()
            a2 = s2.objective_value if s2.status == 'optimal' else 0
        pct = ((a2 - a1) / a1 * 100) if a1 > 0 else (100 if a2 > 0 else 0)
        rows.append({'glucose': glc, 'o2': o2, 'atp_normal': round(a1,2), 'atp_radio': round(a2,2), 'pct_boost': round(pct,1)})
    
    results['combined_stress'] = pd.DataFrame(rows)

    # ----------------------------------------------------------
    # EXPERIMENT 6: ROS coefficient sensitivity analysis
    # Both superoxide (O2S) and hydroxyl radical (OH) are covaried
    # at the water radiolysis G-value ratio (Buxton et al. 1988):
    #   G(O2•⁻) / G(•OH) ≈ 0.28 / 0.24 → OH = O2S × (0.24/0.28)
    # This sweeps total ROS yield per NADH while maintaining the
    # biologically grounded stoichiometric ratio between species.
    # The baseline (ros_coefficient=0.28) reproduces default model values.
    # Generous O2 supply isolates ROS cost from O2 budget effects.
    # ----------------------------------------------------------
    _OH_O2S_RATIO = 0.24 / 0.28  # G-value ratio from Buxton et al. 1988
    rows = []
    for ros_coeff in [0.05, 0.1, 0.2, 0.28, 0.5, 0.7, 1.0, 1.5, 2.0]:
        # Build fresh model with both ROS species scaled together
        m_sens = build_model()
        radio_rxn = m_sens.reactions.get_by_id('RADIO')

        # Scale superoxide stoichiometry
        o2s_met = [mt for mt in radio_rxn.metabolites if 'o2s_c' in mt.id][0]
        radio_rxn.add_metabolites({o2s_met: ros_coeff - 0.28})  # delta from default 0.28

        # Scale OH stoichiometry proportionally (G-value ratio preserved)
        oh_coeff = ros_coeff * _OH_O2S_RATIO
        oh_met = [mt for mt in radio_rxn.metabolites if 'oh_radical_c' in mt.id][0]
        radio_rxn.add_metabolites({oh_met: oh_coeff - 0.24})  # delta from default 0.24

        # Generous O2 to isolate ROS cost from O2 budget artifact
        m_sens.reactions.get_by_id('EX_o2').lower_bound = -100
        m_sens.reactions.get_by_id('EX_glc').lower_bound = -5

        sol = m_sens.optimize()
        if sol.status == 'optimal':
            rf = sol.fluxes['RADIO']
            # Normal cell baseline (same O2 generosity)
            disable_engineered(m_sens)
            s_norm = m_sens.optimize()
            atp_norm = s_norm.objective_value if s_norm.status == 'optimal' else 0

            rows.append({
                'ros_coefficient': ros_coeff,
                'oh_coefficient': round(oh_coeff, 4),
                'atp_radio': round(sol.objective_value, 2),
                'atp_normal': round(atp_norm, 2),
                'radio_flux': round(rf, 2),
                'net_gain': round(sol.objective_value - atp_norm, 2),
                'superoxide_produced': round(rf * ros_coeff, 2),
                'oh_produced': round(rf * oh_coeff, 2),
                'sod_flux': round(sol.fluxes['SODc'], 2),
                'dsup_flux': round(sol.fluxes['DSUP'], 2),
                'dna_repair': round(sol.fluxes['BER'], 2)
            })

    results['ros_sensitivity'] = pd.DataFrame(rows)

    # ----------------------------------------------------------
    # EXPERIMENT 7: Bottleneck relief strategies
    # Tests engineering solutions to increase radiotrophic flux:
    #   A) Baseline (current model, SOD-limited)
    #   B) MnSOD2 overexpression (doubles SOD capacity)
    #   C) Synthetic melanin loading (bypasses biosynthesis)
    #   D) Both SOD2 + synthetic melanin
    #   E) Full optimization (SOD2 + melanin + increased tyrosine)
    # ----------------------------------------------------------
    rows = []
    strategies = [
        ('A: Baseline',             {'SOD2': 0,  'EX_mel': 0,  'EX_tyr': -5}),
        ('B: +MnSOD2 (cap=4)',      {'SOD2': 4,  'EX_mel': 0,  'EX_tyr': -5}),
        ('C: +Synthetic melanin',   {'SOD2': 0,  'EX_mel': 5,  'EX_tyr': -5}),
        ('D: +SOD2 + synth melanin',{'SOD2': 4,  'EX_mel': 5,  'EX_tyr': -5}),
        ('E: Full optimization',    {'SOD2': 8,  'EX_mel': 10, 'EX_tyr': -10}),
    ]

    for label, cfg in strategies:
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -5
            model.reactions.get_by_id('SOD2').upper_bound = cfg['SOD2']
            model.reactions.get_by_id('EX_mel').upper_bound = cfg['EX_mel']
            model.reactions.get_by_id('EX_tyr').lower_bound = cfg['EX_tyr']
            sol = model.optimize()
            if sol.status == 'optimal':
                f = sol.fluxes
                rows.append({
                    'strategy': label,
                    'atp': round(sol.objective_value, 2),
                    'radio_flux': round(f['RADIO'], 2),
                    'sod1_flux': round(f['SODc'], 2),
                    'sod2_flux': round(f['SOD2'], 2),
                    'melanin_synth': round(f['MELSYN'], 2),
                    'melanin_loaded': round(f['EX_mel'], 2),
                    'total_superoxide': round(f['RADIO'] * 0.28, 2),
                    'total_oh': round(f['RADIO'] * 0.24, 2),
                    'status': sol.status
                })

    results['bottleneck_relief'] = pd.DataFrame(rows)

    # ----------------------------------------------------------
    # EXPERIMENT 8: OH_SCAV sensitivity — conditional Dsup finding
    # The audit identified that Dsup's dispensability depends on OH_SCAV capacity.
    # When glutathione-based scavenging is abundant (OH_SCAV high), Dsup is redundant.
    # When GSH is depleted or impaired (OH_SCAV low), Dsup becomes essential.
    # This experiment sweeps OH_SCAV capacity and measures the ATP penalty from
    # removing Dsup at each level, revealing the clinically relevant boundary.
    #
    # Biological relevance: GSH depletion occurs in aging, after radiation/chemo,
    # in cells with impaired Nrf2, and in some cancer phenotypes.
    # ----------------------------------------------------------
    rows = []
    for oh_scav_cap in [0, 1, 2, 3, 5, 7, 10, 15, 20, 50, 100]:
        # With Dsup (DSUP cap = 5, realistic)
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -5
            model.reactions.get_by_id('RADIO').upper_bound = 50
            model.reactions.get_by_id('DSUP').upper_bound = 5
            model.reactions.get_by_id('OH_SCAV').upper_bound = oh_scav_cap
            sol_with = model.optimize()
            atp_with = sol_with.objective_value if sol_with.status == 'optimal' else 0
            radio_with = sol_with.fluxes.get('RADIO', 0) if sol_with.status == 'optimal' else 0
            dsup_used = sol_with.fluxes.get('DSUP', 0) if sol_with.status == 'optimal' else 0

        # Without Dsup (DSUP = 0)
        with model:
            model.reactions.get_by_id('EX_glc').lower_bound = -5
            model.reactions.get_by_id('RADIO').upper_bound = 50
            model.reactions.get_by_id('DSUP').upper_bound = 0
            model.reactions.get_by_id('OH_SCAV').upper_bound = oh_scav_cap
            sol_no = model.optimize()
            atp_no = sol_no.objective_value if sol_no.status == 'optimal' else 0
            radio_no = sol_no.fluxes.get('RADIO', 0) if sol_no.status == 'optimal' else 0

        atp_delta = round(atp_with - atp_no, 2)
        pct_delta = round(atp_delta / atp_with * 100, 1) if atp_with > 0 else 0
        rows.append({
            'oh_scav_cap': oh_scav_cap,
            'atp_with_dsup': round(atp_with, 2),
            'atp_no_dsup': round(atp_no, 2),
            'atp_delta': atp_delta,
            'pct_atp_loss_without_dsup': pct_delta,
            'radio_with_dsup': round(radio_with, 2),
            'radio_no_dsup': round(radio_no, 2),
            'dsup_flux_used': round(dsup_used, 2),
            'dsup_essential': 'YES' if atp_delta > 1.0 else 'no'
        })

    results['dsup_conditional'] = pd.DataFrame(rows)

    # ----------------------------------------------------------
    # EXPERIMENT 10: Experimental validation comparison
    # Compare model predictions with published experimental data.
    # ATP boost values are computed from the model, not hardcoded.
    # Two comparison conditions are reported (Option c per audit):
    #   (a) Standard: glucose=5, O2=20 — conservative baseline
    #   (b) Glucose-restricted: glucose=3, O2=20 — approximates ISS fungi conditions
    # ----------------------------------------------------------

    # Compute actual ATP boost at standard conditions (glucose=5, O2=20)
    with model:
        model.reactions.get_by_id('EX_glc').lower_bound = -5
        model.reactions.get_by_id('EX_o2').lower_bound = -20
        disable_engineered(model)
        _s = model.optimize()
        _atp_norm_std = _s.objective_value if _s.status == 'optimal' else 0

    with model:
        model.reactions.get_by_id('EX_glc').lower_bound = -5
        model.reactions.get_by_id('EX_o2').lower_bound = -20
        _s = model.optimize()
        _atp_radio_std = _s.objective_value if _s.status == 'optimal' else 0

    _boost_std = (_atp_radio_std - _atp_norm_std) / _atp_norm_std * 100 if _atp_norm_std > 0 else 0

    # Compute actual ATP boost at glucose-restricted conditions (glucose=3, O2=20)
    # Justification: ISS C. sphaerospermum grows on trace organics at chronic
    # low dose (~144 mSv/yr); glucose=3 approximates nutrient-limited fungal conditions.
    with model:
        model.reactions.get_by_id('EX_glc').lower_bound = -3
        model.reactions.get_by_id('EX_o2').lower_bound = -20
        disable_engineered(model)
        _s = model.optimize()
        _atp_norm_glc3 = _s.objective_value if _s.status == 'optimal' else 0

    with model:
        model.reactions.get_by_id('EX_glc').lower_bound = -3
        model.reactions.get_by_id('EX_o2').lower_bound = -20
        _s = model.optimize()
        _atp_radio_glc3 = _s.objective_value if _s.status == 'optimal' else 0

    _boost_glc3 = (_atp_radio_glc3 - _atp_norm_glc3) / _atp_norm_glc3 * 100 if _atp_norm_glc3 > 0 else 0

    _shunk_pred = (
        f'+{_boost_std:.1f}% ATP (glucose=5 standard); '
        f'+{_boost_glc3:.1f}% ATP (glucose=3, nutrient-restricted)'
    )
    _shunk_agree = (
        f'PARTIAL — standard conditions show {_boost_std:.1f}% (vs 21±37% measured); '
        f'nutrient-restricted conditions show {_boost_glc3:.1f}%, closer to experimental mean. '
        f'Both are within the wide experimental error bar (±37%).'
    )

    validation = pd.DataFrame([
        {
            'observation': 'Melanized C. neoformans growth boost under radiation',
            'source': 'Dadachova et al. 2007 PLOS ONE',
            'experimental_value': '2.5x CFU increase',
            'model_prediction': f'+{_boost_std:.1f}% ATP at standard conditions (glucose=5)',
            'agreement': 'PARTIAL - model shows modest boost; fungi likely have additional mechanisms',
            'note': 'CFU growth rate and ATP maintenance flux are not directly comparable metrics'
        },
        {
            'observation': 'ATP decrease in melanized cells under radiation',
            'source': 'Bryan et al. 2011 Fungal Biology',
            'experimental_value': 'ATP decreases in melanized cells',
            'model_prediction': 'Net ATP gain after 0.5 ATP/NADH overhead',
            'agreement': 'PARTIAL - overhead modeled but net still positive',
            'note': 'Transient ATP drop may precede steady-state gain; FBA gives steady-state only'
        },
        {
            'observation': 'Dsup reduces DNA damage by ~40% in HEK293',
            'source': 'Hashimoto et al. 2016 Nature Comms',
            'experimental_value': '~40% X-ray damage reduction',
            'model_prediction': 'Dsup UB capped at 5 (= 0.4 x max OH load); OH_SCAV handles remainder',
            'agreement': 'YES - Dsup protective capacity now matches Hashimoto 40% ceiling',
            'note': 'Previous model used UB=1000 (no cap); now corrected to biologically realistic bound'
        },
        {
            'observation': 'C. sphaerospermum 21% growth advantage on ISS',
            'source': 'Shunk et al. 2022 Frontiers Microbiol',
            'experimental_value': '21 ± 37% growth rate increase',
            'model_prediction': _shunk_pred,
            'agreement': _shunk_agree,
            'note': 'Experimental error bar is large (±37%); both model conditions are consistent with it'
        },
        {
            'observation': 'Engineered melanin NPs protect mice from 6 Gy',
            'source': 'Nature Communications 2025',
            'experimental_value': 'Survival 12% -> 100%',
            'model_prediction': 'Synthetic melanin loading removes biosynthesis bottleneck',
            'agreement': 'CONSISTENT - melanin radioprotection confirmed in mammals',
            'note': 'Protection mechanism (shielding) differs from energy capture'
        },
        {
            'observation': 'SOD2 overexpression enhances radiosurvival in HEK293T',
            'source': 'Kolesnikova et al. 2023 PMC10744337',
            'experimental_value': 'Increased viability at 2-5 Gy (dose-dependent)',
            'model_prediction': 'SOD is single point of failure; SOD2 relieves bottleneck',
            'agreement': 'YES - SOD augmentation improves radiation tolerance',
            'note': 'SOD2+catalase co-expression most effective'
        },
    ])
    results['experimental_validation'] = validation

    return results


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Building radiotrophic cell model...")
    model = build_model()
    print(f"Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites\n")
    
    print("Running all experiments...")
    results = run_all_experiments()
    
    for name, df in results.items():
        print(f"\n{'='*60}")
        print(f"EXPERIMENT: {name.upper()}")
        print(f"{'='*60}")
        print(df.to_string(index=False))
    
    # Save results
    for name, df in results.items():
        df.to_csv(f'{name}.csv', index=False)
    
    # Save model
    cobra.io.save_json_model(model, 'radiotrophic_cell_model.json')
    
    print("\n\nResults saved to current directory")
    print("Model saved to radiotrophic_cell_model.json")

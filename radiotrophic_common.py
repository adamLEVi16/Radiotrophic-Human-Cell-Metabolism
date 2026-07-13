"""
Shared physical basis for the radiotrophic cell models
======================================================
Single source of truth that connects the abstract model fluxes to real,
dimensioned quantities: absorbed radiation dose, radiolytic ROS production,
and the energy budget. Both the constraint-based (FBA) and kinetic (ODE)
models import from here so their results are expressed in the same currency.

Why this module exists
----------------------
The original FBA and ODE models used incompatible, largely arbitrary flux
units, and the radiotrophic reaction generated reducing equivalents without
any link to the energy actually deposited by radiation. That decoupling is
what let a "free" NADH source inflate the ATP yield. Everything here is
tied to absorbed dose (Gy = J/kg) so that claims can be checked against
conservation of energy.

Key references
--------------
  - Water radiolysis primary yields (G-values): Buxton et al. (1988)
    J. Phys. Chem. Ref. Data 17:513. Low-LET, neutral, ~microsecond yields.
  - In-vivo ATP hydrolysis free energy ~50-60 kJ/mol: Nicholls & Ferguson,
    Bioenergetics 4 (2013).
  - Mammalian cell mass ~3 ng: Park et al. (2008) PNAS (cell buoyant mass).
  - Whole-body ATP turnover ~ body mass/day: standard bioenergetics estimate.
"""

# ============================================================
# FUNDAMENTAL CONSTANTS
# ============================================================
AVOGADRO = 6.02214076e23           # 1/mol
EV_J = 1.602176634e-19             # J per eV

# 1 molecule / 100 eV  ->  mol / J
#   = 1 / (100 * EV_J) molecules/J / AVOGADRO
MOLEC_PER_100EV_TO_MOL_PER_J = 1.0 / (100.0 * EV_J) / AVOGADRO   # ~1.036e-7 mol/J

# ============================================================
# WATER RADIOLYSIS G-VALUES (Buxton et al. 1988)
# Primary radical yields per 100 eV of absorbed energy, low-LET water.
# In aerated cytosol the hydrated electron and H atom are scavenged by O2:
#     e-_aq + O2 -> O2*-        H* + O2 -> HO2* (<-> O2*- + H+)
# so the effective superoxide yield ~ G(e-_aq) + G(H*).
# ============================================================
G_OH_100EV = 2.7        # hydroxyl radical
G_EAQ_100EV = 2.6       # hydrated electron
G_H_100EV = 0.6         # hydrogen atom
G_H2O2_100EV = 0.7      # hydrogen peroxide
G_O2S_100EV = G_EAQ_100EV + G_H_100EV   # ~3.2, superoxide via O2 scavenging

# Converted to mol/J (SI). At rho_water = 1 kg/L these equal mol/(L*J) too,
# so dose-rate [Gy/s] * G [mol/J] gives production in mol/(L*s) = M/s.
G_OH = G_OH_100EV * MOLEC_PER_100EV_TO_MOL_PER_J          # ~0.28 umol/J
G_O2S = G_O2S_100EV * MOLEC_PER_100EV_TO_MOL_PER_J        # ~0.33 umol/J
G_H2O2 = G_H2O2_100EV * MOLEC_PER_100EV_TO_MOL_PER_J      # ~0.07 umol/J

# Fraction of primary radicals intercepted by melanin before they reach the
# bulk cytosol (Schweitzer et al. 2009 estimate melanin quenches ~half).
MELANIN_QUENCH = 0.5

# ============================================================
# CELLULAR ENERGETICS
# ============================================================
CELL_MASS_KG = 3e-12        # ~3 ng wet mass, typical mammalian cell (Park 2008)
WATER_DENSITY = 1.0         # kg/L (cytosol ~ water for dose purposes)

DG_ATP = 50e3               # J/mol, in-vivo ATP hydrolysis free energy (conservative)
DG_NADH = 220e3             # J/mol, reducing-equivalent free energy stored in NADH
                            # (NAD+ + H+ + 2e- -> NADH, vs O2). Used as the minimum
                            # radiation energy required to drive one NADH radiotrophically.

# ATP demand of a single mammalian cell, molecules/s.
# Whole-body average ~1e7/s; metabolically active cells up to ~1e9/s.
ATP_DEMAND_LOW = 1e7
ATP_DEMAND_HIGH = 1e9

# Dry weight is ~30% of wet mass, so 1 gDW ~ 3.3 g wet ~ 3.3e-3 kg tissue
# (dose is defined per kg of water-equivalent tissue).
WET_PER_GDW_KG = 3.3e-3

# ============================================================
# REFERENCE DOSE-RATE REGIMES (Gy/s = J/kg/s)
# ============================================================
DOSE_REGIMES = {
    "natural_background": 7.6e-11,   # ~2.4 mSv/yr
    "iss_leo":            4.6e-9,     # ~144 mSv/yr, low Earth orbit
    "chernobyl_hotspot":  2.8e-7,    # ~1 mSv/hr contaminated hotspot
    "lab_gamma_source":   1e-1,      # ~0.1 Gy/s radiobiology irradiator
    "acutely_lethal_1s":  5.0,       # ~5 Gy in 1 s (whole-body lethal dose)
}

# ============================================================
# DOSE  <->  RADIOLYTIC ROS
# ============================================================
def ros_production_rates(dose_rate, melanin_quench=MELANIN_QUENCH):
    """
    Radiolytic ROS production rates at a given absorbed dose rate.

    Args:
        dose_rate: absorbed dose rate in Gy/s (= J/kg/s).
        melanin_quench: fraction of primary radicals removed by melanin.

    Returns:
        dict of production rates in M/s (mol/L/s): 'oh', 'o2s', 'h2o2'.
    """
    keep = 1.0 - melanin_quench
    return {
        "oh":   dose_rate * G_OH * keep,
        "o2s":  dose_rate * G_O2S * keep,
        "h2o2": dose_rate * G_H2O2 * keep,
    }

# ============================================================
# DOSE  <->  ENERGY BUDGET
# ============================================================
def radiation_power_per_cell(dose_rate):
    """Absorbed radiation power for one cell, in watts."""
    return dose_rate * CELL_MASS_KG

def atp_demand_power(atp_per_s):
    """Power required to sustain a given ATP turnover rate, in watts."""
    return atp_per_s / AVOGADRO * DG_ATP

def energy_fraction(dose_rate, atp_per_s, efficiency=1.0):
    """
    Fraction of a cell's ATP-demand power that absorbed radiation could
    supply at the given energy-conversion efficiency (1.0 = thermodynamic
    ceiling). Values << 1 mean radiotrophy is energetically negligible.
    """
    return efficiency * radiation_power_per_cell(dose_rate) / atp_demand_power(atp_per_s)

def breakeven_dose_rate(atp_per_s, efficiency=1.0):
    """Dose rate (Gy/s) at which radiation power equals ATP-demand power."""
    return atp_demand_power(atp_per_s) / (efficiency * CELL_MASS_KG)

# ============================================================
# DOSE  <->  FBA RADIO FLUX
# The FBA model expresses fluxes in mmol/gDW/h. Radiotrophic NADH generation
# is bounded by the radiation energy absorbed per gDW per hour divided by the
# energy stored per NADH. This is the constraint the original model lacked.
# ============================================================
def max_radio_flux(dose_rate, efficiency=1.0):
    """
    Upper bound on the RADIO reaction (mmol NADH / gDW / h) permitted by
    energy conservation at a given dose rate.
    """
    energy_per_gdw_per_h = dose_rate * WET_PER_GDW_KG * 3600.0   # J/gDW/h
    mol_nadh_per_h = efficiency * energy_per_gdw_per_h / DG_NADH  # mol/gDW/h
    return mol_nadh_per_h * 1e3                                   # mmol/gDW/h

def dose_for_radio_flux(radio_flux, efficiency=1.0):
    """Inverse of max_radio_flux: dose rate needed to sustain a RADIO flux."""
    mol_nadh_per_h = radio_flux * 1e-3
    energy_per_gdw_per_h = mol_nadh_per_h * DG_NADH / efficiency
    return energy_per_gdw_per_h / (WET_PER_GDW_KG * 3600.0)


if __name__ == "__main__":
    print("Radiotrophic model — shared physical basis")
    print("=" * 55)
    print(f"G(*OH)  = {G_OH*1e6:.3f} umol/J   G(O2*-) = {G_O2S*1e6:.3f} umol/J")
    print(f"Cell mass = {CELL_MASS_KG*1e12:.1f} ng   dG_ATP = {DG_ATP/1e3:.0f} kJ/mol")
    print()
    print("Dose rate needed to reach the original model's RADIO flux (~28.6):")
    for eff in (1.0, 0.1, 0.01):
        print(f"  efficiency {eff:5.2f} -> {dose_for_radio_flux(28.57, eff):.1f} Gy/s")
    print()
    print("Break-even dose rate (radiation power = ATP demand):")
    for atp in (ATP_DEMAND_LOW, ATP_DEMAND_HIGH):
        print(f"  demand {atp:.0e} ATP/s -> {breakeven_dose_rate(atp):.2f} Gy/s")

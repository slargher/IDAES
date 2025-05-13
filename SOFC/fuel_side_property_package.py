import logging
from pyomo.environ import units as pyunits
from idaes.core import VaporPhase, LiquidPhase, Component
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.phase_equil import CubicComplementarityVLE
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.pure import RPP4
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import LogBubbleDew

# Set up logger
_log = logging.getLogger(__name__)

def get_fuel_property_package():
    # Configuration for fuel-side (NH₃, H₂, N₂, H₂O)
    configuration = {
        "components": {
            # Ammonia NH₃ component
            "nh3": {
                "type": Component,  # Define ammonia as a component
                "elemental_composition": {"H": 3, "N": 1},
                "enth_mol_ig_comp": RPP4,  # Ideal gas enthalpy function 
                "entr_mol_ig_comp": RPP4,  # Ideal gas entropy function 
                "pressure_sat_comp": RPP4,  # Saturation pressure function for ammonia
                "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},  # Defines the phase equilibrium equation for vapor and liquid phases
                "parameter_data": {
                    "mw": (17.0305e-3, pyunits.kg / pyunits.mol), #molecular weight data taken from NIST
                    "pressure_crit": (113.57e5, pyunits.Pa), #critical pressure data from engineering toolbox
                    "temperature_crit": (405.56, pyunits.K), #critical temperature data from engineering toolbox
                    "omega": 0.253, # Acentric factor of ammonia (a measure of non-ideal behavior)
                    "cp_mol_ig_comp_coeff": { #check these data and if it is based on the shomate equation # Coefficients for ideal gas heat capacity calculation using a polynomial function
                        "A": (19.99563, pyunits.J / pyunits.mol / pyunits.K),  # Cp° = A + B*t + C*t2 + D*t3 + E/t2
                        "B": (49.77119, pyunits.J / pyunits.mol / pyunits.K**2),  # Gas Phase Heat Capacity (Shomate Equation)
                        "C": (-15.37599, pyunits.J / pyunits.mol / pyunits.K**3),  
                        "D": (1.921168, pyunits.J / pyunits.mol / pyunits.K**4),  
                    },
                    "enth_mol_form_vap_comp_ref": (45.9e3, pyunits.J / pyunits.mol), # Enthalpy of formation of ammonia vapor (J/mol)
                    "entr_mol_form_vap_comp_ref": (-200, pyunits.J / pyunits.mol / pyunits.K), # Entropy of formation of ammonia vapor (J/mol·K)
                    "pressure_sat_comp_coeff": { # Antoine coefficients for calculating the saturation pressure of ammonia
                        "A": (7.36048, None),  
                        "B": (926.13, None),    # data from Physical and Chemical Equilibrium for Chemical Engineers, Second Edition. Noel de Nevers.
                        "C": (240.17, None),
                    },
                },
            },

            # Hydrogen H₂ component
            "h2": {
                "type": Component,
                "enth_mol_ig_comp": RPP4,
                "entr_mol_ig_comp": RPP4,
                "pressure_sat_comp": RPP4,
                "phase_equilibrium_form": {("Vap"): log_fugacity},
                "parameter_data": RPP4.H2,
            },
            # Nitrogen N₂ component
            "n2": {
                "type": Component,
                "enth_mol_ig_comp": RPP4,
                "entr_mol_ig_comp": RPP4,
                "pressure_sat_comp": RPP4,
                "phase_equilibrium_form": {("Vap"): log_fugacity},
                "parameter_data": RPP4.N2,
            },
            # Water H₂O component
            "h2o": {
                "type": Component,
                "enth_mol_ig_comp": RPP4,
                "entr_mol_ig_comp": RPP4,
                "pressure_sat_comp": RPP4,
                "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                "parameter_data": RPP4.H2O,
            },
        },
        "phases": {
            "Vap": {
                "type": VaporPhase,
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},
            },
        },
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 300, 500, pyunits.K),
            "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
        },
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
        "phases_in_equilibrium": [("Vap", "Liq")],
        "phase_equilibrium_state": {("Vap", "Liq"): CubicComplementarityVLE},
        "bubble_dew_method": LogBubbleDew,
        "parameter_data": {
            "PR_kappa": {
                ("nh3", "nh3"): 0.000,
            }
        },
    }

    return configuration

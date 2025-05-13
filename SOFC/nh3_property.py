# Import Python libraries
import logging
import copy
import enum

# Import Pyomo units
import pyomo.environ as pyo
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import VaporPhase, LiquidPhase, Component, PhaseType

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil.forms import (
    log_fugacity,
)
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.pure import (
    NIST,
    RPP4,
    RPP5,
    ChapmanEnskogLennardJones,
    Eucken,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    ConcentrationForm,
)
from idaes.models.properties.modular_properties.transport_properties import (
    ViscosityWilke,
    ThermalConductivityWMS,
    NoMethod,
)
from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import (
    wilke_phi_ij_callback,
)

from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.rate_constant import arrhenius
from idaes.models.properties.modular_properties.reactions.rate_forms import (
    power_law_rate,
)
from idaes.core.util.exceptions import ConfigurationError

# Set up logger
_log = logging.getLogger(__name__)


class EosType(enum.Enum):
    PR = 1
    IDEAL = 2


# Property Sources

# Source: NIST webbook
# Properties: Heat capacity coefficients for all species except ethane,
# propane, and butane. Reference enthalpies and entropies for all species.

# Source: The Properties of Gases and Liquids (1987)
# 4th edition, Chemical Engineering Series - Robert C. Reid
# Properties: Critical temperatures and pressures. Omega.
# Heat capacity coefficients for ethane, propane, and butane.

_phase_dicts_pr = {
    "Vap": {
        "type": VaporPhase,
        "equation_of_state": Cubic,
        "equation_of_state_options": {"type": CubicType.PR},
        "visc_d_phase": ViscosityWilke,
        "therm_cond_phase": ThermalConductivityWMS,
    },
    "Liq": {
        "type": LiquidPhase,
        "equation_of_state": Cubic,
        "equation_of_state_options": {"type": CubicType.PR},
        "visc_d_phase": NoMethod,
        "therm_cond_phase": NoMethod,
    },
}

_phase_dicts_ideal = {
    "Vap": {
        "type": VaporPhase,
        "equation_of_state": Ideal,
        "visc_d_phase": ViscosityWilke,
        "transport_property_options": {
            "viscosity_phi_ij_callback": wilke_phi_ij_callback,
        },
        "therm_cond_phase": ThermalConductivityWMS,
    },
}

_component_params = {
    "H2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"H": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "visc_d_phase_comp": {"Vap": ChapmanEnskogLennardJones},
        "therm_cond_phase_comp": {"Vap": Eucken},
        "parameter_data": {
            "mw": (0.0020159, pyunits.kg / pyunits.mol),
            "pressure_crit": (13e5, pyunits.Pa),
            "temperature_crit": (33.2, pyunits.K),
            "omega": -0.218,
            "cp_mol_ig_comp_coeff": {
                "A": 33.066178,
                "B": -11.363417,
                "C": 11.432816,
                "D": -2.772874,
                "E": -0.158558,
                "F": -9.980797,
                "G": 172.707974,
                "H": 0.0,
            },
            "lennard_jones_sigma": (2.826, pyunits.angstrom),
            "lennard_jones_epsilon_reduced": (59.7, pyunits.K),
            "f_int_eucken": 1,
        },
    },
        "CO": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"C": 1, "O": 1},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "visc_d_phase_comp": {"Vap": ChapmanEnskogLennardJones},
        "therm_cond_phase_comp": {"Vap": Eucken},
        "parameter_data": {
            "mw": (0.0280101, pyunits.kg / pyunits.mol),
            "pressure_crit": (35e5, pyunits.Pa),
            "temperature_crit": (132.9, pyunits.K),
            "omega": 0.066,
            "cp_mol_ig_comp_coeff": {
                "A": 25.56759,
                "B": 6.09613,
                "C": 4.054656,
                "D": -2.671301,
                "E": 0.131021,
                "F": -118.0089,
                "G": 227.3665,
                "H": -110.5271,
            },
            "lennard_jones_sigma": (3.690, pyunits.angstrom),
            "lennard_jones_epsilon_reduced": (91.7, pyunits.K),
            "f_int_eucken": 1,
        },
    },
    "H2O": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase, PhaseType.liquidPhase],
        "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
        "elemental_composition": {"H": 2, "O": 1},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "pressure_sat_comp": NIST,
        "parameter_data": {
            "mw": (0.01801528, pyunits.kg / pyunits.mol),
            "pressure_crit": (221.2e5, pyunits.Pa),
            "temperature_crit": (647.3, pyunits.K),
            "omega": 0.344,
            "cp_mol_ig_comp_coeff": {
                "A": 30.092,
                "B": 6.832514,
                "C": 6.793435,
                "D": -2.53448,
                "E": 0.082139,
                "F": -250.881,
                "G": 223.3967,
                "H": -241.8264,
            },
            "pressure_sat_comp_coeff": {  # NIST <- Stull 1947
                "A": 4.6543,
                "B": 1435.264,
                "C": -64.848,
            },
            "lennard_jones_sigma": (2.641, pyunits.angstrom),
            "lennard_jones_epsilon_reduced": (809.1, pyunits.K),
            "f_int_eucken": 1,
        },
    },
    "CO2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "visc_d_phase_comp": {"Vap": ChapmanEnskogLennardJones},
        "therm_cond_phase_comp": {"Vap": Eucken},
        "elemental_composition": {"C": 1, "O": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.04401, pyunits.kg / pyunits.mol),
            "pressure_crit": (73.8e5, pyunits.Pa),
            "temperature_crit": (304.1, pyunits.K),
            "omega": 0.239,
            "cp_mol_ig_comp_coeff": {
                "A": 24.99735,
                "B": 55.18696,
                "C": -33.69137,
                "D": 7.948387,
                "E": -0.136638,
                "F": -403.6075,
                "G": 228.2431,
                "H": -393.5224,
            },
            "lennard_jones_sigma": (3.941, pyunits.angstrom),
            "lennard_jones_epsilon_reduced": (195.2, pyunits.K),
            "f_int_eucken": 1,
        },
    },
    "O2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "visc_d_phase_comp": {"Vap": ChapmanEnskogLennardJones},
        "therm_cond_phase_comp": {"Vap": Eucken},
        "elemental_composition": {"O": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.031998, pyunits.kg / pyunits.mol),
            "pressure_crit": (50.4e5, pyunits.Pa),
            "temperature_crit": (154.6, pyunits.K),
            "omega": 0.025,
            "cp_mol_ig_comp_coeff": {
                "A": 30.03235,
                "B": 8.772972,
                "C": -3.988133,
                "D": 0.788313,
                "E": -0.741599,
                "F": -11.32468,
                "G": 236.1663,
                "H": 0.0,
            },
            "lennard_jones_sigma": (3.467, pyunits.angstrom),
            "lennard_jones_epsilon_reduced": (106.7, pyunits.K),
            "f_int_eucken": 1,
        },
    },
    "N2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"N": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "visc_d_phase_comp": {"Vap": ChapmanEnskogLennardJones},
        "therm_cond_phase_comp": {"Vap": Eucken},
        "parameter_data": {
            "mw": (0.0280134, pyunits.kg / pyunits.mol),
            "pressure_crit": (33.9e5, pyunits.Pa),
            "temperature_crit": (126.2, pyunits.K),
            "omega": 0.039,
            "cp_mol_ig_comp_coeff": {
                "A": 19.50583,
                "B": 19.88705,
                "C": -8.598535,
                "D": 1.369784,
                "E": 0.527601,
                "F": -4.935202,
                "G": 212.39,
                "H": 0.0,
            },
            "lennard_jones_sigma": (3.798, pyunits.angstrom),
            "lennard_jones_epsilon_reduced": (71.4, pyunits.K),
            "f_int_eucken": 1,
        },
    },
    }

_water_visc_d = {"Vap": ChapmanEnskogLennardJones, "Liq": None}
_water_therm_cond = {"Vap": Eucken, "Liq": None}
# returns a configuration dictionary for the list of specified components
def get_property(components=None, phases="Vap", eos=EosType.PR, scaled=False):
    if components is None:
        components = list(_component_params.keys())
    configuration = {
        "components": {},  # fill in later based on selected components
        "parameter_data": {},
        "phases": {},
        # Set base units of measurement
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        # Specifying state definition
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (0, 8000, 50000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 500, 2500, pyunits.K),
            "pressure": (5e4, 1.3e5, 1e8, pyunits.Pa),
        },
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
    }

    c = configuration["components"]
    if isinstance(phases, str):
        phases = [phases]
    for comp in components:
        c[comp] = copy.deepcopy(_component_params[comp])
        if comp == "H2O":
            c["H2O"]["visc_d_phase_comp"] = copy.deepcopy(
                {p: _water_visc_d[p] for p in phases}
            )
            c["H2O"]["therm_cond_phase_comp"] = copy.deepcopy(
                {p: _water_therm_cond[p] for p in phases}
            )
    for k in phases:
        if eos == EosType.PR:
            configuration["phases"][k] = copy.deepcopy(_phase_dicts_pr[k])
        elif eos == EosType.IDEAL:
            if k == "Liq":
                raise ConfigurationError(
                    "This parameter set does not support Ideal EOS with liquid"
                )
            configuration["phases"][k] = copy.deepcopy(_phase_dicts_ideal[k])
        else:
            raise ValueError("Invalid EoS.")
    if len(phases) > 1:
        p = tuple(phases)
        configuration["phases_in_equilibrium"] = [p]
        configuration["phase_equilibrium_state"] = {p: SmoothVLE}

    # Fill the binary parameters with zeros.
    d = configuration["parameter_data"]
    d["PR_kappa"] = {(a, b): 0 for a in c for b in c}

    # Change to scaled units if specified
    if scaled:
        configuration["base_units"]["mass"] = pyunits.Mg
        configuration["base_units"]["amount"] = pyunits.kmol

    return configuration

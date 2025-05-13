from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.phase_equil import Ideal
from idaes.models.properties.modular_properties.pure import RPP4
from idaes.core import VaporPhase, Component

def get_air_property_package():
    air_components = {
        "O2": {"type": Component, "valid_phase_types": VaporPhase, "parameter_data": RPP4.O2},
        "N2": {"type": Component, "valid_phase_types": VaporPhase, "parameter_data": RPP4.N2},
        "H2O": {"type": Component, "valid_phase_types": VaporPhase, "parameter_data": RPP4.H2O},
    }

    air_config = {
        "components": air_components,
        "phases": {"Vap": {"type": VaporPhase, "equation_of_state": Ideal}},
        "base_units": {
            "time": "s",
            "length": "m",
            "mass": "kg",
            "amount": "mol",
            "temperature": "K",
        },
    }

    return GenericParameterBlock(default=air_config)

# Save the updated script with DOF prints added
# Import coerenti con il tuo stile
from pyomo.environ import ConcreteModel, SolverFactory, value
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.models.unit_models import Translator, Separator, Mixer
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models_extra.power_generation.properties.natural_gas_PR import get_prop, get_rxn, EosType
import pyomo.environ as pyo
import idaes.core.util.scaling as iscale
import idaes.core.util.model_statistics as stattools
from idaes.core.util.initialization import propagate_state
import idaes.logger as idaeslog
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.models.unit_models as um
from pyomo.network import Arc
from idaes.models.properties.modular_properties.base.generic_reaction import GenericReactionParameterBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom, 
    report_statistics,
    activated_constraints_set,
    activated_equalities_set,
    unfixed_variables_set,
    fixed_variables_set,
    variables_set
)

# Logger silenzioso
idaeslog.getLogger("idaes.init").setLevel(idaeslog.ERROR)
idaeslog.getLogger("idaes.solve").setLevel(idaeslog.ERROR)

# Modello
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

# Proprietà
fuel_comps = {"H2": 0.99, "H2O": 0.01}
air_comps = {"O2": 0.21, "N2": 0.79}
reaction_comps = {"H2", "H2O", "O2"}

m.fs.fuel_props = GenericParameterBlock(**get_prop(fuel_comps, phases={"Vap"}, eos=EosType.PR))
m.fs.air_props = GenericParameterBlock(**get_prop(air_comps, phases={"Vap"}, eos=EosType.PR))
m.fs.reaction_props = GenericParameterBlock(**get_prop(reaction_comps, phases={"Vap"}, eos=EosType.PR))

def print_stream(unit, label, stream):
    separator = '=' * 60
    subsection = '-' * 60

    print(f"\n{separator}")
    print(f"{label:^60}")
    print(f"{separator}")

    print("\nSTREAM PROPERTIES:")
    print(f"{'Flow [mol/s]':<20}: {value(stream.flow_mol[0]):.6e}")
    print(f"{'Temperature [K]':<20}: {value(stream.temperature[0]):.2f}")
    print(f"{'Pressure [Pa]':<20}: {value(stream.pressure[0]):.2f}")

    print("\nMOLE FRACTIONS:")
    for t, c in stream.mole_frac_comp.keys():
        print(f"  x_{c:<10}: {value(stream.mole_frac_comp[t, c]):.6e}")

    print(f"\nDOF after {label}: {degrees_of_freedom(unit)}")

    print(f"\n{subsection}")
    print("MODEL STATISTICS REPORT:")
    print(subsection)
    report_statistics(unit)

    print(f"\n{subsection}")
    print("VARIABLES:")
    print(subsection)
    for v in variables_set(unit):
        print(f"  - {v}")

    print(f"\n{subsection}")
    print("ACTIVATED CONSTRAINTS:")
    print(subsection)
    for c in activated_constraints_set(unit):
        print(f"  - {c}")

    print(f"\n{subsection}")
    print("ACTIVATED EQUALITIES:")
    print(subsection)
    for eq in activated_equalities_set(unit):
        print(f"  - {eq}")

    print(f"\n{subsection}")
    print("UNFIXED VARIABLES:")
    print(subsection)
    for ufv in unfixed_variables_set(unit):
        print(f"  - {ufv}")

    print(f"\n{subsection}")
    print("FIXED VARIABLES:")
    print(subsection)
    for fv in fixed_variables_set(unit):
        print(f"  - {fv}")

    print(separator)
    print("\n")

# ----------------------------------------------------------------------------------------------------------------
# UNITÀ: Translator H2
# ----------------------------------------------------------------------------------------------------------------
m.fs.translator_h2 = Translator(
    inlet_property_package=m.fs.fuel_props,
    outlet_property_package=m.fs.reaction_props,
    outlet_state_defined=False,
)
m.fs.translator_h2.inlet.flow_mol[0].fix(0.00015)
m.fs.translator_h2.inlet.temperature[0].fix(1000)
m.fs.translator_h2.inlet.pressure[0].fix(101325)
m.fs.translator_h2.inlet.mole_frac_comp[:, "H2"].fix(0.99)
m.fs.translator_h2.inlet.mole_frac_comp[:, "H2O"].fix(0.01)

@m.fs.translator_h2.Constraint(m.fs.time)
def temp_eqn_h2(b, t): 
    return b.properties_out[t].temperature == b.properties_in[t].temperature

@m.fs.translator_h2.Constraint(m.fs.time)
def press_eqn_h2(b, t): 
    return b.properties_out[t].pressure == b.properties_in[t].pressure

@m.fs.translator_h2.Constraint(m.fs.time)
def flow_eqn_h2(b, t): 
    return b.properties_out[t].flow_mol == b.properties_in[t].flow_mol

@m.fs.translator_h2.Constraint(m.fs.time, list(m.fs.translator_h2.properties_in[0].mole_frac_comp.keys()))
def mole_frac_eqn_h2(b, t, j): 
    return b.properties_out[t].mole_frac_comp[j] == b.properties_in[t].mole_frac_comp[j]


# ----------------------------------------------------------------------------------------------------------------
# UNITÀ: Separator Air in o2_rich e o2_poor
# ----------------------------------------------------------------------------------------------------------------
m.fs.separator = Separator(
    property_package=m.fs.air_props,
    outlet_list=["o2_rich_strm", "o2_poor_strm"],
    split_basis=um.SplittingType.componentFlow,
)
m.fs.separator.inlet.flow_mol[0].fix(0.00015 * 10)
m.fs.separator.inlet.temperature[0].fix(1000)
m.fs.separator.inlet.pressure[0].fix(101325)
m.fs.separator.inlet.mole_frac_comp[0, "O2"].fix(0.21)
m.fs.separator.inlet.mole_frac_comp[0, "N2"].fix(0.79)
m.fs.separator.split_fraction[0, "o2_rich_strm", "O2"].fix(1.0)
m.fs.separator.split_fraction[0, "o2_rich_strm", "N2"].fix(0.0)


# ----------------------------------------------------------------------------------------------------------------
# UNITÀ: translator o2_rich da separator a reaction_props
# ----------------------------------------------------------------------------------------------------------------
m.fs.translator_o2 = Translator(
    inlet_property_package=m.fs.air_props,
    outlet_property_package=m.fs.reaction_props,
    outlet_state_defined=False,
)
m.fs.o2_to_translator = Arc(source=m.fs.separator.o2_rich_strm, destination=m.fs.translator_o2.inlet)

@m.fs.translator_o2.Constraint(m.fs.time)
def flow_eqn_o2(b, t): 
    return b.properties_out[t].flow_mol == b.properties_in[t].flow_mol

@m.fs.translator_o2.Constraint(m.fs.time)
def temp_eqn_o2(b, t): 
    return b.properties_out[t].temperature == b.properties_in[t].temperature

@m.fs.translator_o2.Constraint(m.fs.time)
def press_eqn_o2(b, t): 
    return b.properties_out[t].pressure == b.properties_in[t].pressure

m.fs.translator_o2.outlet.mole_frac_comp[:, "O2"].fix(1.0)


# ----------------------------------------------------------------------------------------------------------------
# UNITÀ: Mixer per H2 e O2 inlet
# ----------------------------------------------------------------------------------------------------------------
m.fs.mixer = um.Mixer(
    inlet_list=["fuel_strm","o2_rich_strm"],
    momentum_mixing_type=um.MomentumMixingType.none,
    property_package=m.fs.reaction_props,
)
m.fs.mix_h2 = Arc(source=m.fs.translator_h2.outlet, destination=m.fs.mixer.fuel_strm)
m.fs.mix_o2 = Arc(source=m.fs.translator_o2.outlet, destination=m.fs.mixer.o2_rich_strm)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

@m.fs.mixer.Constraint(m.fs.time)
def pressure_eqn(b, t):
    return b.mixed_state[t].pressure == b.o2_rich_strm_state[t].pressure


# ----------------------------------------------------------------------------------------------------------------
# UNITÀ: Stochiometric Reactor
# ----------------------------------------------------------------------------------------------------------------
m.fs.rxn_props = GenericReactionParameterBlock(**get_rxn(m.fs.reaction_props, {"h2_cmb"}))
m.fs.reactor = um.StoichiometricReactor(
    property_package=m.fs.reaction_props,
    reaction_package=m.fs.rxn_props,
    has_pressure_change=False,
    has_heat_transfer=True,
)

m.fs.mix_to_reactor = Arc(source=m.fs.mixer.outlet, destination=m.fs.reactor.inlet)
pyo.TransformationFactory("network.expand_arcs").apply_to(m)


# ----------------------------------------------------------------------------------------------------------------
# UNITÀ: Stochiometric Reactor
# ----------------------------------------------------------------------------------------------------------------
m.fs.reactor_separator = Separator(
    property_package=m.fs.reaction_props,
    outlet_list=['water_strm', 'o2_strm'],
    split_basis=um.SplittingType.componentFlow,
)
m.fs.react_to_sep = Arc(source=m.fs.reactor.outlet, destination=m.fs.reactor_separator.inlet)
pyo.TransformationFactory("network.expand_arcs").apply_to(m)

for j in ["H2", "H2O"]:
    m.fs.reactor_separator.split_fraction[0, "water_strm", j].fix(1.0)
m.fs.reactor_separator.split_fraction[0, "o2_strm", "O2"].fix(1.0)


# ----------------------------------------------------------------------------------------------------------------
# UNITÀ: Translator for H2O e H2 stream
# ----------------------------------------------------------------------------------------------------------------
m.fs.translator_h2_out = Translator(
    inlet_property_package=m.fs.reaction_props,
    outlet_property_package=m.fs.fuel_props,
    outlet_state_defined=False,
)

# Connessione dallo stream "water_strm" del separatore
m.fs.water_to_translator = Arc(
    source=m.fs.reactor_separator.water_strm,
    destination=m.fs.translator_h2_out.inlet
)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

# Constraint FTP standard (totale)
@m.fs.translator_h2_out.Constraint(m.fs.time)
def flow_eqn_h2_out(b, t):
    return b.properties_out[t].flow_mol == b.properties_in[t].flow_mol

@m.fs.translator_h2_out.Constraint(m.fs.time)
def temp_eqn_h2_out(b, t):
    return b.properties_out[t].temperature == b.properties_in[t].temperature

@m.fs.translator_h2_out.Constraint(m.fs.time)
def press_eqn_h2_out(b, t):
    return b.properties_out[t].pressure == b.properties_in[t].pressure

# Componenti comuni
comps_out = set(m.fs.translator_h2_out.properties_out[0].mole_frac_comp.keys())
comps_out.remove("H2")
# SOLO constraint flussi molari individuali per componenti comuni
@m.fs.translator_h2_out.Constraint(m.fs.time, comps_out)
def component_flow_eqn_h2_out(b, t, j):
    return (
                b.properties_out[t].mole_frac_comp[j]
                == b.properties_in[t].mole_frac_comp[j]
            )
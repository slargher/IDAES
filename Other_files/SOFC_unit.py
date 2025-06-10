import pyomo.environ as pyo
from pyomo.network import Arc
from pyomo.environ import Var
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
import idaes.core.util.constants as iconst

from idaes.core.util.config import is_physical_parameter_block
from idaes.core import UnitModelBlockData, declare_process_block_class
from idaes.core import MaterialBalanceType, MomentumBalanceType
import idaes.models.unit_models as um
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import GenericReactionParameterBlock
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop, get_rxn, EosType,
)
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom, report_statistics
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.logger as idaeslog

@declare_process_block_class("SofcUnit", doc="Simple SOFC model for process design.")
class SofcDesignData(UnitModelBlockData):

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    CONFIG.declare(
        "cathode_side_prop_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package for the oxygen side.",
            doc=(
                "Property package for the oxygen side, using "
                "idaes.models_extra.power_generation.properties.natural_gas_PR is "
                "strongly recommended, either Peng-Robinson or Ideal is okay"
            ),
        ),
    )
    CONFIG.declare(
        "anode_side_prop_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package for the hydrogen side.",
            doc=(
                "Property package for the hydrogen side, using "
                "idaes.models_extra.power_generation.properties.natural_gas_PR is "
                "strongly recommended, either Peng-Robinson or Ideal is okay"
            ),
        ),
    )
    CONFIG.declare(
        "cathode_side_prop_package_args",
        ConfigBlock(
            implicit=True,
            description="Property package arguments for the oxygen side.",
            doc="Property package arguments for the oxygen side.",
        ),
    )
    CONFIG.declare(
        "anode_side_prop_package_args",
        ConfigBlock(
            implicit=True,
            description="Property package arguments for the hydrogen side.",
            doc="Property package arguments for the hydrogen side.",
        ),
    )
    CONFIG.declare(
        "reaction_eos",
        ConfigValue(
            default=EosType.PR,
            domain=In(EosType),
            description="Physical properties for electrolysis reactions",
            doc=(
                "Reaction properties equation of state in: "
                "{EosType.PR, EosType.IDEAL}."
            ),
        ),
    )
    CONFIG.declare(
        "has_heat_transfer",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat transfer term construction flag",
            doc="Indicates whether the SOEC is adiabatic. Default=False (adiabatic)",
        ),
    )

    def build(self):
        """Construct a unit model data block"""
        super().build()

        if set(self.config.anode_side_prop_package.component_list) != {
            "H2",
            "H2O",
        }:
            raise ConfigurationError(
                "SOFC anode side must contain exactly H2 and H2O"
            )
        
        if set(self.config.cathode_side_prop_package.component_list) != {
            'O2',
            'N2',
        }:
            raise ConfigurationError(
                'SOFC cathode side must contain exactly O2 and N2'
            )

        self._add_reaction()
        self._add_unit_models()
        self._add_arcs()
        self._add_variables()
        self._add_constraints()
        self._translator_h2_in_constraints()
        self._translator_o2_in_constraints()
        self._translator_h2_out_constraints()
        self._translator_o2_out_constraints()
        self._add_ports()

    def _add_reaction(self):
        reaction_comps  = {"H2", "H2O", "O2"}

        self.reaction_props = GenericParameterBlock(
            **get_prop(reaction_comps, 
                       phases={"Vap"}, 
                       eos=self.config.reaction_eos,
                       ),
                       doc = 'Physical property package for the reaction'
        )

        self.rxn_props = GenericReactionParameterBlock(
            **get_rxn(self.reaction_props, 
                      {"h2_cmb"},
                      ),
                      doc = 'Reaction parameters'
        )

    def _add_unit_models(self):

        # 1) Translator for H2 stream
        self.translator_h2 = um.Translator(
            inlet_property_package       = self.config.anode_side_prop_package,
            inlet_property_package_args  = self.config.anode_side_prop_package_args,
            outlet_property_package      = self.reaction_props,
            outlet_state_defined         = False,
        )

        # 2) Air separator
        self.separator = um.Separator(
            property_package             = self.config.cathode_side_prop_package,
            outlet_list                  = ["o2_rich_strm", "o2_poor_strm"],
            split_basis                  = um.SplittingType.componentFlow,
        )

        # 3) Translator for O2-rich stream
        self.translator_o2 = um.Translator(
            inlet_property_package       = self.config.cathode_side_prop_package,
            outlet_property_package      = self.reaction_props,
            outlet_state_defined         = False,
        )

        # 4) Mixer for H2 and O2-rich stream
        self.mixer = um.Mixer(
            inlet_list                   = ["fuel_strm", "o2_rich_strm"],
            property_package             = self.reaction_props,
            momentum_mixing_type         = um.MomentumMixingType.none,
            energy_mixing_type           = um.MixingType.none,
        )

        # 5) SOEC reactor
        self.reactor = um.StoichiometricReactor(
            property_package             = self.reaction_props,
            reaction_package             = self.rxn_props,
            has_pressure_change          = False,
            has_heat_transfer            = True,
        )

        # 6) Separator after reactor (water/O2 separation)
        self.reactor_separator = um.Separator(
            property_package             = self.reaction_props,
            outlet_list                  = ["water_strm", "o2_strm"],
            split_basis                  = um.SplittingType.componentFlow,
        )

        # 7) Translator: water-rich stream → H2 property basis
        self.translator_h2_out = um.Translator(
            inlet_property_package       = self.reaction_props,
            outlet_property_package      = self.config.anode_side_prop_package,
            outlet_property_package_args = self.config.anode_side_prop_package_args,
            outlet_state_defined         = False,
        )

        # 8) Translator: O2-rich stream → air property basis
        self.translator_o2_out = um.Translator(
            inlet_property_package       = self.reaction_props,
            outlet_property_package      = self.config.cathode_side_prop_package,
            outlet_state_defined         = False,
        )

        # 9) Mixer for recycled O2 and N2 bleed
        self.mixer_out = um.Mixer(
            inlet_list                   = ["o2_strm_final", "o2_poor_strm"],
            momentum_mixing_type         = um.MomentumMixingType.none,
            property_package             = self.config.cathode_side_prop_package,
        )

        # 10) Heater for final gas conditioning
        self.heater = um.Heater(
            property_package             = self.config.cathode_side_prop_package,
            has_pressure_change          = False,
        )

    def _add_arcs(self):
        self.o2_to_translator                = Arc(source=self.separator.o2_rich_strm,
                                             destination=self.translator_o2.inlet)
        
        self.mix_h2                          = Arc(source=self.translator_h2.outlet,
                                                destination=self.mixer.fuel_strm)
        
        self.mix_o2                          = Arc(source=self.translator_o2.outlet,
                                                destination=self.mixer.o2_rich_strm)
        
        self.mix_to_reactor                  = Arc(source=self.mixer.outlet,
                                                destination=self.reactor.inlet)
        
        self.react_to_sep                    = Arc(source=self.reactor.outlet,
                                                destination=self.reactor_separator.inlet)
        
        self.separator_fuel_to_translator    = Arc(source=self.reactor_separator.water_strm,
                                                destination=self.translator_h2_out.inlet)
        
        self.separator_o2rich_to_translator  = Arc(source=self.reactor_separator.o2_strm,
                                                destination=self.translator_o2_out.inlet)
        
        self.o2_poor_to_heater               = Arc(source=self.separator.o2_poor_strm,
                                                destination=self.heater.inlet)
        
        self.heater_to_mixer_out             = Arc(source=self.heater.outlet,
                                                destination=self.mixer_out.o2_poor_strm)
        
        self.mix_o2_final                    = Arc(source=self.translator_o2_out.outlet,
                                                destination=self.mixer_out.o2_strm_final)
        
        # finally expand arcs
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

    def _add_variables(self):

        self.fuel_util = Var(
            self.flowsheet().time, 
            initialize=0.7, 
            bounds=(0, 1),
            units= pyo.units.dimensionless,
            doc="Fuel-utilisation factor",
        )
        self.i = Var(
            self.flowsheet().time, 
            initialize=0, 
            bounds=(0, 0.4),
            units=pyo.units.ampere / pyo.units.m**2, 
            doc="Current density",
        )
        self.cell_voltage = Var(
            self.flowsheet().time, 
            initialize=1.18, 
            bounds=(0.6, 1.3),
            units=pyo.units.volt, 
            doc="Cell voltage",
        )
        self.number_of_cells = Var(
            domain=pyo.NonNegativeIntegers,
            units=pyo.units.dimensionless, 
            doc="Number of cells in stack",
        )
        self.cell_area = Var(
            domain=pyo.NonNegativeReals,
            units= pyo.units.m**2, 
            doc="Single-cell active area",
        )
        self.sofc_power_dc = Var(
            self.flowsheet().time,
            initialize = 0,
            units= pyo.units.W, 
            doc="Gross DC power",
        )

    def _add_constraints(self):
        "Contraints for the model, the order is same as arcs excluded the translators"
        time = self.flowsheet().time
        t0 = self.flowsheet().time.first()

        self.separator.split_fraction[:, "o2_rich_strm", "O2"].fix(1)
        self.separator.split_fraction[:, "o2_rich_strm", "N2"].fix(0)

        self.reactor_separator.split_fraction[:, "o2_strm", "O2"].fix(1)
        self.reactor_separator.split_fraction[:, "o2_strm", "H2"].fix(0)
        self.reactor_separator.split_fraction[:, "o2_strm", "H2O"].fix(0)

        @self.mixer.Constraint(time)
        def pressure_eqn_mixer_in(b, t):
            return b.mixed_state[t].pressure == b.fuel_strm_state[t].pressure
        

        @self.reactor.Constraint(time)
        def fuel_utilization_eqn(b, t):
            return (
                b.control_volume.properties_out[t].flow_mol_comp["H2"]
                == b.control_volume.properties_in[t].flow_mol_comp["H2"]
                * (1 - self.fuel_util[t])
            )
        
        @self.Expression(time)
        def current_density_expr(b, t):
            return (
                (b.reactor.control_volume.properties_in[t].flow_mol_comp['H2']-
                b.reactor.control_volume.properties_out[t].flow_mol_comp['H2'] )
                * 2 * iconst.Constants.faraday_constant / (b.number_of_cells * b.cell_area)  
            )
        
        @self.Constraint(time)
        def current_density_eqn(b, t):
            return b.i[t] == b.current_density_expr[t]
        
        @self.Constraint(time)
        def voltage_eqn(b, t):
            return b.cell_voltage[t] == (
                   - 434.47 * b.i[t]**5
                   + 492.56 * b.i[t]**4
                   - 211.6  * b.i[t]**3
                   + 44.464 * b.i[t]**2
                   - 5.4157 * b.i[t]
                   + 1.184375
            )
        
        @self.Constraint(time)
        def sofc_power_dc_constraint(b, t):
            return b.sofc_power_dc[t] == b.i[t] * b.cell_voltage[t] * b.number_of_cells
        
        @self.Constraint(time)
        def sofc_energy_balance(b, t):
            #  LATO SINISTRO: calore “prodotto” dal reattore (segno cambiato)
            lhs = - pyo.units.convert(b.reactor.heat_duty[t], to_units=pyo.units.W)

            #  LATO DESTRO: potenza elettrica + calore scambiato con l’Heater
            rhs = b.sofc_power_dc[t] + pyo.units.convert(
                    b.heater.control_volume.heat[t], to_units=pyo.units.W)

            return lhs == rhs

        @self.Constraint(time)
        def heater_inequality(b, t):
            return b.heater.control_volume.heat[t] <= 0
    

        @self.Constraint(time)
        def dT_upper(b, t):
            return ( b.mixer_out.mixed_state[t].temperature
                - b.translator_h2_out.properties_out[t].temperature ) <= 10

        @self.Constraint(time)
        def dT_lower(b, t):
            return ( b.mixer_out.mixed_state[t].temperature
                - b.translator_h2_out.properties_out[t].temperature ) >= -10
        
        
        @self.mixer_out.Constraint(time)
        def pressure_eqn_mix_out(b, t):
            return b.mixed_state[t].pressure == b.o2_poor_strm_state[t].pressure

    @staticmethod
    def _general_translators_constraints(translator):
        
        @translator.Constraint(translator.flowsheet().time)
        def temperature_eqn(b, t):
            return b.properties_out[t].temperature == b.properties_in[t].temperature

        @translator.Constraint(translator.flowsheet().time)
        def flow_mol_eqn(b, t):
            return b.properties_out[t].flow_mol == b.properties_in[t].flow_mol

        @translator.Constraint(translator.flowsheet().time)
        def pressure_eqn(b, t):
            return b.properties_out[t].pressure == b.properties_in[t].pressure
         
    def _translator_h2_in_constraints(self):
        time = self.flowsheet().time
        t0 = self.flowsheet().time.first()

        comp_in = set(self.translator_h2.properties_in[t0].mole_frac_comp.keys())

        @self.translator_h2.Constraint(time, comp_in)
        def mole_frac_eqn_h2(b, t, j): 
            return b.properties_out[t].mole_frac_comp[j] == b.properties_in[t].mole_frac_comp[j]
        
        self._general_translators_constraints(self.translator_h2)

    def _translator_o2_in_constraints(self):
        time = self.flowsheet().time
        t0 = self.flowsheet().time.first()
            
        # OUTLET MOLAR FRACTIONS "BALANCE"
        comps_out = set(self.translator_o2.properties_out[t0].mole_frac_comp.keys())
        comps_out.remove('O2')
        @self.translator_o2.Constraint(time, comps_out)
        def mole_frac_eqn_o2_out_IN(b, t, j):
            return b.properties_out[t].mole_frac_comp[j] == 1e-19
        
        self._general_translators_constraints(self.translator_o2)

    def _translator_h2_out_constraints(self):
        time = self.flowsheet().time
        t0 = self.flowsheet().time.first()

        comps_out = set(self.translator_h2_out.properties_out[t0].mole_frac_comp.keys())
        comps_out.remove("H2")

        @self.translator_h2_out.Constraint(time, comps_out)
        def component_flow_eqn_h2_out(b, t, j):
            return (
                        b.properties_out[t].mole_frac_comp[j]
                        == b.properties_in[t].mole_frac_comp[j]
                    )
        
        self._general_translators_constraints(self.translator_h2_out)

    def _translator_o2_out_constraints(self):
        time = self.flowsheet().time
        t0 = self.flowsheet().time.first()

        # OUTLET MOLAR FRACTIONS "BALANCE"
        @self.translator_o2_out.Constraint(time)
        def mole_frac_eqn_o2_out_OUT(b, t):
            return b.properties_out[t].mole_frac_comp['N2'] == 1e-19
        
        self._general_translators_constraints(self.translator_o2_out)

    def _add_ports(self):
        """Add unit level ports"""
        self.add_inlet_port(
            name="anode_side_inlet",
            block=self.translator_h2.properties_in,
            doc="Hydrogen side inlet port",
        )
        self.add_inlet_port(
            name="cathode_side_inlet",
            block=self.separator.mixed_state,
            doc="Oxygen side inlet port",
        )
        self.add_outlet_port(
            name="anode_side_outlet",
            block=self.translator_h2.properties_out,
            doc="Hydrogen side outlet port",
        )
        self.add_outlet_port(
            name="cathode_side_outlet",
            block=self.mixer_out.mixed_state,
            doc="Oxygen side outlet port",
        )

    def initialize_build(
            self,
            outlvl=idaeslog.NOTSET,
            solver=None,
            optarg=None,
    ):

        import logging
        from pyomo.common.collections import ComponentSet
        from idaes.core.util.model_statistics import degrees_of_freedom
        from pyomo.util.infeasible import log_infeasible_constraints
        from idaes.core.util import from_json, to_json, StoreSpec
        from idaes.core.solvers import get_solver

        logging.getLogger("pyomo.core").setLevel(logging.INFO)
        logging.getLogger("pyomo.util.infeasible").setLevel(logging.INFO)

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        init_log.info_high("SOEC Initialization Starting")
        #solver_obj = get_solver(solver, optarg)
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        self.anode_side_inlet.fix()
        self.cathode_side_inlet.fix()
        self.anode_side_outlet.unfix()
        self.cathode_side_outlet.unfix()

        self.translator_h2.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.mix_h2)
        
        self.separator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.o2_to_translator)

        self.translator_o2.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.mix_o2)

        self.mixer.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.mix_to_reactor)

        self.reactor.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.react_to_sep)

        self.reactor_separator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        propagate_state(self.separator_fuel_to_translator)
        self.translator_h2_out.initialize(outlvl=outlvl, solver=solver, optarg=optarg) # FINE LATO ANODO  

        propagate_state(self.separator_o2rich_to_translator)
        self.translator_o2_out.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        propagate_state(self.o2_poor_to_heater)
        self.heater.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        propagate_state(self.heater_to_mixer_out)
        propagate_state(self.mix_o2_final)
        self.mixer_out.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        from_json(self, sd=istate, wts=sp)
        init_log.info_high("SOEC Initialization Complete")





    # def initialize_build(
    #         self,
    #         outlvl=idaeslog.NOTSET,
    #         solver=None,
    #         optarg=None,
    # ):

    #     import logging
    #     from idaes.core.util.initialization import propagate_state

    #     logging.getLogger("pyomo.core").setLevel(logging.INFO)
    #     logging.getLogger("pyomo.util.infeasible").setLevel(logging.INFO)

    #     init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
    #     init_log.info_high("SOFC Initialization Starting")

    #     units_sequence = [
    #         (self.translator_h2, self.mix_h2),
    #         (self.separator, self.o2_to_translator),
    #         (self.translator_o2, self.mix_o2),
    #         (self.mixer, self.mix_to_reactor),
    #         (self.reactor, self.react_to_sep),
    #         (self.reactor_separator, None),
    #         (self.translator_h2_out, self.separator_fuel_to_translator),
    #         (self.translator_o2_out, self.separator_o2rich_to_translator),
    #         (self.heater, self.o2_poor_to_heater),
    #         (self.mixer_out, [self.heater_to_mixer_out, self.mix_o2_final])
    #     ]

    #     for unit, arcs in units_sequence:
    #         init_log.info(f"Initializing {unit.name}, DOF before: {degrees_of_freedom(unit)}")
            
    #         # initialize unit
    #         unit.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
            
    #         init_log.info(f"DOF after initializing {unit.name}: {degrees_of_freedom(unit)}")
            
    #         # propagate states if necessary
    #         if arcs:
    #             if isinstance(arcs, list):
    #                 for arc in arcs:
    #                     propagate_state(arc)
    #             else:
    #                 propagate_state(arcs)

    #     # DOF check finale per tutto il modello
    #     total_dof = degrees_of_freedom(self)
    #     init_log.info_high(f"Initialization complete. Total DOF for the model: {total_dof}")



        
        
        








































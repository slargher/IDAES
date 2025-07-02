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
from ammonia_PR import (
    get_prop, get_rxn, EosType,
)
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom, report_statistics
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.logger as idaeslog
from pyomo.environ import Suffix

@declare_process_block_class("SofcUnit", doc="Simple SOFC model for process design.")
class SofcDesignData(UnitModelBlockData):

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    CONFIG.declare(
        "cathode_side_prop_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package for the air side.",
            doc=(
                "Property package for the air side, using "
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

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    CONFIG.declare(
        "anode_side_prop_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package for the fuel side.",
            doc=(
                "Property package for the fuel side, using "
                "idaes.models_extra.power_generation.properties.natural_gas_PR is "
                "strongly recommended, either Peng-Robinson or Ideal is okay"
            ),
        ),
    )

    CONFIG.declare(
        "anode_side_prop_package_args",
        ConfigBlock(
            implicit=True,
            description="Property package arguments for the fuel side.",
            doc="Property package arguments for the fuel side.",
        ),
    )
    

    CONFIG.declare(
    "recycle_side_prop_package",
    ConfigValue(
        domain=is_physical_parameter_block,
        description="Property package for the recycle side."
        ),
    )

    CONFIG.declare(
        "recycle_side_prop_package_args",
        ConfigBlock(
            implicit=True,
            description="Property package arguments for the recycle side.",
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
            doc="Indicates whether the SOFC is adiabatic. Default=False (adiabatic)",
        ),
    )

    def build(self):
        """Construct a unit model data block"""
        super().build()

        if set(self.config.anode_side_prop_package.component_list) != {
            "NH3", "H2", "N2",
        }:
            raise ConfigurationError(
                "SOFC anode side must contain exactly NH3, N2, H2"
            )
        
        if set(self.config.cathode_side_prop_package.component_list) != {
            'O2',
            'N2',
        }:
            raise ConfigurationError(
                'SOFC cathode side must contain exactly O2 and N2'
            )

        #self._add_recycle_package()
        self._add_combustion_reaction()
        self._add_cracking_reaction()
        self._add_unit_models()
        self._add_arcs()
        self._add_variables()
        self._add_constraints()
        self._translator_h2_constraints()
        self._translator_o2_constraints()
        self._translator_h2_out_constraints()
        self._translator_o2_out_constraints()
        self._translator_n2_constraints()
        #self._translator_nh3_constraints()
        self._add_ports()
        self._setup_scaling_factors()

    #def _add_recycle_package(self):
        #"""Add anode outlet property package"""
        #recycle_prop_package_comps  = {"H2", "H2O","N2","NH3"}
        #self.recycle_side_prop_package = GenericParameterBlock(
            #**get_prop(
                #recycle_prop_package_comps,
                #phases={"Vap"},
                #eos=self.config.recycle_prop_package_args.eos,
            #),
            #doc="Physical property package for the anode outlet",
        #)

    def _add_combustion_reaction(self):
        cmb_reaction_comps  = {"H2", "H2O", "O2"}

        self.cmb_reaction_props = GenericParameterBlock(
            **get_prop(cmb_reaction_comps, 
                       phases={"Vap"}, 
                       eos=self.config.reaction_eos,
                       ),
                       doc = 'Physical property package for the reaction'
        )

        self.cmb_rxn_props = GenericReactionParameterBlock(
            **get_rxn(self.cmb_reaction_props, 
                      {"h2_cmb"},
                      ),
                      doc = 'Reaction parameters'
        )

    def _add_cracking_reaction(self):
        # crk_reaction_comps  = {"NH3", "H2", "N2"}

        # self.crk_reaction_props = GenericParameterBlock(
        #     **get_prop(crk_reaction_comps, 
        #                phases={"Vap"}, 
        #                eos=self.config.reaction_eos,
        #                ),
        #                doc = 'Physical property package for the reaction'
        # )

        self.crk_rxn_props = GenericReactionParameterBlock(
            **get_rxn(self.config.anode_side_prop_package, 
                      {"nh3_crk"},
                      ),
                      doc = 'Reaction parameters'
        )

    def _add_unit_models(self):

        #1) Translator for NH3 stream add reaction package
        # self.translator_nh3 = um.Translator(
        #     inlet_property_package       = self.config.anode_side_prop_package, 
        #     inlet_property_package_args  = self.config.anode_side_prop_package_args,  
        #     outlet_property_package      = self.crk_reaction_props,
        #     outlet_state_defined         = False,
        # )

        # 2) NH3 cracker
        self.cracker = um.StoichiometricReactor(
            property_package             = self.config.anode_side_prop_package,
            reaction_package             = self.crk_rxn_props,
            has_pressure_change          = False,
            has_heat_transfer            = True,
            has_heat_of_reaction         = False,
        )

        # 3) NH3 separator into H2 and N2
        self.nh3_separator = um.Separator(
            property_package             = self.config.anode_side_prop_package, 
            outlet_list                  = ["h2_strm", "n2_strm"],
            split_basis                  = um.SplittingType.componentFlow,
        )

        # 4) Translator for N2 stream add prop package
        self.translator_n2 = um.Translator(
            inlet_property_package       = self.config.anode_side_prop_package,  
            outlet_property_package      = self.config.recycle_prop_package,
            outlet_state_defined         = False,
        )

        # 4) Translator for H2 stream add prop package
        self.translator_h2 = um.Translator(
            inlet_property_package       = self.config.anode_side_prop_package,  
            outlet_property_package      = self.cmb_reaction_props,
            outlet_state_defined         = False,
        )


        # 5) Air separator
        self.air_separator = um.Separator(
            property_package             = self.config.cathode_side_prop_package,
            outlet_list                  = ["o2_rich_strm", "o2_poor_strm"],
            split_basis                  = um.SplittingType.componentFlow,
        )

        # 6) Translator for O2-rich stream
        self.translator_o2 = um.Translator(
            inlet_property_package       = self.config.cathode_side_prop_package,
            outlet_property_package      = self.cmb_reaction_props,
            outlet_state_defined         = False,
        )

        # 7) Mixer for H2 and O2-rich stream
        self.mixer = um.Mixer(
            inlet_list                   = ["h2_strm", "o2_rich_strm"],
            property_package             = self.cmb_reaction_props,
            momentum_mixing_type         = um.MomentumMixingType.none,             
        )

        # 8) SOFC reactor
        self.reactor = um.StoichiometricReactor(
            property_package             = self.cmb_reaction_props,
            reaction_package             = self.cmb_rxn_props,
            has_pressure_change          = False,
            has_heat_transfer            = True,
            has_heat_of_reaction         = False,
        )

        # 9) Separator after reactor (water/O2 separation)
        self.reactor_separator = um.Separator(
            property_package             = self.cmb_reaction_props,
            outlet_list                  = ["water_strm", "o2_strm"],
            split_basis                  = um.SplittingType.componentFlow,
        )

        # 12) Heater for recycle side
        self.heater_n2 = um.Heater(
            property_package             = self.config.recycle_prop_package,
            has_pressure_change          = False,
        )

        # 10) Translator: water-rich stream → H2 property basis
        self.translator_h2_out = um.Translator(
            inlet_property_package       = self.cmb_reaction_props,
            outlet_property_package      = self.config.recycle_prop_package,
            outlet_state_defined         = False,
        )

        # 14) Mixer for recycled O2 and N2 mix
        self.mixer_out = um.Mixer(
            inlet_list                   = ["n2_strm", "h2o_strm"],
            property_package             = self.config.recycle_prop_package,
            momentum_mixing_type         = um.MomentumMixingType.none,
        )

        # 12) Heater for cathode side
        self.heater_air = um.Heater(
            property_package             = self.config.cathode_side_prop_package,
            has_pressure_change          = False,
        )

        # 13) Translator: O2-rich stream → air property basis
        self.translator_o2_out = um.Translator(
            inlet_property_package       = self.cmb_reaction_props,
            outlet_property_package      = self.config.cathode_side_prop_package,
            outlet_state_defined         = False,
        )


        # 11) Mixer for cathode air side
        self.mixer_air = um.Mixer(
            inlet_list                   = ["o2_strm", "o2_poor_strm"],
            property_package             = self.config.cathode_side_prop_package,
            momentum_mixing_type         = um.MomentumMixingType.none,             
        )


    def _add_arcs(self):
        # self.fuel_to_cracker                = Arc(source=self.translator_nh3.outlet,
        #                                      destination=self.cracker.inlet)
        
        self.cracker_to_sep               = Arc(source=self.cracker.outlet,
                                                destination=self.nh3_separator.inlet)
        
        self.sep_n2_to_trans               = Arc(source=self.nh3_separator.n2_strm,
                                                destination=self.translator_n2.inlet)
        
        self.sep_h2_to_trans               = Arc(source=self.nh3_separator.h2_strm,
                                                destination=self.translator_h2.inlet)
        
        self.trans_n2_to_heater_n2             = Arc(source=self.translator_n2.outlet,
                                                destination=self.heater_n2.inlet)
        
        self.heat_n2_to_mix_out             = Arc(source=self.heater_n2.outlet,
                                                  destination=self.mixer_out.h2o_strm)

        self.n2_air_to_heater_air              = Arc(source=self.air_separator.o2_poor_strm,
                                                destination=self.heater_air.inlet)

        self.heat_o2_to_mix_air               = Arc(source=self.heater_air.outlet,
                                                destination=self.mixer_air.o2_strm)
        
        self.sep_o2_to_trans                  = Arc(source=self.air_separator.o2_rich_strm,
                                                destination=self.translator_o2.inlet)
        
        self.o2_to_mixer                          = Arc(source=self.translator_o2.outlet,
                                                destination=self.mixer.o2_rich_strm)
        
        self.h2_to_mixer                          = Arc(source=self.translator_h2.outlet,
                                                destination=self.mixer.h2_strm)
        
        self.mix_to_reactor                  = Arc(source=self.mixer.outlet,
                                                destination=self.reactor.inlet)
        
        self.react_to_sep                    = Arc(source=self.reactor.outlet,
                                                destination=self.reactor_separator.inlet)
        
        self.sep_h2o_to_trans               = Arc(source=self.reactor_separator.water_strm,
                                                destination=self.translator_h2_out.inlet)
        
        self.react_sep_o2_to_trans                 = Arc(source=self.reactor_separator.o2_strm,
                                                destination=self.translator_o2_out.inlet)

        self.h2o_to_mixer_out                  = Arc(source=self.translator_h2_out.outlet,
                                                destination=self.mixer_out.h2o_strm)
        
        self.o2_to_mixer_air                  = Arc(source=self.translator_o2_out.outlet,
                                                destination=self.mixer_air.o2_strm)

        
        # finally expand arcs
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

    def _add_variables(self):

        self.fuel_util = Var(
            self.flowsheet().time, 
            initialize=0.77, 
            bounds=(0, 1),
            units= pyo.units.dimensionless,
            doc="Fuel-utilisation factor",
        )

        # Sara: does it make sense?
        self.ammonia_util = Var(
            self.flowsheet().time, 
            initialize=1.0, 
            bounds=(0, 1),
            units= pyo.units.dimensionless,
            doc="Ammonia-utilisation factor",
        )

        self.number_of_cells = Var(
            domain=pyo.NonNegativeIntegers,
            units=pyo.units.dimensionless, 
            doc="Number of cells in stack",
        )
        self.cell_area = Var(
            domain=pyo.NonNegativeReals,
            units= pyo.units.cm**2, 
            doc="Single-cell active area",
        )
        self.dT_Catout_Anout = Var(
            self.flowsheet().time,
            domain=pyo.NonNegativeReals, 
            units=pyo.units.K,
            doc="Max temperature differente btw Anode_outlet and Cathode_outlet",
        )
        self.i = Var(
            self.flowsheet().time, 
            initialize=0, 
            bounds=(0, 0.4),
            units=pyo.units.ampere / pyo.units.cm**2, 
            doc="Current density",
        )
        self.cell_voltage = Var(
            self.flowsheet().time, 
            initialize=1.18, 
            bounds=(0.6, 1.3),
            units=pyo.units.volt, 
            doc="Cell voltage",
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

        self.air_separator.split_fraction[:, "o2_rich_strm", "O2"].fix(1)
        self.air_separator.split_fraction[:, "o2_rich_strm", "N2"].fix(0)

        self.reactor_separator.split_fraction[:, "o2_strm", "O2"].fix(1)
        self.reactor_separator.split_fraction[:, "o2_strm", "H2"].fix(0)
        self.reactor_separator.split_fraction[:, "o2_strm", "H2O"].fix(0)

        self.nh3_separator.split_fraction[:, "h2_strm", "H2"].fix(1)
        self.nh3_separator.split_fraction[:, "h2_strm", "N2"].fix(0)
        self.nh3_separator.split_fraction[:, "h2_strm", "NH3"].fix(0)


        #SARA:does it make sense?
        @self.cracker.Constraint(time)
        def ammonia_utilization_eqn(b, t):
             return (
                  b.control_volume.properties_out[t].flow_mol_comp["NH3"]
                  == b.control_volume.properties_in[t].flow_mol_comp["NH3"]
                  * (1.0 - self.ammonia_util[t])
              )

    
        #sara:added new
        @self.mixer_air.Constraint(time)
        def pressure_eqn_mixer_air(b, t):
            return b.mixed_state[t].pressure == b.o2_poor_strm_state[t].pressure
        
        @self.mixer.Constraint(time)
        def pressure_eqn_mixer(b, t):
            return b.mixed_state[t].pressure == b.o2_rich_strm_state[t].pressure
        
        @self.mixer_out.Constraint(time)
        def pressure_eqn_mixer_out(b, t):
            return b.mixed_state[t].pressure == b.n2_strm_state[t].pressure


        
        @self.reactor.Constraint(time)
        def fuel_utilization_eqn(b, t):
            return (
                b.control_volume.properties_out[t].flow_mol_comp["H2"]
                == b.control_volume.properties_in[t].flow_mol_comp["H2"]
                * (1.0 - self.fuel_util[t])
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
        
        #From fit of polarization curve just for now use polynmial
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
            return b.sofc_power_dc[t] == b.i[t] * b.cell_voltage[t] * b.number_of_cells * b.cell_area

        @self.Constraint(time)
        def sofc_energy_balance(b, t):
            return (pyo.units.convert(b.heater_n2.control_volume.heat[t], to_units=pyo.units.W) +  pyo.units.convert(b.heater_air.control_volume.heat[t], to_units=pyo.units.W) 
                   == - b.sofc_power_dc[t] - pyo.units.convert(b.reactor.control_volume.heat[t], to_units=pyo.units.W) + pyo.units.convert(b.cracker.control_volume.heat[t], to_units=pyo.units.W))
  
        @self.Constraint(time)
        def dT_upper_an(b, t):
            return ( b.cracker.control_volume.properties_in[t].temperature
                - b.translator_h2_out.properties_out[t].temperature ) <= b.dT_Catout_Anout[t]

        @self.Constraint(time)
        def dT_lower_an(b, t):
            return ( b.cracker.control_volume.properties_in[t].temperature
                - b.translator_h2_out.properties_out[t].temperature ) >= -b.dT_Catout_Anout[t]
        
        @self.Constraint(time)
        def dT_upper_cat(b, t):
            return ( b.air_separator.mixed_state[t].temperature
                - b.mixer_air.mixed_state[t].temperature ) <= b.dT_Catout_Anout[t]

        @self.Constraint(time)
        def dT_lower_cat(b, t):
            return ( b.air_separator.mixed_state[t].temperature
                - b.mixer_air.mixed_state[t].temperature ) >= -b.dT_Catout_Anout[t]
  



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
         
    # def _translator_nh3_constraints(self):
    #      time = self.flowsheet().time
    #      t0 = self.flowsheet().time.first()
    
    #      comp_in = set(self.translator_nh3.properties_in[t0].mole_frac_comp.keys())
    #      @self.translator_nh3.Constraint(time, comp_in)
    #      def mole_frac_eqn_nh3(b, t, j): 
    #          return b.properties_out[t].mole_frac_comp[j] == b.properties_in[t].mole_frac_comp[j]
         
    #      self._general_translators_constraints(self.translator_nh3)

    def _translator_o2_constraints(self):
        time = self.flowsheet().time
        t0 = self.flowsheet().time.first()
            
        # OUTLET MOLAR FRACTIONS "BALANCE"
        comps_out = set(self.translator_o2.properties_out[t0].mole_frac_comp.keys())
        comps_out.remove('O2')
        @self.translator_o2.Constraint(time, comps_out)
        def mole_frac_eqn_o2_out_IN(b, t, j):
            return b.properties_out[t].mole_frac_comp[j] == 1e-19
        
        self._general_translators_constraints(self.translator_o2)

    #to check
    def _translator_n2_constraints(self):
        time = self.flowsheet().time
        t0 = self.flowsheet().time.first()

        comp_in = set(self.translator_n2.properties_in[t0].mole_frac_comp.keys())
        @self.translator_n2.Constraint(time, comp_in)
        def mole_frac_eqn_n2(b, t, j): 
            return b.properties_out[t].mole_frac_comp[j] == b.properties_in[t].mole_frac_comp[j]
        
        self._general_translators_constraints(self.translator_n2)


    def _translator_h2_constraints(self):
        time = self.flowsheet().time
        t0 = self.flowsheet().time.first()

        comps_in = set(self.translator_h2.properties_in[t0].mole_frac_comp.keys())
        comps_out = set(self.translator_h2.properties_out[t0].mole_frac_comp.keys())

        # Remove unwanted components from output side
        comps_out -= {"N2", "NH3"}

        # Take intersection to get only components present in both sides
        comps = comps_in.intersection(comps_out)

        @self.translator_h2.Constraint(time, comps)
        def component_flow_eqn_h2_reactor(b, t, j):
            return b.properties_out[t].mole_frac_comp[j] == b.properties_in[t].mole_frac_comp[j]

        self._general_translators_constraints(self.translator_h2)

    def _translator_h2_out_constraints(self):
        time = self.flowsheet().time
        t0 = self.flowsheet().time.first()
    
        @self.translator_h2_out.Constraint(time)
        def mole_frac_eqn_h2_0(b, t): 
            return b.properties_out[t].mole_frac_comp['H2O'] == b.properties_in[t].mole_frac_comp['H2O']
        
        comp_out = set(self.translator_h2_out.properties_out[t0].mole_frac_comp.keys())
        comp_out.remove("H2O")
        comp_out.remove("H2")
        
        @self.translator_h2_out.Constraint(time, comp_out)
        def mole_frac_eqn_h2(b, t, j): 
            return b.properties_out[t].mole_frac_comp[j] == 1e-19
        
        self._general_translators_constraints(self.translator_h2_out)



    def _translator_o2_out_constraints(self):
        time = self.flowsheet().time
        t0 = self.flowsheet().time.first()
            
        # OUTLET MOLAR FRACTIONS "BALANCE"
        comps_out = set(self.translator_o2_out.properties_out[t0].mole_frac_comp.keys())
        comps_out.remove('O2')
        @self.translator_o2_out.Constraint(time, comps_out)
        def mole_frac_eqn_o2_out(b, t, j):
            return b.properties_out[t].mole_frac_comp[j] == 1e-19
        
        self._general_translators_constraints(self.translator_o2_out)


    def _add_ports(self):
        """Add unit level ports"""
        # self.add_inlet_port(
        #     name="anode_side_inlet",
        #     block=self.translator_nh3.properties_in,
        #     doc="Ammonia side inlet port",
        # )

        self.add_inlet_port(
            name="anode_side_inlet",
            block=self.cracker.control_volume.properties_in,
            doc="Ammonia side inlet port",
        )

        self.add_inlet_port(
            name="cathode_side_inlet",
            block=self.air_separator.mixed_state,
            doc="Oxygen side inlet port",
        )
        self.add_outlet_port(
            name="anode_side_outlet",
            block=self.mixer_out.mixed_state,
            doc="Hydrogen side outlet port",
        )
        self.add_outlet_port(
            name="cathode_side_outlet",
            block=self.mixer_air.mixed_state,
            doc="Oxygen side outlet port",
        )

    def _setup_scaling_factors(self):
        """
        Simple one-stop scaling routine for the SOFC example flowsheet.

        Call *after* you create `model = sofc_example_flowsheet()` and
        *before* you solve.
        """

        # ------------------------------------------------------------------
        # 1)  Default scaling for every property package we use
        # ------------------------------------------------------------------
        for pp in [
            self.config.anode_side_prop_package,     # H2/H2O
            self.config.cathode_side_prop_package, # O2/N2 (air)
            self.config.recycle_prop_package,  
            self.cmb_reaction_props,
            self.crk_rxn_props,         # common reaction basis
        ]:
            pp.set_default_scaling("flow_mol", 1e1)
            pp.set_default_scaling("flow_mol_phase", 1e1)
            pp.set_default_scaling("flow_mol_phase_comp", 1e1)
            pp.set_default_scaling("temperature", 1e-2)
            pp.set_default_scaling("pressure", 1e-2)
            pp.set_default_scaling("mole_frac_comp", 1e2)
            pp.set_default_scaling("mole_frac_phase_comp", 1e2)
            pp.set_default_scaling("enth_mol_phase", 1e-5)   # ≈ 2 × 10⁵ J mol⁻¹
            pp.set_default_scaling("entr_mol_phase", 1e-5)

        # ------------------------------------------------------------------
        # 2)  Specific variables that the defaults don’t catch
        # ------------------------------------------------------------------
        # -- reaction extent (mol s⁻¹) in the SOFC reactor
        # Scaling reaction extents in the SOFC reactor (combustion)
        for (t, rxn), var in self.reactor.control_volume.rate_reaction_extent.items():
            if rxn == "h2_cmb":
                iscale.set_scaling_factor(var, 50)  # scaling for hydrogen combustion

        # Scaling reaction extents in the ammonia cracker
        for (t, rxn), var in self.cracker.control_volume.rate_reaction_extent.items():
            if rxn == "nh3_crk":
                iscale.set_scaling_factor(var, 100)  # scaling for ammonia cracking

        # ------------------------------------------------------------------
        # 3)  Compute & apply the scaling – and silence the warnings
        # ------------------------------------------------------------------
        idaeslog.getLogger("idaes.core.util.scaling",
                        level=idaeslog.ERROR)              

        iscale.calculate_scaling_factors(self)                   

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

        logging.getLogger("pyomo.core").setLevel(logging.INFO)
        logging.getLogger("pyomo.util.infeasible").setLevel(logging.INFO)

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")

        init_log.info_high("SOEC Initialization Starting")

        self.anode_side_inlet.fix()
        self.cathode_side_inlet.fix()
        self.anode_side_outlet.unfix()
        self.cathode_side_outlet.unfix()


        # Initialize translator for NH3 stream
        #self.translator_nh3.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate state NH3 → Cracker inlet
        #propagate_state(self.fuel_to_cracker)

        # Initialize NH3 cracker
        self.cracker.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate state Cracker → NH3 separator
        propagate_state(self.cracker_to_sep)

        # Initialize NH3 separator
        self.nh3_separator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate separator N2 stream → N2 translator inlet
        propagate_state(self.sep_n2_to_trans)

        # Initialize N2 translator
        self.translator_n2.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate separator H2 stream → H2 translator inlet
        propagate_state(self.sep_h2_to_trans)

        # Initialize H2 translator
        self.translator_h2.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate N2 translator → Heater N2 inlet
        propagate_state(self.trans_n2_to_heater_n2)

        # Initialize Heater for N2
        self.heater_n2.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate heater N2 outlet → Mixer out (N2 stream)
        propagate_state(self.heat_n2_to_mix_out)

        # Initialize air separator
        self.air_separator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate O2-poor stream → Heater for air
        propagate_state(self.n2_air_to_heater_air)

        # Initialize heater for cathode side air
        self.heater_air.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate heated O2-poor → mixer_air
        propagate_state(self.heat_o2_to_mix_air)

        # Propagate air separator O2-rich → Translator O2 inlet
        propagate_state(self.sep_o2_to_trans)

        # Initialize translator for O2-rich stream
        self.translator_o2.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate translator O2 outlet → Mixer O2-rich inlet
        propagate_state(self.o2_to_mixer)

        # Propagate translator H2 outlet → Mixer H2 inlet
        propagate_state(self.h2_to_mixer)

        # Initialize mixer for H2 and O2-rich stream
        self.mixer.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate mixer outlet → Reactor inlet
        propagate_state(self.mix_to_reactor)

        # Initialize SOFC reactor
        self.reactor.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate reactor outlet → Reactor separator inlet
        propagate_state(self.react_to_sep)

        # Initialize reactor separator
        self.reactor_separator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate reactor separator water stream → Translator H2 out inlet
        propagate_state(self.sep_h2o_to_trans)

        # Initialize translator for H2O → recycle basis
        self.translator_h2_out.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate translator H2 out → mixer_out
        propagate_state(self.h2o_to_mixer_out)

        # Propagate reactor separator O2 stream → Translator O2 out inlet
        propagate_state(self.react_sep_o2_to_trans)

        # Initialize translator for O2 out
        self.translator_o2_out.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Propagate translator O2 out → mixer_air
        propagate_state(self.o2_to_mixer_air)

        # Initialize mixer for cathode air side
        self.mixer_air.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        # Initialize mixer for recycled N2 and H2O stream
        self.mixer_out.initialize(outlvl=outlvl, solver=solver, optarg=optarg)





    
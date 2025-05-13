
import pyomo.environ as pyo
from pyomo.network import Arc
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

import idaes.core.util.constants as iconst
from idaes.core import UnitModelBlockData, declare_process_block_class
from idaes.core.util.config import is_physical_parameter_block
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)
import idaes.models.unit_models as um  # um = unit models
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver


#Peng Robinson EOS
# The Peng Robinson EOS is a cubic equation of state that is commonly used to model the
# thermodynamic behavior of natural gas mixtures.           
from ammonia_PR import (
    get_prop,
    get_rxn,
    EosType,
)
"""check these file on github IDAES to see how to adapt them for sofc"""
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.logger as idaeslog
from idaes.core.util import from_json, to_json, StoreSpec

@declare_process_block_class("Nh3SofcDesign", doc="Simple NH3-SOFC model for process design.")
class Nh3SofcDesignData(UnitModelBlockData):
    """Simple 0D NH3-SOFC model. This is used for design point flowsheets or surrogate model basis."""

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    CONFIG.declare(
        "air_side_property_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package for the air side (O₂, N₂, H₂O).",
        ),
    )
    CONFIG.declare(
        "fuel_side_property_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package for the fuel side (NH₃, H₂, N₂, H₂O).",
        ),
    )
    CONFIG.declare(
        "air_side_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Property package arguments for the air side.",
        ),
    )
    CONFIG.declare(
        "fuel_side_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Property package arguments for the fuel side.",
        ),
    )
    CONFIG.declare(
        "reaction_eos",
        ConfigValue(
            default=EosType.PR,  # Use PR (Peng-Robinson) EOS for real gas behavior
            domain=In(EosType),
            description="Physical properties for fuel cell reactions",
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

        # Validation for the fuel side (NH₃) 
        if set(self.config.fuel_side_property_package.component_list) != {
            "NH3",
            "H2",
            "N2",
            "H2O",
        }:
            raise ConfigurationError(
                "NH₃-SOFC fuel side must contain exactly NH3, H2, N2, and H2O"
            )
        # Validation for the air side (must have oxygen in it)
        if not "O2" in self.config.air_side_property_package.component_list:
            raise ConfigurationError("NH₃-SOFC air side must contain O2")

        # build the block, use smaller methods to make it easier to follow

        #self._add_electrolysis_properties()
        self._add_fuel_cell_properties() 
        self._add_unit_models()
        self._add_arcs()
        self._add_variables()
        self._add_constraints()
        self._add_air_translator_constraints()
        self._add_fuel_inlet_translator_constraints()  
        self._add_fuel_outlet_translator_constraints()
        self._add_ports()
        self._scaling_guess()


    """IMPORTANT: create a document like natural_gas_PR.py where I have the components, the reactions and 
    I can use the get_props and get_rxn """

    def _add_fuel_cell_properties(self):
        """Add the fuel cell property package. This includes NH₃, H₂, N₂, H₂O, and O₂
        and is used internally by the unit for the SOFC reaction.
        """
        # Property parameters for the NH₃-SOFC fuel side (NH₃, H₂, N₂, H₂O) and air side (O₂)
        self.fuel_cell_prop_params = GenericParameterBlock(
            **get_prop(
                {"NH3", "H2", "N2", "H2O", "O2"},  # Include NH₃, H₂, N₂, H₂O, and O₂ components
                phases={"Vap"},  # Phase is typically vapor for gases in fuel cells
                eos=self.config.reaction_eos,  # Use the chosen equation of state (e.g., PR or IDEAL)
            ),
            doc="Physical property parameters for the NH₃-SOFC fuel cell reaction",
        )

        # Reaction parameters for NH₃ oxidation (anode) and oxygen reduction (cathode)
        self.fuel_cell_rxn_params = GenericReactionParameterBlock(
            **get_rxn(self.fuel_cell_prop_params, {"nh3_crk", "h2_cmb"}),  # Define NH₃ oxidation and O₂ reduction reactions
            doc="Reaction parameters for NH₃-SOFC",
        )
        """can I call two function?"""

    def _add_unit_models(self):
        """Add the unit models that make up the NH₃-SOFC system.
        """
         # Translate NH₃ fuel side stream to fuel cell property package
        self.fuel_inlet_translator = um.Translator(
            doc="Translate nh3 fuel side properties to fuel cell property package",
            inlet_property_package=self.config.fuel_side_property_package,
            inlet_property_package_args=self.config.fuel_side_property_package_args,
            outlet_property_package=self.fuel_cell_prop_params,
            outlet_state_defined=False,
        )

        # Translate air/oxygen side stream to fuel cell property package
        self.air_inlet_translator = um.Translator(
            doc="Translate air/oxygen side properties to SOFC property package",
            inlet_property_package=self.config.air_side_property_package,
            inlet_property_package_args=self.config.air_side_property_package_args,
            outlet_property_package=self.fuel_cell_prop_params,
            outlet_state_defined=False,
        )

        # SOFC reactor: handles both NH₃ oxidation and O₂ reduction
        self.fuel_cell_reactor = um.StoichiometricReactor(
            doc="NH₃-SOFC reactor (anode and cathode)",
            property_package=self.fuel_cell_prop_params,
            reaction_package=self.fuel_cell_rxn_params,
            has_pressure_change=False,
            has_heat_transfer=True,
        )

        self.air_separator = um.Separator(
            property_package=self.fuel_cell_prop_params,
            split_basis=um.SplittingType.componentFlow,
            outlet_list=["fuel_strm", "air_strm"],
        )
        self.fuel_outlet_translator = um.Translator(
            doc="Translate sofc properties to nh3 properties",
            inlet_property_package=self.fuel_cell_prop_params,
            outlet_property_package=self.config.fuel_side_property_package,
            outlet_property_package_args=self.config.fuel_side_property_package_args,
            outlet_state_defined=False,
        )
        self.air_translator = um.Translator(
            doc="Translate electrolysis properties to air properties",
            inlet_property_package=self.fuel_cell_prop_params,
            outlet_property_package=self.config.air_side_property_package,
            outlet_property_package_args=self.config.air_side_property_package_args,
            outlet_state_defined=False,
        )
        self.air_mixer = um.Mixer(
            property_package=self.config.air_side_property_package,
            property_package_args=self.config.air_side_property_package_args,
            momentum_mixing_type=um.MomentumMixingType.none,
            inlet_list=["sweep_strm", "air_strm"],
        )
        self.sweep_heater = um.Heater(
            property_package=self.config.air_side_property_package,
            property_package_args=self.config.air_side_property_package_args,
        )
        

    def _add_arcs(self):
        """Add streams to connect internal units in the NH₃-SOFC model."""
    
        self.strm1 = Arc(
            doc="Fuel NH₃ side inlet to SOFC reactor",
            source=self.fuel_inlet_translator.outlet,
            destination=self.sofc_reactor.inlet,
        )

        self.strm2 = Arc(
            doc="sofc reactor to air separation",
            source=self.sofc_reactor.outlet,
            destination=self.air_separator.inlet,
        )
        self.strm3 = Arc(
            doc="air separation to fuel side outlet translator",
            source=self.air_separator.fuel_strm,
            destination=self.fuel_outlet_translator.inlet,
        )
        self.strm4 = Arc(
            doc="air from sofc to translator to sweep properties",
            source=self.air_separator.o2_strm,
            destination=self.o2_translator.inlet,
        )
        self.strm5 = Arc(
            doc="air translator to air/sweep mixer",
            source=self.air_translator.outlet,
            destination=self.oair_mixer.air_strm,
        )
        self.strm6 = Arc(
            doc="air/sweep mixer to sweep heater",
            source=self.air_mixer.outlet,
            destination=self.sweep_heater.inlet,
        )

        # Probably don't need this and may need to call it again for dynamic
        # models, but this will make it easy for a user to have a flowsheet
        # with one steady state model in it.
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)
    
        
    def _add_variables(self):
        """Add variables for the NH₃-SOFC unit model level constraints."""
        
        self.current = pyo.Var(
            self.flowsheet().time,
            initialize=0,
            units=pyo.units.ampere,
            doc="SOFC electric current",
        )

        self.cell_potential = pyo.Var(
            self.flowsheet().time,
            initialize=0.8,
            units=pyo.units.volts,
            doc="Electric potential of a single SOFC cell",
        )

        self.fuel_utilization = pyo.Var(
            self.flowsheet().time,
            units=pyo.units.dimensionless,
            doc="Fraction of NH₃ fuel utilized in the SOFC",
        )

        self.heat = pyo.Var(
            self.flowsheet().time,
            units=pyo.units.W,
            initialize=0,
            doc="Net heat transferred in the SOFC",
        )

        self.fuel_side_outlet_temperature = pyo.Reference(
            self.sofc_reactor.control_volume.properties_out[:].temperature
        )

        self.air_side_outlet_temperature = pyo.Reference(
            self.sweep_heater.control_volume.properties_out[:].temperature
        )

    def _add_constraints(self):
        """Add unit model constraints for nh3-sofc"""

        @self.Constraint(self.flowsheet().time)
        def fuel_utilization_eqn(b, t):
            return b.sofc_reactor.control_volume.properties_out[t].flow_mol_comp["NH3"] == (
                b.sofc_reactor.control_volume.properties_in[t].flow_mol_comp["NH3"]
                * (1.0 - b.fuel_utilization[t])
            )

        @self.Expression(self.flowsheet().time)
        def current_expr(b, t):
            nh3_in = b.sofc_reactor.control_volume.properties_in[t].flow_mol_comp["NH3"]
            nh3_out = b.sofc_reactor.control_volume.properties_out[t].flow_mol_comp["NH3"]
            nh3_reacted = nh3_in - nh3_out
            return 3 * nh3_reacted * iconst.Constants.faraday_constant

        @self.Constraint(self.flowsheet().time)
        def current_eqn(b, t):
            return b.current[t] == b.current_expr[t]

        # Fix O2 separator to remove O2 from fuel side (after NH₃ oxidation)
        self.air_separator.split_fraction[:, "air_strm", "O2"].fix(1)
        self.air_separator.split_fraction[:, "air_strm", "H2O"].fix(0)
        self.air_separator.split_fraction[:, "air_strm", "H2"].fix(0)
        self.air_separator.split_fraction[:, "air_strm", "NH3"].fix(0)
        self.air_separator.split_fraction[:, "air_strm", "N2"].fix(0)

        @self.air_mixer.Constraint(self.flowsheet().time)
        def pressure_eqn(b, t):
            return b.mixed_state[t].pressure == b.sweep_strm_state[t].pressure

        if not self.config.has_heat_transfer:

            @self.Constraint(self.flowsheet().time)
            def heat_transfer_eqn(b, t):
                return b.heat[t] == 0

        @self.Expression(self.flowsheet().time)
        def delta_enth(b, t):
            return (
                b.air_separator.fuel_strm_state[t].get_enthalpy_flow_terms("Vap")
                + b.sweep_heater.control_volume.properties_out[t].get_enthalpy_flow_terms("Vap")
                - b.sofc_reactor.control_volume.properties_in[t].get_enthalpy_flow_terms("Vap")
                - b.air_mixer.sweep_strm_state[t].get_enthalpy_flow_terms("Vap")
                + b.heat[t]
            )

        @self.Constraint(self.flowsheet().time)
        def cell_potential_eqn(b, t):
            return b.cell_potential[t] == b.delta_enth[t] / b.current[t]

    @staticmethod
    def _translator_ftp_constraints(translator):
        """Generalize some translator block constraints"""

        @translator.Constraint(translator.flowsheet().time)
        def temperature_eqn(b, t):
            return b.properties_out[t].temperature == b.properties_in[t].temperature

        @translator.Constraint(translator.flowsheet().time)
        def flow_mol_eqn(b, t):
            return b.properties_out[t].flow_mol == b.properties_in[t].flow_mol

        @translator.Constraint(translator.flowsheet().time)
        def pressure_eqn(b, t):
            return b.properties_out[t].pressure == b.properties_in[t].pressure

    def _add_air_translator_constraints(self):
        """Add constraints to translate from fuel cell O₂ stream to sweep properties."""
        t0 = self.flowsheet().time.first()
        comps = set(self.air_translator.properties_out[t0].mole_frac_comp.keys())
        comps.discard("O2")

        @self.air_translator.Constraint(self.flowsheet().time, comps)
        def mole_frac_comp_eqn(b, t, c):
            return b.properties_out[t].mole_frac_comp[c] == 1e-19

        self.air_translator.properties_out[:].mole_frac_comp["O2"] = 1

        self._translator_ftp_constraints(self.air_translator)

    def _add_fuel_inlet_translator_constraints(self):
        """Translate the inlet hydrogen stream to reactor side properties. This
        stream contains NH3, H2, N2, and possibly H2O.
        """
        t0 = self.flowsheet().time.first()
        comps = set(self.fuel_inlet_translator.properties_in[t0].mole_frac_comp.keys())

        @self.fuel_inlet_translator.Constraint(self.flowsheet().time, comps)
        def mole_frac_comp_eqn(b, t, c):
            return (
                b.properties_out[t].mole_frac_comp[c]
                == b.properties_in[t].mole_frac_comp[c]
            )

        self._translator_ftp_constraints(self.nh3_inlet_translator)

    def _add_anode_outlet_translator_constraints(self):
        """Translate SOFC anode outlet properties to ammonia-side outlet properties.
        Includes H2, H2O, N2.
        """
        t0 = self.flowsheet().time.first()
        comps = set(self.anode_outlet_translator.properties_out[t0].mole_frac_comp.keys())
        comps.remove("H2")  # sum = 1 constraint

        @self.anode_outlet_translator.Constraint(self.flowsheet().time, comps)
        def mole_frac_comp_eqn(b, t, c):
            return (
                b.properties_out[t].mole_frac_comp[c]
                == b.properties_in[t].mole_frac_comp[c]
            )

        self._translator_ftp_constraints(self.anode_outlet_translator)


    def _add_ports(self):
        """Add unit level ports"""
        self.add_inlet_port(
            name="fuel_side_inlet",
            block=self.fuel_inlet_translator.properties_in,
            doc="nh3 side inlet port",
        )
        self.add_inlet_port(
            name="air_side_inlet",
            block=self.air_mixer.sweep_strm_state,
            doc="Oxygen side inlet port (air/sweep gas)",
        )
        self.add_outlet_port(
            name="fuel_side_outlet",
            block=self.anode_outlet_translator.properties_out,
            doc="fuel side outlet port",
        )
        self.add_outlet_port(
            name="oxygen_side_outlet",
            block=self.sweep_heater.control_volume.properties_out,
            doc="Oxygen side outlet port",
        )


    def set_flow_scale(self, scale=1):
        """Set default flow scaling in the property packages based on the
        expected magnitude of typical flows

        Args:
            scale (float): Expected molar flow variable scale

        Returns:
            None
        """

        self.config.air_side_property_package.set_default_scaling("flow_mol", scale)
        self.config.fuel_side_property_package.set_default_scaling("flow_mol", scale)
        self.sofc_prop_params.set_default_scaling("flow_mol", scale)

        self.config.air_side_property_package.set_default_scaling("flow_mol_phase", scale)
        self.config.fuel_side_property_package.set_default_scaling("flow_mol_phase", scale)
        self.sofc_prop_params.set_default_scaling("flow_mol_phase", scale)

        for v in self.sofc_reactor.control_volume.rate_reaction_extent.values():
            iscale.set_scaling_factor(v, scale)
        for v in self.sofc_reactor.control_volume.rate_reaction_generation.values():
            iscale.set_scaling_factor(v, scale)


    def set_heat_scale(self, scale=1e-5):
        """Set the heat transfer scale factor roughly based on the size of the
        process.

        Args:
            scale (float): Expected heat transfer scale

        Returns:
            None
        """
        iscale.set_scaling_factor(self.sofc_reactor.control_volume.heat, scale)
        iscale.set_scaling_factor(self.sweep_heater.control_volume.heat, scale)
        iscale.set_scaling_factor(self.heat, scale)


    def set_current_scale(self, scale=1e-6):
        """Set the expected electrical current scale

        Args:
            scale (float): Expected electrical current scale

        Returns:
            None
        """
        iscale.set_scaling_factor(self.current, scale)


    def _scaling_guess(self):
        self.set_flow_scale()
        self.set_heat_scale()
        self.set_current_scale()

        # For components that exist we don't expect particularly low
        # concentrations, so change the default to 10, which should be good
        # for mole fractions > 0.01 which is what we expect here.
        self.config.air_side_property_package.set_default_scaling(
            "mole_frac_comp", 10
        )
        self.config.fuel_side_property_package.set_default_scaling(
            "mole_frac_comp", 10
        )
        self.sofc_prop_params.set_default_scaling("mole_frac_comp", 10)

        self.config.air_side_property_package.set_default_scaling(
            "mole_frac_phase_comp", 10
        )
        self.config.fuel_side_property_package.set_default_scaling(
            "mole_frac_phase_comp", 10
        )
        self.sofc_prop_params.set_default_scaling("mole_frac_phase_comp", 10)

        # Set some other scale factors that we have a good guess for
        # Set scaling for enthalpy terms across major components
        unt = self.sofc_reactor
        iscale.set_scaling_factor(
            unt.control_volume.properties_in[0.0].enth_mol_phase["Vap"], 1e-4
        )
        iscale.set_scaling_factor(
            unt.control_volume.properties_out[0.0].enth_mol_phase["Vap"], 1e-4
        )

        unt = self.air_separator
        iscale.set_scaling_factor(unt.h2_strm_state[0.0].enth_mol_phase["Vap"], 1e-4)

        unt = self.air_mixer
        iscale.set_scaling_factor(unt.o2_strm_state[0.0].enth_mol_phase["Vap"], 1e-4)
        iscale.set_scaling_factor(unt.mixed_state[0.0].enth_mol_phase["Vap"], 1e-4)
        iscale.set_scaling_factor(unt.sweep_strm_state[0.0].enth_mol_phase["Vap"], 1e-4)

        unt = self.sweep_heater
        iscale.set_scaling_factor(
            unt.control_volume.properties_in[0.0].enth_mol_phase["Vap"], 1e-4
        )
        iscale.set_scaling_factor(
            unt.control_volume.properties_out[0.0].enth_mol_phase["Vap"], 1e-4
        )

    def calculate_scaling_factors(self):
        """Calculate scale factors for the unit model equations"""
        for t in self.flowsheet().time:
            iscale.constraint_scaling_transform(
                self.fuel_utilization_eqn[t],  # assumes you've renamed this constraint
                iscale.get_scaling_factor(
                    self.sofc_reactor.control_volume.properties_in[t].flow_mol
                ),
            )
            iscale.constraint_scaling_transform(
                self.current_eqn[t], iscale.get_scaling_factor(self.current[t])
            )
            iscale.constraint_scaling_transform(
                self.air_mixer.pressure_eqn[t],
                iscale.get_scaling_factor(self.air_mixer.mixed_state[t].pressure),
            )
            iscale.constraint_scaling_transform(
                self.air_translator.pressure_eqn[t],
                iscale.get_scaling_factor(self.air_translator.properties_in[t].pressure),
            )
            iscale.constraint_scaling_transform(
                self.air_translator.temperature_eqn[t],
                iscale.get_scaling_factor(self.air_translator.properties_in[t].temperature),
            )
            iscale.constraint_scaling_transform(
                self.air_translator.flow_mol_eqn[t],
                iscale.get_scaling_factor(self.air_translator.properties_in[t].flow_mol),
            )
            iscale.constraint_scaling_transform(
                self.fuel_inlet_translator.pressure_eqn[t],
                iscale.get_scaling_factor(self.fuel_inlet_translator.properties_in[t].pressure),
            )
            iscale.constraint_scaling_transform(
                self.fuel_inlet_translator.temperature_eqn[t],
                iscale.get_scaling_factor(self.fuel_inlet_translator.properties_in[t].temperature),
            )
            iscale.constraint_scaling_transform(
                self.fuel_inlet_translator.flow_mol_eqn[t],
                iscale.get_scaling_factor(self.fuel_inlet_translator.properties_in[t].flow_mol),
            )
            iscale.constraint_scaling_transform(
                self.fuel_outlet_translator.pressure_eqn[t],
                iscale.get_scaling_factor(self.fuel_outlet_translator.properties_in[t].pressure),
            )
            iscale.constraint_scaling_transform(
                self.fuel_outlet_translator.temperature_eqn[t],
                iscale.get_scaling_factor(self.fuel_outlet_translator.properties_in[t].temperature),
            )
            iscale.constraint_scaling_transform(
                self.fuel_outlet_translator.flow_mol_eqn[t],
                iscale.get_scaling_factor(self.fuel_outlet_translator.properties_in[t].flow_mol),
            )
            iscale.constraint_scaling_transform(self.cell_potential_eqn[t], 1)
            if not self.config.has_heat_transfer:
                iscale.constraint_scaling_transform(
                    self.heat_transfer_eqn[t], iscale.get_scaling_factor(self.heat[t])
                )

        for (t, i), c in self.air_translator.mole_frac_comp_eqn.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.air_translator.properties_out[t].mole_frac_comp[i]
                ),
            )
        for (t, i), c in self.fuel_inlet_translator.mole_frac_comp_eqn.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.fuel_inlet_translator.properties_out[t].mole_frac_comp[i]
                ),
            )
        for (t, i), c in self.fuel_outlet_translator.mole_frac_comp_eqn.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.fuel_outlet_translator.properties_out[t].mole_frac_comp[i]
                ),
            )



    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        init_log.info_high("SOFC Initialization Starting")
        solver_obj = get_solver(solver, optarg)

        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        self.air_side_inlet.fix()
        self.fuel_side_inlet.fix()
        self.air_side_outlet.unfix()
        self.fuel_side_outlet.unfix()

        self.fuel_inlet_translator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.strm1)

        self.sofc_reactor.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.strm2)

        self.air_separator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.strm3)

        self.fuel_outlet_translator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.strm4)

        self.air_translator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.strm5)

        self.air_mixer.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.strm6)

        self.sweep_heater.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        for t in self.current.keys():
            self.current[t] = pyo.value(self.current_expr[t])

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self, tee=slc.tee)
        if not pyo.check_optimal_termination(res):
            raise InitializationError("SOFC failed to initialize.")

        from_json(self, sd=istate, wts=sp)
        init_log.info_high("SOFC Initialization Complete")
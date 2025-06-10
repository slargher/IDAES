#################################################################################
# Simple 0D SOFC model for process design using IDAES IP
# Adapted style from SoecDesign (SOEC) to a Fuel Cell (SOFC)
#################################################################################
from pyomo.environ import Constraint, Var, Reference, units, Param, value
from pyomo.common.config import ConfigValue, ConfigBlock, In, Bool
from pyomo.network import Arc
import idaes.core.util.constants as iconst
from idaes.core import UnitModelBlockData, declare_process_block_class
from idaes.core.util.config import is_physical_parameter_block
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import GenericReactionParameterBlock
import idaes.models.unit_models as um
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.models_extra.power_generation.properties.natural_gas_PR import get_prop, get_rxn, EosType
from idaes.core.util.exceptions import ConfigurationError
import pyomo.environ as pyo
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from pyomo.util.infeasible import log_infeasible_constraints
import logging
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.model_statistics import (
    degrees_of_freedom, 
    report_statistics,
    activated_constraints_set,
    activated_equalities_set,
    unfixed_variables_set,
    fixed_variables_set,
    variables_set
)


@declare_process_block_class("SofcDesign", doc="Simple 0D SOFC model for process design.")
class SofcDesignData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    # Property packages for fuel (H2/H2O) and air (O2/N2)
    CONFIG.declare(
        "fuel_side_property_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package for the fuel side (H2/H2O)",
        ),
    )
    CONFIG.declare(
        "air_side_property_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package for the air side (O2/N2)",
        ),
    )
    CONFIG.declare(
        "reaction_eos",
        ConfigValue(
            default=EosType.PR,
            domain=In(EosType),
            description="EoS for reaction properties",
        ),
    )
    CONFIG.declare(
        "has_heat_transfer",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Include heat transfer?",
        ),
    )

    def build(self):
        super().build()
        # Validate components
        fcomps = set(self.config.fuel_side_property_package.component_list)
        if 'H2' not in fcomps:
            raise ConfigurationError("Fuel side must contain H2")
        
        acomps = set(self.config.air_side_property_package.component_list)
        if "O2" not in acomps or "N2" not in acomps:
            raise ConfigurationError("Air side must contain O2 and N2")

        # Build internal reaction properties
        self._add_reaction_properties()

        # Add translators, mixer, reactor
        self._add_unit_models()

        # Connect with arcs
        self._add_arcs()

        # Add variables and constraints
        self._add_variables()
        self._add_constraints()
        #self._add_air_separator_split_constraints()
        self._add_air_separator_split_constraints()
        self._add_fuel_in_trans_constraints()
        self._add_fuel_out_trans_constraints()
        self._add_o2rich_in_trans_constraints()
        self._add_catode_out_trans_constraints()

        # Expose ports
        self._add_ports()
        # Scaling
        self._scaling_guess()


    def _add_reaction_properties(self):
        # Generic property block for reaction (H2, H2O, O2)
        self.reaction_prop_params = GenericParameterBlock(
            **get_prop({"H2","H2O","O2"}, phases={"Vap"}, eos=self.config.reaction_eos)
        )
        # Reaction: H2 + 0.5 O2 -> H2O
        self.reaction_rxn_params = GenericReactionParameterBlock(
            **get_rxn(self.reaction_prop_params, {"h2_cmb"})
        )


    def _add_unit_models(self):

        # Translators for inlets
        self.fuel_in_trans = um.Translator(
            inlet_property_package=self.config.fuel_side_property_package,
            outlet_property_package=self.reaction_prop_params,
            outlet_state_defined=False,
        )

        self.air_separator = um.Separator(
            property_package=self.config.air_side_property_package,
            split_basis=um.SplittingType.componentFlow,  # <--- chiave
            outlet_list=["o2_poor_strm", "o2_rich_strm"],
        )

        self.o2rich_in_trans = um.Translator(
            inlet_property_package = self.config.air_side_property_package,
            outlet_property_package=self.reaction_prop_params,
            outlet_state_defined=False,
        )
        
        # Mixer to combine H2 and O2 streams
        self.mix = um.Mixer(
            inlet_list=["fuel_strm","o2_rich_strm"],
            property_package=self.reaction_prop_params,
            momentum_mixing_type=um.MomentumMixingType.none,
        )

        # Stoichiometric reactor simulating SOFC reaction
        self.fuel_cell_reactor = um.StoichiometricReactor(
            property_package=self.reaction_prop_params,
            reaction_package=self.reaction_rxn_params,
            has_heat_transfer=self.config.has_heat_transfer,
            has_pressure_change=False,
        )

        self.water_air_separator = um.Separator(
            property_package=self.reaction_prop_params,
            split_basis=um.SplittingType.componentFlow,
            outlet_list=['water_strm', 'o2_strm']
        )

        # Translators for outlets
        self.fuel_out_trans = um.Translator(
            inlet_property_package=self.reaction_prop_params,
            outlet_property_package=self.config.fuel_side_property_package,
            outlet_state_defined=False,
        )

        self.catode_out_trans = um.Translator(
            inlet_property_package=self.reaction_prop_params,
            outlet_property_package=self.config.air_side_property_package,
            outlet_state_defined=False,
        )

        self.mix_with_fake_air = um.Mixer(
            inlet_list=["o2_rich_strm", "fake_air_strm"],
            property_package=self.config.air_side_property_package,
            momentum_mixing_type=um.MomentumMixingType.none,
        )


       
    def _add_arcs(self):
        self.arc_fuel_to_mixer= Arc(
            source=self.fuel_in_trans.outlet,
            destination=self.mix.fuel_strm,
        )
        
        self.arc_o2RichSeparator_to_trans  = Arc(
            source=self.air_separator.o2_rich_strm,
            destination=self.o2rich_in_trans.inlet,
        )

        self.arc_o2RichTrans_to_mixer= Arc(
            source=self.o2rich_in_trans.outlet,
            destination=self.mix.o2_rich_strm,
        )

        self.arc_mixer_to_reactor = Arc(
            source=self.mix.outlet,
            destination=self.fuel_cell_reactor.inlet,
        )

        self.arc_reactor_to_separator = Arc(
            source=self.fuel_cell_reactor.outlet,
            destination=self.water_air_separator.inlet,
        )
        
        self.water_to_translator = Arc(
            source=self.water_air_separator.water_strm,
            destination=self.fuel_out_trans.inlet,
        )

        self.air_to_translator = Arc(
            source=self.water_air_separator.o2_strm,
            destination = self.catode_out_trans.inlet,
        )

        self.trans_to_mixer = Arc(
            source = self.catode_out_trans.outlet,
            destination = self.mix_with_fake_air.o2_rich_strm
        )

        self.separator_to_mixer = Arc(
            source=self.air_separator.o2_poor_strm,
            destination=self.mix_with_fake_air.fake_air_strm,
        )
        


        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

    def _add_variables(self):
        # Time index
        t = self.flowsheet().time

        # Current [A]
        self.current = Var(
            t,
            initialize=140,
            units=units.ampere,
            doc = 'Current'
        )

        # Fuel and air utilization
        self.fuel_util = Var(
            t,
            initialize=0.7, 
            units=units.dimensionless,
            doc = 'Fuel Utilizaiont Factor'
        )

        # Cell voltage
        self.cell_voltage = Var(
            t,
            initialize=1.3,
            units=units.volts,
            doc = 'Total Cell Voltage'
        )
        
        # Heat transfer
        self.thermal_power = Var(
            t, 
            initialize = 0,
            units=units.W,
            doc = 'Heat')
        
        self.number_of_cells = Var(
            domain = pyo.NonNegativeIntegers,
            units=units.dimensionless,
            doc="Numero di celle in stack"
        )
        self.cell_area = Var(
            domain=pyo.NonNegativeReals,
            units=units.m**2,
            doc="Area di una singola cella [m2]"
        )
        
        self.exchange_current_density = Param(initialize=2e4, mutable=True, units=units.ampere / units.m**2, doc="i0")
        self.activation_coeff = Param(initialize=0.05, mutable=True, units=units.volt, doc="Activation coefficient A")
        self.ohmic_resistance = Param(initialize=0.2, mutable=True, units=units.ohm * units.m**2, doc="Ohmic area resistance")


        self.hydrogen_side_outlet_temperature = pyo.Reference(
            self.fuel_cell_reactor.control_volume.properties_out[:].temperature
        )
        self.oxygen_side_outlet_temperature = pyo.Reference(
            self.mix_with_fake_air.mixed_state[:].temperature
        )

    def _add_constraints(self):
        
        @self.Constraint(self.flowsheet().time)
        def fuel_utilization_eqn(b, t):
            return ( b.fuel_cell_reactor.control_volume.properties_out[t].flow_mol_comp['H2']
                    == b.fuel_cell_reactor.control_volume.properties_in[t].flow_mol_comp['H2'] 
                    * (1 - b.fuel_util[t]) 
                    )
        
        @self.Expression(self.flowsheet().time)
        def current_expr(b, t):
            return (
                (b.fuel_cell_reactor.control_volume.properties_in[t].flow_mol_comp['H2']-
                b.fuel_cell_reactor.control_volume.properties_out[t].flow_mol_comp['H2'] )
                * 2 * iconst.Constants.faraday_constant / b.number_of_cells
            )
        
        @self.Constraint(self.flowsheet().time)
        def current_eqn(b, t):
            return b.current[t] == b.current_expr[t]
        
        @self.Expression(self.flowsheet().time)
        def power_expr(b, t):
            """
            Compute the reversible Nernst potential for a hydrogen-fueled SOFC at time index t,
            including an empirical linear temperature correction around 1073 K and the
            thermodynamic gas-composition adjustment.

            The reaction is:
                H2(g) + 1/2 O2(g) -> H2O(g)

            Returns:
                E_Nernst (V): Reversible (open-circuit) potential at time t
            """
            # 1. Read local operating conditions
            # Temperature at the reactor outlet (K)
            T = b.fuel_cell_reactor.control_volume.properties_out[t].temperature
            # Partial pressures (bar) from mole fraction * total pressure
            P_H2 = b.fuel_cell_reactor.control_volume.properties_in[t].mole_frac_comp["H2"]  * b.fuel_cell_reactor.control_volume.properties_in[t].pressure
            P_O2 = b.fuel_cell_reactor.control_volume.properties_in[t].mole_frac_comp["O2"]  * b.fuel_cell_reactor.control_volume.properties_in[t].pressure
            P_H2O = b.fuel_cell_reactor.control_volume.properties_out[t].mole_frac_comp["H2O"] * b.fuel_cell_reactor.control_volume.properties_out[t].pressure

            # 2. Empirical reference potential, linearized around T_ref = 1073 K:
            #    E0(T) ≈ E0(T_ref) + (dE0/dT) * (T - T_ref)
            # For H2/O2 -> H2O:
            #    dE0/dT = ΔS0 / (n F) ≈ -44.5 J/(mol·K) / (2 × 96485 C/mol) ≈ -0.0002 V/K
            #    E0(1073 K) ≈ 1.0 V (approximate standard reversible potential at 1073 K)
            E_empirical_ref = 1.0 - 0.0002 * (T - 1073)

            # 3. Thermodynamic Nernst correction for reactant/product activities:
            #    E_corr = (R T) / (n F) * ln(a_H2 * a_O2^0.5 / a_H2O)
            # We rewrite it as:
            #    E_corr = - (R T) / (2 F) * ln(P_H2O / (P_H2 * sqrt(P_O2)))
            R = iconst.Constants.gas_constant        # J/(mol·K)
            F = iconst.Constants.faraday_constant    # C/mol
            n = 2                          # electrons transferred per H2 -> H2O reaction
            nernst_correction = (R * T) / (n * F) * \
                                pyo.log(P_H2O / (P_H2 * pyo.sqrt(P_O2)))

            # 4. Combine to get the reversible (Nernst) potential
            #    E_Nernst = E_empirical_ref - E_corr
            return ((E_empirical_ref - nernst_correction) - b.cell_voltage[t]) * b.current[t] * b.number_of_cells

        @self.Constraint(self.flowsheet().time)
        def power_eqn(b, t):
            return b.thermal_power[t] == b.power_expr[t]
        
        @self.Constraint(self.flowsheet().time)
        def voltage_eqn(b, t):
            i = b.current[t] / b.cell_area  # [A/m2]

            # Nernst voltage già incluso in power_expr: E_Nernst = V + losses
            eta_act = b.activation_coeff * pyo.asinh(i / b.exchange_current_density)
            eta_ohm = i * b.ohmic_resistance

            V_loss = eta_act + eta_ohm
            E_nernst = b.power_expr[t] / (b.current[t] * b.number_of_cells) + b.cell_voltage[t]

            return b.cell_voltage[t] == E_nernst - V_loss

        for t in self.flowsheet().time:
            for comp in self.air_separator.mixed_state[t].flow_mol_comp:
                print(f"[DEBUG] Inlet flow of {comp}: {pyo.value(self.air_separator.mixed_state[t].flow_mol_comp[comp])}")

        
        # Fix the O2 separator split fractions to just remove O2 from fuel side
        self.water_air_separator.split_fraction[:, "o2_strm", "O2"].fix(1)
        self.water_air_separator.split_fraction[:, "o2_strm", "H2O"].fix(0)
        self.water_air_separator.split_fraction[:, "o2_strm", "H2"].fix(0)
        
    def _add_air_separator_split_constraints(self):
        self.air_separator.inlet.mole_frac_comp[0, "O2"].fix(0.21)
        self.air_separator.inlet.mole_frac_comp[0, "N2"].fix(0.79)

        # Fix split fractions: send all O2 to o2_outlet, all N2 to n2_outlet
        self.air_separator.split_fraction[0, "o2_rich_strm", "O2"].fix(1.0)
        self.air_separator.split_fraction[0, "o2_rich_strm", "N2"].fix(0.0)



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


    def _add_fuel_in_trans_constraints(self):
        t0 = self.flowsheet().time.first()

        print("[DEBUG][fuel_in_trans] IN components and values:")
        for k, v in self.fuel_in_trans.properties_in[t0].mole_frac_comp.items():
            print(f"   {k}: {pyo.value(v):.5g}")
        print("[DEBUG][fuel_in_trans] OUT components and values:")
        for k, v in self.fuel_in_trans.properties_out[t0].mole_frac_comp.items():
            print(f"   {k}: {pyo.value(v):.5g}")

        self.fuel_in_trans.inlet.mole_frac_comp[:, "H2"].fix(0.99)
        self.fuel_in_trans.inlet.mole_frac_comp[:, "H2O"].fix(0.01)

        comps = set(self.fuel_in_trans.properties_in[t0].mole_frac_comp.keys())

        @self.fuel_in_trans.Constraint(self.flowsheet().time, comps)
        def mole_frac_comp_eqn(b, t, c):
            return (
                b.properties_out[t].mole_frac_comp[c]
                == b.properties_in[t].mole_frac_comp[c]
            )

        self._translator_ftp_constraints(self.fuel_in_trans)


    def _add_fuel_out_trans_constraints(self):
        t0 = self.flowsheet().time.first()
        print("[DEBUG][fuel_out_trans] IN components and values:")
        for k, v in self.fuel_out_trans.properties_in[t0].mole_frac_comp.items():
            print(f"   {k}: {pyo.value(v):.5g}")
        print("[DEBUG][fuel_out_trans] OUT components and values:")
        for k, v in self.fuel_out_trans.properties_out[t0].mole_frac_comp.items():
            print(f"   {k}: {pyo.value(v):.5g}")

        comps = set(self.fuel_out_trans.properties_out[t0].mole_frac_comp.keys())
        comps.remove("H2")

        @self.fuel_out_trans.Constraint(self.flowsheet().time, comps)
        def mole_frac_comp_eqn(b, t, c):
            return (
                b.properties_out[t].mole_frac_comp[c]
                == b.properties_in[t].mole_frac_comp[c]
            )

        self._translator_ftp_constraints(self.fuel_out_trans)


    def _add_o2rich_in_trans_constraints(self):
        t0 = self.flowsheet().time.first()
        print("[DEBUG][o2rich_in_trans] IN components and values:")
        for k, v in self.o2rich_in_trans.properties_in[t0].mole_frac_comp.items():
            print(f"   {k}: {pyo.value(v):.5g}")
        print("[DEBUG][o2rich_in_trans] OUT components and values:")
        for k, v in self.o2rich_in_trans.properties_out[t0].mole_frac_comp.items():
            print(f"   {k}: {pyo.value(v):.5g}")

        comps = set(self.o2rich_in_trans.properties_out[t0].mole_frac_comp.keys())
        comps.remove("O2")

        @self.o2rich_in_trans.Constraint(self.flowsheet().time, comps)
        def mole_frac_zero_eqn(b, t, c):
            return b.properties_out[t].mole_frac_comp[c] == 1e-19

        for t in self.flowsheet().time:
            self.o2rich_in_trans.properties_out[t].mole_frac_comp["O2"].fix(1)

        self._translator_ftp_constraints(self.o2rich_in_trans)


    def _add_catode_out_trans_constraints(self):
        t0 = self.flowsheet().time.first()
        print("[DEBUG][catode_out_trans] IN components and values:")
        for k, v in self.catode_out_trans.properties_in[t0].mole_frac_comp.items():
            print(f"   {k}: {pyo.value(v):.5g}")
        print("[DEBUG][catode_out_trans] OUT components and values:")
        for k, v in self.catode_out_trans.properties_out[t0].mole_frac_comp.items():
            print(f"   {k}: {pyo.value(v):.5g}")

        comps = set(self.catode_out_trans.properties_out[t0].mole_frac_comp.keys())
        comps.remove("O2")

        @self.catode_out_trans.Constraint(self.flowsheet().time, comps)
        def mole_frac_comp_eqn(b, t, c):
            return b.properties_out[t].mole_frac_comp[c] == 1e-19

        for t in self.flowsheet().time:
            self.catode_out_trans.properties_out[t].mole_frac_comp["O2"].fix(1)

        self._translator_ftp_constraints(self.catode_out_trans)


    def _add_ports(self):
        self.add_inlet_port(
            name="fuel_side_inlet",
            block=self.fuel_in_trans.properties_in,
        )

        self.add_inlet_port(
            name="air_side_inlet",
            block=self.air_separator.mixed_state,
        )

        self.add_outlet_port(
            name="fuel_side_outlet",
            block=self.fuel_out_trans.properties_out,
        )

        self.add_outlet_port(
            name="air_side_outlet",
            block=self.mix_with_fake_air.mixed_state,
        )

    def _scaling_guess(self,
                    flow_scale_fuel=1e-4,
                    flow_scale_air=1e-3,
                    enth_scale=1e-4,
                    current_scale=1e2,
                    voltage_scale=1,
                    heat_scale=1e0,
                    mole_frac_scale=10):
        
        """
        Standard scaling guess routine for SOFC/SOEC units.

        Parameters:
            flow_scale_fuel (float): scaling for fuel molar flow (mol/s)
            flow_scale_air (float): scaling for air molar flow (mol/s)
            enth_scale (float): scaling for molar enthalpy [J/mol]
            current_scale (float): scaling for current [A]
            voltage_scale (float): scaling for cell voltage [V]
            heat_scale (float): scaling for heat [W]
            mole_frac_scale (float): scaling for mole fractions (dimensionless)
        """

        # Set default flow mol scaling
        self.config.fuel_side_property_package.set_default_scaling("flow_mol", flow_scale_fuel)
        self.config.air_side_property_package.set_default_scaling("flow_mol", flow_scale_air)
        self.reaction_prop_params.set_default_scaling("flow_mol", flow_scale_fuel)

        self.config.fuel_side_property_package.set_default_scaling("flow_mol_phase", flow_scale_fuel)
        self.config.air_side_property_package.set_default_scaling("flow_mol_phase", flow_scale_air)
        self.reaction_prop_params.set_default_scaling("flow_mol_phase", flow_scale_fuel)

        # Scaling for mole fractions
        self.config.fuel_side_property_package.set_default_scaling("mole_frac_comp", mole_frac_scale)
        self.config.air_side_property_package.set_default_scaling("mole_frac_comp", mole_frac_scale)
        self.reaction_prop_params.set_default_scaling("mole_frac_comp", mole_frac_scale)
        self.reaction_prop_params.set_default_scaling("mole_frac_phase_comp", mole_frac_scale)

        # Scaling for key variables
        for v in self.current.values():
            iscale.set_scaling_factor(v, current_scale)
        for v in self.cell_voltage.values():
            iscale.set_scaling_factor(v, voltage_scale)
        for v in self.thermal_power.values():
            iscale.set_scaling_factor(v, heat_scale)
        for v in self.fuel_util.values():
            iscale.set_scaling_factor(v, 1)

        # Set enthalpy scaling
        t0 = self.flowsheet().time.first()
        try:
            iscale.set_scaling_factor(self.fuel_cell_reactor.control_volume.properties_in[t0].enth_mol_phase["Vap"], enth_scale)
            iscale.set_scaling_factor(self.fuel_cell_reactor.control_volume.properties_out[t0].enth_mol_phase["Vap"], enth_scale)
        except AttributeError:
            pass  # Optional: not all models have enth_mol_phase

        # Apply additional scaling if separator/mixer used
        if hasattr(self, "air_separator"):
            try:
                iscale.set_scaling_factor(self.air_separator.o2_poor_strm_state[t0].enth_mol_phase["Vap"], enth_scale)
            except AttributeError:
                pass

        if hasattr(self, "mix"):
            try:
                iscale.set_scaling_factor(self.mix.outlet_state[t0].enth_mol_phase["Vap"], enth_scale)
            except AttributeError:
                pass

    def apply_constraint_scaling(self):
        """
        Apply appropriate scaling factors to each key constraint based on
        previously scaled variables.
        """
        tset = self.flowsheet().time

        # Utility function to apply scaling safely
        def scale_constraint(constraint, var, default=1):
            try:
                sf = iscale.get_scaling_factor(var, default)
                iscale.constraint_scaling_transform(constraint, sf)
            except Exception:
                pass

        for t in tset:
            # Main equation constraints
            scale_constraint(self.fuel_utilization_eqn[t], self.fuel_util[t])
            scale_constraint(self.current_eqn[t], self.current[t])
            scale_constraint(self.voltage_eqn[t], self.cell_voltage[t])
            scale_constraint(self.power_eqn[t], self.thermal_power[t])
            scale_constraint(self.cell_potential_eqn[t], self.cell_voltage[t])

            # Translator: fuel_in_trans
            scale_constraint(self.fuel_in_trans.flow_mol_eqn[t], self.fuel_in_trans.properties_in[t].flow_mol)
            scale_constraint(self.fuel_in_trans.temperature_eqn[t], self.fuel_in_trans.properties_in[t].temperature)
            scale_constraint(self.fuel_in_trans.pressure_eqn[t], self.fuel_in_trans.properties_in[t].pressure)

            # Translator: o2rich_in_trans
            scale_constraint(self.o2rich_in_trans.flow_mol_eqn[t], self.o2rich_in_trans.properties_in[t].flow_mol)
            scale_constraint(self.o2rich_in_trans.temperature_eqn[t], self.o2rich_in_trans.properties_in[t].temperature)
            scale_constraint(self.o2rich_in_trans.pressure_eqn[t], self.o2rich_in_trans.properties_in[t].pressure)

            # Translator: fuel_out_trans
            scale_constraint(self.fuel_out_trans.flow_mol_eqn[t], self.fuel_out_trans.properties_in[t].flow_mol)
            scale_constraint(self.fuel_out_trans.temperature_eqn[t], self.fuel_out_trans.properties_in[t].temperature)
            scale_constraint(self.fuel_out_trans.pressure_eqn[t], self.fuel_out_trans.properties_in[t].pressure)

        # Mole fraction equations
        for (t, j), con in self.fuel_in_trans.mole_frac_comp_eqn.items():
            var = self.fuel_in_trans.properties_out[t].mole_frac_comp[j]
            sf = iscale.get_scaling_factor(var, default=10)
            iscale.constraint_scaling_transform(con, sf)

        for (t, j), con in self.o2rich_in_trans.mole_frac_comp_eqn.items():
            var = self.o2rich_in_trans.properties_out[t].mole_frac_comp[j]
            sf = iscale.get_scaling_factor(var, default=10)
            iscale.constraint_scaling_transform(con, sf)

        for (t, j), con in self.fuel_out_trans.mole_frac_comp_eqn.items():
            var = self.fuel_out_trans.properties_out[t].mole_frac_comp[j]
            sf = iscale.get_scaling_factor(var, default=10)
            iscale.constraint_scaling_transform(con, sf)

    def initialize_build(self, outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        """
        Generalized initializer for SOFC model, initializing internal units in logical order,
        propagating states, solving the unit, and restoring input states.
        """
        logging.getLogger("pyomo.core").setLevel(logging.INFO)
        logging.getLogger("pyomo.util.infeasible").setLevel(logging.INFO)
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        init_log.info_high("SOFC Initialization Starting")
        solver_obj = get_solver(solver, optarg)

        # Save current fixed state to restore later
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # Fix inputs, unfix outputs
        self.fuel_side_inlet.fix()
        self.air_side_inlet.fix()
        self.fuel_util.fix()
        self.fuel_side_outlet.unfix()
        self.air_side_outlet.unfix()

        # Initialize unit models in order and propagate states
        self.fuel_in_trans.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.arc_fuel_to_mixer)

        self.air_separator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.arc_o2RichSeparator_to_trans)

        self.o2rich_in_trans.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.arc_o2RichTrans_to_mixer)

        self.mix.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.arc_mixer_to_reactor)

        self.fuel_cell_reactor.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.arc_reactor_to_separator)

        self.water_air_separator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.water_to_translator)
        propagate_state(self.air_to_translator)

        self.fuel_out_trans.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        self.catode_out_trans.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        self.mix_with_fake_air.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.trans_to_mixer)
        propagate_state(self.separator_to_mixer)



        # Compute initial current estimate from flow diff
        for t in self.current:
            self.current[t] = pyo.value(self.current_expr[t])
            print("Initial current guess:", pyo.value(self.current[self.flowsheet().time.first()]))



        dof = degrees_of_freedom(self)
        print("\n[DEBUG] DOF =", dof)

        from pyomo.core.expr import current as EXPR
        from pyomo.common.collections import ComponentSet

        print("\n[DEBUG] Variabili non fissate:")
        for v in self.component_data_objects(ctype=pyo.Var, active=True, descend_into=True):
            if not v.fixed:
                try:
                    print(f" - {v.name} = {pyo.value(v):.4g}")
                except:
                    print(f" - {v.name} = [None]")

        print("\n[DEBUG] Vincoli attivi:")
        for c in self.component_data_objects(ctype=pyo.Constraint, active=True, descend_into=True):
            try:
                print(f" - {c.name} = {pyo.value(c):.4g}")

            except:
                print(f" - {c.name} = [None]")
            

        # Solve the entire block
        # with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
        #     res = solver_obj.solve(self, tee=slc.tee)
        res = solver_obj.solve(self, tee=True)
        log_infeasible_constraints(self, log_expression=True, log_variables=True)
        print(res)

        if not pyo.check_optimal_termination(res):
            raise InitializationError("SOFC failed to initialize.")

        # Restore previous input states
        from_json(self, sd=istate, wts=sp)
        init_log.info_high("SOFC Initialization Complete")

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
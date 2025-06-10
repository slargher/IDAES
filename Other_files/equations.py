#Heater new unit     

@m.fs.Constraint()
    def sofc_power_dc_constraint(fs):
        return fs.sofc_power_dc == fs.stack_current * fs.sofc.stack_voltage # from excel, multiply * number of cells

# gives the outlet heat   

    @m.fs.Constraint()
    def SOFC_energy_balance(fs):
        return -1 * pyunits.convert(
            fs.anode.heat_duty[0], pyunits.MW
        ) == fs.sofc_power_dc + pyunits.convert(
            fs.cathode_heat.heat_duty[0], pyunits.MW
        )


# Constants
F = 96485.3329  # Faraday's constant, C/mol
FU = 0.77  # Fuel utilization factor
n_H2_in = 0.075397  # mol/s
A_cell_cm2 = 320  # cm^2
A_cell_m2 = A_cell_cm2 * 1e-4  # convert to m^2
N_cells = 80

T_fuel_in = 668.9 + 273.15  # = 942.05 K
T_fuel_out = 740.0 + 273.15  # = 1013.15 K
T_tot_in = 675.3 + 273.15  # = 948.45 K
T_air_in = 675.3 + 273.15  # = 948.45 K
T_air_out = 740.0 + 273.15  # = 1013.15 K
T_tot_out = 740.0 + 273.15  # = 1013.15 K

# If you want, you can also store the summed values directly:
T_fuel_in = 942.05
T_fuel_out = 1013.15
T_tot_in = 948.45
T_air_in = 948.45
T_air_out = 1013.15
T_tot_out = 1013.15

# Output variables
j_model_t
V
Lambda = n_O2_in / (0.5 * n_H2_in * FU)  # Lambda = ṅ_O2,in(t) / 0.5 * ṅ_H2,in(t) * FU

# Equations

# 1. Fuel Utilization Constraint
n_H2_out = n_H2_in * (1 - FU)  # ṅ_H2,out(t) = ṅ_H2,in(t) * (1 - FU(t))

# 2. Current (j) Expression
# Assuming A_attiva is same as A_cell
A_attiva = A_cell_m2
j_t = ((n_H2_in - n_H2_out) * 2 * F * N_cells) / A_attiva  # j(t)

# 3. Current Equality Constraint
j_model_t = j_t  # j_model(t) = j(t)

# 4. Voltage Polynomial Fit Equation
# V = V(j) = -434.47 * I^5 + 492.56 * I^4 - 211.6 * I^3 + 44.464 * I^2 - 5.4157 * I + 1.168
V = -434.47 * j_model_t**5 + 492.56 * j_model_t**4 - 211.6 * j_model_t**3 + 44.464 * j_model_t**2 - 5.4157 * I + 1.168

# 5. Electrical Power
P = V * j_model_t * A_cell_cm2 * N_cells

# Thermal Balance
Q_rxn = (+ (dHstd_T25C + Cp_avg * dT) * (n_H2_in_t * FU))

# Cp_avg definition
Cp_avg = (cp_avg_H2 * x_H2) - (cp_avg_H2O * x_H2O) - 0.5 * (cp_avg_O2 * x_O2)

# cp at average temperature
cp_avg_x = cp((T1 + T2) / 2)



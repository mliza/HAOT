'''
    Date:   08/27/2023
    Author: Martin E. Liza
    File:   aerodynamics.py
    Def:    Contains aerodynamics helper functions. 
'''
import os 
import molmass
import numpy as np 
import scipy.constants as s_consts
from haot.src import constants_tables

def sutherland_law_viscoity(temperature_K, molecule='Air'):
    """
    Calculates the Sutherland's law

    Parameters:
        temperature_K (float): reference temperature
        molecule (str) : Air (default), Argon, N2, O2 
    
    Returns:
        dynamic_viscosity (float): [kg/ms]

    """
    const = constants_tables.sutherland_constants(molecule)

    dynamic_viscosity = const['temperature_ref'] + const['sutherland_visc']
    dynamic_viscosity /= (temperature_K + const['sutherland_visc'])
    dynamic_viscosity *= (temperature_K / const['temperature_ref'])**(3/2)

    return const['viscosity_ref'] * dynamic_viscosity # [kg/ms]

def sutherland_law_conductivity(temperature_K, molecule='Air'):
    """
    Calculates the Sutherland's law

    Parameters:
        temperature_K (float): reference temperature
        molecule (str) : Air (default), Argon, N2, O2 
    
    Returns:
        thermal_conductivity (float): [W/mK]

    """
    const = constants_tables.sutherland_constants(molecule)
    thermal_conductivity = const['sutherland_cond']
    thermal_conductivity += const['temperature_ref']
    thermal_conductivity /= (temperature_K + const['sutherland_cond'])
    thermal_conductivity *= (temperature_K / const['temperature_ref'])**(3/2)

    return const['conductivity_ref'] * thermal_conductivity # [W/mK]


# Air atomic mass
def air_atomic_mass():
    molecules = ['N+', 'O+', 'NO+', 'N2+', 'O2+', 'N', 'O', 'NO', 'N2', 'O2'] 
    air_atomic_dict = { }
    for i in molecules:
        air_atomic_dict[i] = molmass.Formula(i).mass

    return air_atomic_dict #[g/mol]


# Speed of sound 
def speed_of_sound(temperature_K, adiabatic_indx=1.4):
    gas_const          = s_consts.R                         # [J/mol*K]
    air_atomic_mass    = air_atomic_mass()                  # [g/mol]
    air_molecular_mass = (0.7803 * air_atomic_mass['N2'] +  # [kg/mol]
                          0.2099 * air_atomic_mass['O2'] + 
                          0.0003 * air_atomic_mass['CO2']) * 1E-3 
    spd_of_sound       = np.sqrt( adiabatic_indx * temperature_K * 
                                  gas_const / air_molecular_mass )  
    return spd_of_sound # [m/s]


# Normal shock relations  
def normal_shock_relations(mach_1, adiabatic_indx=1.4):
    # REF: https://www.grc.nasa.gov/www/k-12/airplane/normal.html
    # NOTE: var_r = var_1 / var_2 = var_preshock / var_postshock = [ ] 
    gamma_minus   = adiabatic_indx - 1
    gamma_plus    = adiabatic_indx + 1
    mach_11       = mach_1 ** 2
    mach_2        = np.sqrt( (gamma_minus * mach_11 + 2) /
                             (2 * adiabatic_indx * mach_11 - gamma_minus) )
    pressure_r    = (2 * adiabatic_indx * mach_11 - gamma_minus) / gamma_plus
    temperature_r = ( (2 * adiabatic_indx * mach_11 - gamma_minus) *
                      (gamma_minus * mach_11 + 2) / (gamma_plus**2 * mach_11) )
    density_r     = gamma_plus * mach_11 / (gamma_minus * mach_11 + 2)

    # Return Dictionary 
    normal_shock_dict = { 'mach_2'        : mach_2,
                          'pressure_r'    : pressure_r, 
                          'temperature_r' : temperature_r,
                          'density_r'     : density_r }
    return normal_shock_dict # [ ]

# Oblique shock relations
def oblique_shock_relations(mach_1, shock_angle_deg, adiabatic_indx=1.4):
# REF : Modern Compressible Flows With Historical Ref., eq 4.7 - 4.11 
# NOTE: Equations only work for weak shocks 
# Note ratio = var_1 / var_2
    shock_angle = np.radians(shock_angle_deg)  # radians
    sin2_shock  = np.sin(shock_angle)**2
    mach_n1     = mach_1 * np.sin(shock_angle) # normal mach number 
    mach_n11    = mach_n1 ** 2  # normal mach number square 
    # Calculates Deflection angle (Eq. 4.17) 
    tan_deflection_ang   = ( (2 / np.tan(shock_angle)) *  ( (mach_n11 - 1) / 
                           (mach_1**2 * (adiabatic_indx + 
                           np.cos(2 * shock_angle)) + 2) ) ) 
    deflection_angle_deg = np.degrees(np.arctan(1 / tan_deflection_ang)) 
    # Calculates properties downstream the shock  
    density_r     = ( ((adiabatic_indx + 1) * mach_n1**2) / 
                     ((adiabatic_indx - 1) * mach_n1**2 + 2) )
    pressure_r    = (1 + 2 * adiabatic_indx * (mach_n1**2 - 1) / 
                     (adiabatic_indx + 1) )
    temperature_r = pressure_r * 1 / density_r 
    # Calculates mach 2 
    mach_n2 = np.sqrt( (mach_n1**2 + (2 / (adiabatic_indx - 1))) / 
            ((2 * adiabatic_indx / (adiabatic_indx - 1)) * mach_n1**2 - 1) )
    mach_2  = mach_n2 / np.sin(np.radians(shock_angle_deg - 
                              deflection_angle_deg)) 
    # Dictionary 
    oblique_shock_dict = { 'mach_2'                : mach_2,
                           'pressure_r'            : pressure_r, 
                           'temperature_r'         : temperature_r,
                           'density_r'             : density_r, 
                           'deflection_angle_degs' : deflection_angle_deg } 
    return oblique_shock_dict



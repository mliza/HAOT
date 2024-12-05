"""
    Date:   08/27/2023
    Author: Martin E. Liza
    File:   aerodynamics.py
    Def:    Contains aerodynamics helper functions.
"""

import molmass
import numpy as np
import scipy.constants as s_consts
from haot import constants_tables


def sutherland_law_viscosity(temperature_K, molecule="Air"):
    """
    Calculates the Sutherland's law of viscosity

    Parameters:
        temperature_K (float): reference temperature
        molecule (str): Air (default), Argon, N2, O2

    Returns:
        dynamic_viscosity (float): [kg/ms]

    """
    const = constants_tables.sutherland_constants(molecule)

    dynamic_viscosity = const["temperature_ref"] + const["sutherland_visc"]
    dynamic_viscosity /= temperature_K + const["sutherland_visc"]
    dynamic_viscosity *= (temperature_K / const["temperature_ref"]) ** (3 / 2)

    return const["viscosity_ref"] * dynamic_viscosity  # [kg/ms]


def sutherland_law_conductivity(temperature_K, molecule="Air"):
    """
    Calculates the Sutherland's law of thermal conductivity

    Parameters:
        temperature_K (float): reference temperature
        molecule (str): Air (default), Argon, N2, O2

    Returns:
        thermal_conductivity (float): [W/mK]

    """
    const = constants_tables.sutherland_constants(molecule)
    thermal_conductivity = const["sutherland_cond"]
    thermal_conductivity += const["temperature_ref"]
    thermal_conductivity /= temperature_K + const["sutherland_cond"]
    thermal_conductivity *= temperature_K / const["temperature_ref"]
    thermal_conductivity **= 3 / 2

    return const["conductivity_ref"] * thermal_conductivity  # [W/mK]


def air_atomic_molar_mass():
    """Returns the atomic molar mass of an 11-species air in [g/mol]"""
    molecules = ["N+", "O+", "NO+", "N2+", "O2+", "N", "O", "NO", "N2", "O2"]
    air_atomic_dict = {}
    for i in molecules:
        air_atomic_dict[i] = molmass.Formula(i).mass

    return air_atomic_dict  # [g/mol]


def speed_of_sound(temperature_K, adiabatic_indx=1.4):
    """
    Calculates the speed of sound

    Parameters:
        temperature_K (float): reference temperature
        adiabatic_indx (double): adiabatic index, 1.4 (default)

    Returns:
        spd_of_sound (float): speed of sound in [m/s]

    """
    gas_const = s_consts.R  # [J/mol*K]
    air_atomic_mass = air_atomic_molar_mass()  # [g/mol]
    air_molecular_mass = (
        0.7803 * air_atomic_mass["N2"]  # [kg/mol]
        + 0.2099 * air_atomic_mass["O2"]
        + 0.0003 * air_atomic_mass["CO2"]
    ) * 1e-3
    spd_of_sound = np.sqrt(
        adiabatic_indx * temperature_K * gas_const / air_molecular_mass
    )
    return spd_of_sound  # [m/s]


def normal_shock_relations(mach_1, adiabatic_indx=1.4):
    """
    Calculates normal shock relations

    Reference:
        https://www.grc.nasa.gov/www/k-12/airplane/normal.html

    Parameters:
        mach_1 (float): pre-shock mach number
        adiabatic_indx (double): adiabatic index, 1.4 (default)

    Returns:
        normal_shock_dic: dictionary with normal shock properties
            mach_2 (float): post-shock mach number [ ]
            pressure_r (float): pressure ratio (post-shock / pre-shock) [ ]
            temperature_r (float): temperature ratio (post-shock / pre-shock) [ ]
            density_r (float): density ratio (post-shock / pre-shock) [ ]
            pressure_tr (float): stagnation pressure ratio (post-shock / pre-shock) [ ]
            temperature_tr (float): stagnation temperature ratio (post-shock / pre-shock) [ ]

    """
    gamma_minus = adiabatic_indx - 1
    gamma_plus = adiabatic_indx + 1
    mach_11 = mach_1**2
    mach_2 = gamma_minus * mach_11 + 2
    mach_2 /= 2 * adiabatic_indx * mach_11 - gamma_minus
    mach_2 **= 0.5

    pressure_r = (2 * adiabatic_indx * mach_11 - gamma_minus) / gamma_plus

    temperature_r = 2 * adiabatic_indx * mach_11 - gamma_minus
    temperature_r *= gamma_minus * mach_11 + 2
    temperature_r /= gamma_plus**2 * mach_11

    density_r = gamma_plus * mach_11 / (gamma_minus * mach_11 + 2)

    pressure_tr1 = gamma_plus / (2 * adiabatic_indx * mach_11 - gamma_minus)
    pressure_tr1 **= 1 / gamma_minus
    pressure_tr2 = gamma_plus * mach_11 / (gamma_minus * mach_11 + 2)
    pressure_tr2 **= adiabatic_indx / gamma_minus
    pressure_tr = pressure_tr1 * pressure_tr2

    # Return Dictionary
    normal_shock_dict = {
        "mach_2": mach_2,
        "pressure_r": pressure_r,
        "temperature_r": temperature_r,
        "density_r": density_r,
        "pressure_tr": pressure_tr,
        "temperature_tr": 1.0,
    }
    return normal_shock_dict  # [ ]


# Oblique shock relations
# TODO: Update this doc string on this
def oblique_shock_relations(mach_1, shock_angle_deg, adiabatic_indx=1.4):
    # REF : Modern Compressible Flows With Historical Ref., eq 4.7 - 4.11
    # NOTE: Equations only work for weak shocks
    # Note ratio = var_1 / var_2
    shock_angle = np.radians(shock_angle_deg)  # radians
    mach_n1 = mach_1 * np.sin(shock_angle)  # normal mach number
    mach_n11 = mach_n1**2  # normal mach number square
    # Calculates Deflection angle (Eq. 4.17)
    tan_deflection_ang = (2 / np.tan(shock_angle)) * (
        (mach_n11 - 1) / (mach_1**2 * (adiabatic_indx + np.cos(2 * shock_angle)) + 2)
    )
    deflection_angle_deg = np.degrees(np.arctan(1 / tan_deflection_ang))
    # Calculates properties downstream the shock
    density_r = ((adiabatic_indx + 1) * mach_n1**2) / (
        (adiabatic_indx - 1) * mach_n1**2 + 2
    )
    pressure_r = 1 + 2 * adiabatic_indx * (mach_n1**2 - 1) / (adiabatic_indx + 1)
    temperature_r = pressure_r * 1 / density_r
    # Calculates mach 2
    mach_n2 = np.sqrt(
        (mach_n1**2 + (2 / (adiabatic_indx - 1)))
        / ((2 * adiabatic_indx / (adiabatic_indx - 1)) * mach_n1**2 - 1)
    )
    mach_2 = mach_n2 / np.sin(np.radians(shock_angle_deg - deflection_angle_deg))
    # Dictionary
    oblique_shock_dict = {
        "mach_2": mach_2,
        "pressure_r": pressure_r,
        "temperature_r": temperature_r,
        "density_r": density_r,
        "deflection_angle_degs": deflection_angle_deg,
    }
    return oblique_shock_dict

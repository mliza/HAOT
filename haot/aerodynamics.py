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


def sutherland_law_viscosity(temperature_K: float, molecule: str = "Air") -> float:
    """
    Calculates the Sutherland's law of viscosity

    Parameters:
        temperature_K: reference temperature in [K]
        molecule: Air (default), Argon, N2, O2

    Returns:
        dynamic viscosity in [kg/ms]

    Examples:
        >> sutherland_law_viscosity(300.0)

    """
    const = constants_tables.sutherland_constants(molecule)

    dynamic_viscosity = const["temperature_ref"] + const["sutherland_visc"]
    dynamic_viscosity /= temperature_K + const["sutherland_visc"]
    dynamic_viscosity *= (temperature_K / const["temperature_ref"]) ** (3 / 2)

    return const["viscosity_ref"] * dynamic_viscosity  # [kg/ms]


def sutherland_law_conductivity(temperature_K: float, molecule: str = "Air") -> float:
    """
    Calculates the Sutherland's law of thermal conductivity

    Parameters:
        temperature_K: reference temperature in [K]
        molecule: Air (default), Argon, N2, O2

    Returns:
        thermal conductivity in [W/mK]

    Examples:
        >> sutherland_law_conductivity(300.0)

    """
    const = constants_tables.sutherland_constants(molecule)
    thermal_conductivity = const["sutherland_cond"]
    thermal_conductivity += const["temperature_ref"]
    thermal_conductivity /= temperature_K + const["sutherland_cond"]
    thermal_conductivity *= temperature_K / const["temperature_ref"]
    thermal_conductivity **= 3 / 2

    return const["conductivity_ref"] * thermal_conductivity  # [W/mK]


def air_atomic_molar_mass(molecules: str = None) -> dict[str, float]:
    """
    Returns the atomic molar mass

    Paremeters:
        molecule: Molecules that need the molar mass (11 species air is the default)

    Returns:
        species in [g/mol]

    Examples:
        >> air_atomic_molar_mass()

    """
    if not molecules:
        molecules = ["N+", "O+", "NO+", "N2+", "O2+", "N", "O", "NO", "N2", "O2"]
        
    air_atomic_dict = {i: molmass.Formula(i).mass for i in molecules}
    
    return air_atomic_dict  # [g/mol]


def speed_of_sound(temperature_K: float, adiabatic_indx: float = 1.4) -> float:
    """
    Calculates the speed of sound

    Parameters:
        temperature_K: reference temperature in [K]
        adiabatic_indx: adiabatic index, 1.4 (default)

    Returns:
        speed of sound in [m/s]

    Examples:
        >> speed_of_sound(300.0)

    """
    gas_const = s_consts.R  # [J/mol*K]
    air_atomic_mass = air_atomic_molar_mass(["N2", "O2", "Ar", "CO2"])  # [g/mol]

    air_molecular_mass = (
        78 * air_atomic_mass["N2"]
        + 21 * air_atomic_mass["O2"]
        + 0.93 * air_atomic_mass["Ar"]
        + 0.07 * air_atomic_mass["CO2"]
    ) * 1e-5 # [kg/mol]
    spd_of_sound = np.sqrt(
        adiabatic_indx * temperature_K * gas_const / air_molecular_mass
    )
    return spd_of_sound  # [m/s]


def normal_shock_relations(mach_1: float, adiabatic_indx: float = 1.4) -> dict[str, float]:
    """
    Calculates normal shock relations

    Parameters:
        mach_1: pre-shock mach number
        adiabatic_indx: adiabatic index, 1.4 (default)

    Returns:
        dict: A dictionary containing:
            - mach_2: post-shock mach number
            - pressure_r: pressure ratio (post-shock / pre-shock)
            - temperature_r: temperature ratio (post-shock / pre-shock)
            - density_r: density ratio (post-shock / pre-shock)
            - pressure_tr: stagnation pressure ratio (post-shock / pre-shock)
            - temperature_tr: stagnation temperature ratio (post-shock / pre-shock)
        
    Reference:
        Normal Shock Wave - NASA (https://www.grc.nasa.gov/www/k-12/airplane/normal.html)
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

"""
Date:   08/27/2023
Author: Martin E. Liza
File:   aerodynamics.py
Def:    Contains aerodynamics helper functions.
"""

import molmass
import numpy as np
import scipy.constants as s_consts
from haot import constants as constants_tables


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

    Reference:
        Viscous Fluid Flow, International Edition, 4th (White F., ISBN 978 1 260 59786)
    """
    # Checking cases
    if type(temperature_K) is float and temperature_K < 0:
        raise ValueError("Temperature must be greater than 0 Kelvin!")
    if type(temperature_K) is np.ndarray and (temperature_K < 0).any():
        raise ValueError("Temperature must be greater than 0 Kelvin!")
    if molecule not in ["Air", "Argon", "N2", "O2"]:
        raise ValueError("This function only supports Air, Argon, N2 or O2")

    const = constants_tables.sutherland_constants(molecule)

    # Eq 1-34
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

    Reference:
        Viscous Fluid Flow, International Edition, 4th (White F., ISBN 978 1 260 59786)
    """
    # Checking cases
    if type(temperature_K) is float and temperature_K < 0:
        raise ValueError("Temperature must be greater than 0 Kelvin!")
    if type(temperature_K) is np.ndarray and (temperature_K < 0).any():
        raise ValueError("Temperature must be greater than 0 Kelvin!")
    if molecule not in ["Air", "Argon", "N2", "O2"]:
        raise ValueError("This function only supports Air, Argon, N2 or O2")

    const = constants_tables.sutherland_constants(molecule)

    # Eq 1-41b
    thermal_conductivity = const["sutherland_cond"]
    thermal_conductivity += const["temperature_ref"]
    thermal_conductivity /= temperature_K + const["sutherland_cond"]
    thermal_conductivity *= temperature_K / const["temperature_ref"]
    thermal_conductivity **= 3 / 2

    return const["conductivity_ref"] * thermal_conductivity  # [W/mK]


def air_atomic_molar_mass(molecules: str = None) -> dict[str, float]:
    """
    Returns the atomic molar mass

    Parameters:
        molecule: Molecules that need the molar mass (11 species air is the default)

    Returns:
        species in [g/mol]

    Examples:
        >> air_atomic_molar_mass(["N+", "N2"])
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
    # Checking cases
    if type(temperature_K) is float and temperature_K < 0:
        raise ValueError("Temperature must be greater than 0 Kelvin!")
    if type(temperature_K) is np.ndarray and (temperature_K < 0).any():
        raise ValueError("Temperature must be greater than 0 Kelvin!")

    gas_const = s_consts.R  # [J/mol*K]
    air_atomic_mass = air_atomic_molar_mass(["N2", "O2", "Ar", "CO2"])  # [g/mol]

    air_molecular_mass = (
        78 * air_atomic_mass["N2"]
        + 21 * air_atomic_mass["O2"]
        + 0.93 * air_atomic_mass["Ar"]
        + 0.07 * air_atomic_mass["CO2"]
    ) * 1e-5  # [kg/mol]
    spd_of_sound = np.sqrt(
        adiabatic_indx * temperature_K * gas_const / air_molecular_mass
    )
    return spd_of_sound  # [m/s]


def isentropic_relations(
    mach_1: float, adiabatic_indx: float = 1.4
) -> dict[str, float]:
    """
    Calculates isentropic relations

    Parameters:
        mach_1: pre-shock mach number

    Returns:
        dict: A dictionary containing:
            - pressure_s: pressure ratio (post-shock stagnation / pre-shock)
            - temperature_s: temperature ratio (post-shock stagnation / pre-shock)
            - density_s: density ratio (post-shock stagnation / pre-shock)
    Examples:
        >> isentropic_relations(3.0)

    Reference:
        Modern Compressible Flow With Historic Perspective, International
        Edition 4th (Anderson J., ISBN 978 1 260 57082 3)
    """
    # Checking cases
    if mach_1 <= 0:
        raise ValueError("Mach number has to be greater than 0!")
    gamma_minus = adiabatic_indx - 1
    gamma_ratio = gamma_minus / 2

    # Stagnation temperature (Eq. 3.28)
    temperature_s = 1 + gamma_ratio * mach_1**2

    # Stagnation pressure (Eq. 3.29)
    pressure_s = temperature_s ** (adiabatic_indx / gamma_minus)

    # Stagnation density (Eq. 3.30)
    density_s = temperature_s ** (1.0 / gamma_minus)

    isentropic_dict = {
        "pressure_s": pressure_s,
        "temperature_s": temperature_s,
        "density_s": density_s,
    }
    return isentropic_dict


def normal_shock_relations(
    mach_1: float, adiabatic_indx: float = 1.4
) -> dict[str, float]:
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
            - pressure_s: stagnation pressure ratio (post-shock / pre-shock)

    Reference:
        Normal Shock Wave - NASA (https://www.grc.nasa.gov/www/k-12/airplane/normal.html)
    """
    gamma_minus = adiabatic_indx - 1
    gamma_plus = adiabatic_indx + 1
    mach_11 = mach_1**2

    # Mach post-shock
    mach_2 = gamma_minus * mach_11 + 2
    mach_2 /= 2 * adiabatic_indx * mach_11 - gamma_minus
    mach_2 **= 0.5

    # Pressure ratio
    pressure_r = (2 * adiabatic_indx * mach_11 - gamma_minus) / gamma_plus

    # Temperature ratio
    temperature_r = 2 * adiabatic_indx * mach_11 - gamma_minus
    temperature_r *= gamma_minus * mach_11 + 2
    temperature_r /= gamma_plus**2 * mach_11

    # Density ratio
    density_r = gamma_plus * mach_11 / (gamma_minus * mach_11 + 2)

    # Stagnation pressure ratio
    pressure_s1 = gamma_plus / (2 * adiabatic_indx * mach_11 - gamma_minus)
    pressure_s1 **= 1 / gamma_minus
    pressure_s2 = gamma_plus * mach_11 / (gamma_minus * mach_11 + 2)
    pressure_s2 **= adiabatic_indx / gamma_minus
    pressure_s = pressure_s1 * pressure_s2

    normal_shock_dict = {
        "mach_2": mach_2,
        "pressure_r": pressure_r,
        "temperature_r": temperature_r,
        "density_r": density_r,
        "pressure_s": pressure_s,
    }
    return normal_shock_dict  # [ ]


def oblique_shock_relations(
    mach_1: float, shock_angle_deg: float, adiabatic_indx: float = 1.4
) -> dict[str, float]:
    """
    Calculates oblique shock relations for weak shocks

    Parameters:
        mach_1: pre-shock mach number
        shock_angle_deg: shock angle in degrees
        adiabatic_indx: adiabatic index, 1.4 (default)

    Returns:
        dict: A dictionary containing:
            - mach_2: post-shock mach number
            - pressure_r: pressure ratio (post-shock / pre-shock)
            - temperature_r: temperature ratio (post-shock / pre-shock)
            - density_r: density ratio (post-shock / pre-shock)
            - deflection_angle_degs: deflection angle in [degs]
            - mach_n1: normal pre-shock mach number
            - mach_n2: normal post-shock mach number
    Examples:
        >> oblique_shock_relations(3.0, 45.0)

    Reference:
        Modern Compressible Flow With Historic Perspective, International
        Edition 4th (Anderson J., ISBN 978 1 260 57082 3)
    """
    shock_angle = np.radians(shock_angle_deg)
    gamma_minus = adiabatic_indx - 1
    gamma_plus = adiabatic_indx + 1

    # Normal pre-shock mach number (Eq. 4.7)
    mach_n1 = mach_1 * np.sin(shock_angle)
    mach_n11 = mach_n1**2

    # Check if the normal mach is greater than 1
    if not mach_n1 > 1:
        return print("No shock, found")

    # Deflection angle (Eq. 4.17)
    tan_deflection_ang = 2 / np.tan(shock_angle)
    tan_deflection_ang *= mach_n11 - 1
    tan_deflection_ang /= mach_1**2 * (adiabatic_indx + np.cos(2 * shock_angle)) + 2
    deflection_angle_deg = np.degrees(np.arctan(tan_deflection_ang))

    # Density ratio (Eq. 4.8)
    density_r = gamma_plus * mach_n11
    density_r /= gamma_minus * mach_n11 + 2

    # Pressure ratio (Eq. 4.9)
    pressure_r = 2 * adiabatic_indx * (mach_n11 - 1) / gamma_plus + 1

    # Normal post-shock mach number (Eq. 4.10)
    mach_n22 = mach_n11 + 2 / gamma_minus
    mach_n22 /= 2 * adiabatic_indx * mach_n11 / gamma_minus - 1

    # Temperature ratio (Eq. 4.11)
    temperature_r = pressure_r * 1 / density_r

    # Post-shock mach number (Eq. 4.12)
    mach_2 = mach_n22**0.5
    mach_2 /= np.sin(np.radians(shock_angle_deg - deflection_angle_deg))

    oblique_shock_dict = {
        "mach_2": mach_2,
        "pressure_r": pressure_r,
        "temperature_r": temperature_r,
        "density_r": density_r,
        "deflection_angle_degs": deflection_angle_deg,
        "mach_n2": mach_n22**0.5,
        "mach_n1": mach_n1,
    }
    return oblique_shock_dict

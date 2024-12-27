"""
    Date:   03/26/2023
    Author: Martin E. Liza
    File:   optics.py
    Def:    Contains aero optics functions.
"""

from ambiance import Atmosphere
import molmass
import numpy as np
import scipy.constants as s_consts
from haot import constants as constants_tables
from haot import quantum_mechanics as quantum
from haot import conversions


def index_of_refraction(mass_density_dict: dict[str, float]) -> dict[str, float]:
    """
    Calculates dilute and dense index of refraction

    Parameters:
        mass density dictionary in [kg/m^3]

    Returns:
        dict: A dictionary containing
            - dilute: dilute index of refraction
            - dense: dense index of refraction
    """
    pol_consts = constants_tables.polarizability()  # [cm3]
    molar_density = {
        key: conversions.mass_density_to_molar_density(value, key)
        for key, value in mass_density_dict.items()
    }
    # a_i * N_i
    n_const = {
        key: conversions.polarizability_cgs_to_si(pol_consts[key]) * molar_density[key]
        for key in mass_density_dict.keys()
    }
    # Sum (a_i N_i)
    tot_pol_molar = sum(n_const.values())

    # Calculates dilute and dense index of refraction
    n_return = {}
    n_return["dilute"] = 1 + tot_pol_molar / (2 * s_consts.epsilon_0)
    n_temp = tot_pol_molar / (3 * s_consts.epsilon_0)
    n_return["dense"] = ((2 * n_temp + 1) / (1 - n_temp)) ** 0.5

    return n_return


def dielectric_material_const(n_dict: dict[str, float]) -> dict[str, float]:
    """
    Calculates the dielectric medium's constant

    Parameters:
        n_dict: dilute and dense formulation

    Returns:
        dict: A dictionary containing
            - dilute: dilute dielectric constant
            - dense: dense dielectric constant
    """
    # n ~ sqrt(e_r)
    dielectric = {key: s_consts.epsilon_0 * n_dict[key] ** 2 for key in n_dict.keys()}
    return dielectric


def optical_path_length(n_solution, distance):
    OPL = {}
    OPL["dilute"] = n_solution["dilute"] * distance
    OPL["dense"] = n_solution["dense"] * distance
    # TODO: Missing implementation
    print("TODO: Missing this implementation")


def tropina_aproximation(vibrational_number, rotational_number, molecule):
    electron_mass = s_consts.m_e
    electron_charge = s_consts.e
    spectroscopy_const = constants_tables.spectroscopy_constants(molecule)
    # resonance_distance = omega_gi - omega
    # TODO: Missing implementation
    print("TODO: Missing this implementation")


def buldakov_expansion(
    vibrational_number: int, rotational_number: int, molecule: str
) -> float:
    """
    Calculates the Buldakov expansion

    Parameters:
        vibrational_number: vibrational quantum number (has to be positive)
        rotational_number: rotational quantum number (has to be positive)
        molecule: H2, N2, O2

    Returns:
        buldakov expansion in [m^3]

    Reference:
        Temperature Dependence of Polarizability of Diatomic Homonuclear
        Molecules (https://doi.org/10.1134/BF03355985)
    """
    # Load constants
    spectroscopy_const = constants_tables.spectroscopy_constants(molecule)
    derivative_const = constants_tables.buldakov_polarizability_derivatives_2016(
        molecule
    )
    be_we = spectroscopy_const["B_e"] / spectroscopy_const["omega_e"]

    # Dunham potential energy constants
    (a_0, a_1, a_2) = quantum.potential_dunham_coef_012(molecule)
    a_3 = quantum.potential_dunham_coeff_m(a_1, a_2, 3)

    rotational_degeneracy = rotational_number * (rotational_number + 1)
    vibrational_degeneracy = 2 * vibrational_number + 1

    # Split in terms
    tmp_1 = be_we
    tmp_1 *= -3 * a_1 * derivative_const["first"] + derivative_const["second"]
    tmp_1 *= vibrational_degeneracy
    tmp_1 *= 1 / 2

    tmp_2 = be_we**2
    tmp_2 *= derivative_const["first"]
    tmp_2 *= rotational_degeneracy
    tmp_2 *= 4

    tmp_31a = 7
    tmp_31a += 15 * vibrational_degeneracy**2
    tmp_31a *= a_1**3
    tmp_31a *= -3 / 8

    tmp_31b = 23
    tmp_31b += 39 * vibrational_degeneracy**2
    tmp_31b *= a_2
    tmp_31b *= a_1
    tmp_31b *= 1 / 4

    tmp_31c = 5
    tmp_31c += vibrational_degeneracy**2
    tmp_31c *= a_3
    tmp_31c *= -15 / 4

    tmp_31 = derivative_const["first"] * (tmp_31a + tmp_31b + tmp_31c)

    tmp_32a = 7
    tmp_32a += 15 * vibrational_degeneracy**2
    tmp_32a *= a_1**2
    tmp_32a *= 1 / 8

    tmp_32b = 5
    tmp_32b += vibrational_degeneracy**2
    tmp_32b *= a_2
    tmp_32b * --3 / 4

    tmp_32 = derivative_const["second"] * (tmp_32a + tmp_32b)

    tmp_33 = 7
    tmp_33 += 15 * vibrational_degeneracy**2
    tmp_33 *= a_1
    tmp_33 *= derivative_const["third"]
    tmp_33 *= -1 / 24

    tmp_3 = (tmp_31 + tmp_32 + tmp_33) * be_we**2

    tmp_41 = 1 - a_2
    tmp_41 *= 24
    tmp_41 += 27 * a_1 * (1 + a_1)
    tmp_41 *= derivative_const["first"]

    tmp_42 = 1 + 3 * a_1
    tmp_42 *= derivative_const["second"]
    tmp_42 *= -3

    tmp_43 = 1 / 8 * derivative_const["third"]

    tmp_4 = tmp_41 + tmp_42 + tmp_43
    tmp_4 *= rotational_degeneracy
    tmp_4 *= vibrational_degeneracy
    tmp_4 *= be_we**3

    return derivative_const["zeroth"] + tmp_1 + tmp_2 + tmp_3 + tmp_4


def kerl_polarizability_temperature(
    temperature_K: float, molecule: str, wavelength_nm: float
) -> float:
    """
    Calculates the polarizability using Kerl's extrapolation

    Parameters:
        temperature_K: reference temperature in [K]
        molecule: H2, N2, O2, Air
        wavelength_nm: signal's wavelength in [nm]

    Returns:
        polarizability in [m^3]

    Reference:
        Polarizability a(w,T,rho) of Small Molecules in the Gas Phase
        (https://doi.org/10.1002/bbpc.19920960517)

    Examples:
        >> kerl_polarizability_temperature(600.0, 'N2', 533.0)
    """
    # Checking cases
    if type(temperature_K) is float and temperature_K < 0:
        raise ValueError("Temperature must be greater than 0 Kelvin!")
    if type(temperature_K) is np.ndarray and (temperature_K < 0).any():
        raise ValueError("Temperature must be greater than 0 Kelvin!")
    if wavelength_nm <= 0:
        raise ValueError("Wavelength must be greater than 0 nanometers!")
    if molecule not in ["Air", "H2", "N2", "O2"]:
        raise ValueError("This function only supports Air, H2, N2 or O2")
    # Check sizes
    mean_const = constants_tables.kerl_interpolation(molecule)
    angular_frequency = 2 * np.pi * s_consts.speed_of_light / (wavelength_nm * 1e-9)

    tmp = mean_const["c"] * temperature_K**2
    tmp += mean_const["b"] * temperature_K
    tmp += 1
    tmp *= mean_const["groundPolarizability"]
    tmp /= 1 - (angular_frequency / mean_const["groundFrequency"]) ** 2

    return tmp  # [m^3]


def atmospheric_index_of_refraction(
    altitude_m: float, vapor_pressure: float = 0.0
) -> float:
    """
    Calculates the atmospheric index of refraction as a function of altitude

    Parameters:
        altitude_m: altitude in [m]
        vapor_pressure: vapor pressure at given altitude in [mbar], 0.0 (default)
        temperature_K: reference temperature in [K]

    Returns:
        index of refraction in [ ]

    Reference:
        The constants in the equation for atmospheric refractive index at radio frequencies (https://ieeexplore.ieee.org/document/4051437)
    """
    atmospheric_prop = Atmosphere(altitude_m)
    temperature = atmospheric_prop.temperature  # [K]
    pressure = atmospheric_prop.pressure * 0.01  # [mbar]
    [K_1, K_2] = constants_tables.smith_atmospheric_constants()

    refractivity = K_2 * vapor_pressure / temperature
    refractivity += pressure
    refractivity *= K_1 / temperature
    refractivity *= 10**-6

    return refractivity + 1


def gladstone_dale_constant(
    mass_density_dict: dict[str, float] = None
) -> dict[str, float]:
    """
    Calculates Gladstone-Dale constants, returns constants haot.constants if
    mass_density_dict is not provided

    Parameters:
        mass density dictionary in [kg/m^3], None (default)

    Returns:
        dict: A dictionary containing
            - Species Gladstone-Dale constants in [m3/kg]
    """
    pol_consts = constants_tables.polarizability()  # [cm^3]

    # Convert polarizability CGS to SI
    pol_SI = {
        key: conversions.polarizability_cgs_to_si(pol_consts[key])
        for key in pol_consts.keys()
    }
    species_GD = {}

    # Calculates species GD
    for key, val in pol_SI.items():
        species_GD[key] = val * s_consts.N_A / molmass.Formula(key).mass
        species_GD[key] /= 2 * s_consts.epsilon_0
        species_GD[key] *= 1e3  # converts [1/g] to [1/kg]

    # Calculate total GD
    if not mass_density_dict:
        return species_GD
    else:
        tot_density = sum(mass_density_dict.values())
        for key in mass_density_dict.keys():
            species_GD[key] *= mass_density_dict[key] / tot_density

        species_GD["gladstone_dale"] = sum(species_GD.values())
        return species_GD  # [m3/kg]

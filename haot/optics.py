"""
    Date:   03/26/2023
    Author: Martin E. Liza
    File:   aero_optics.py
    Def:    Contains aero optics functions.
"""

from ambiance import Atmosphere
import numpy as np
import scipy.constants as s_consts
from haot import aerodynamics as aero
from haot import constants as constants_tables
from haot import quantum_mechanics as quantum


def polarizability_cgs_to_si(polarizability_cgs: float) -> float:
    """
    Converts volumetric polarizability (CGS) to atomic polarizability (SI)

    Parameters:
        polarizability_cgs: volumetric polarizability in [cm^3]

    Returns
        atomic polarizability in [Fm^2]
    """
    return polarizability_cgs * 4 * np.pi * s_consts.epsilon_0 * 1e-6


def polarizability_si_to_cgs(polarizability_si: float) -> float:
    """
    Converts atomic polarizability (SI) to volumetric polarizability (CGS)

    Parameters:
        polarizability_si: atomic polarizability in [Fm^2]

    Returns:
        volumetric polarizability in [cm^3]

    """
    return polarizability_si * 1e6 / (4 * np.pi * s_consts.epsilon_0)


def gas_density(density_dict: dict[str, float]) -> dict[str, float]:
    """
    Calculates gas density in [particles/m^3]

    Parameters:
        species in [kg/m^3]

    Returns:
        species density in [particles/m^3]

    """
    gas_amu_weight = aero.air_atomic_molar_mass()  # [g/mol]
    gas_density = {
        i: (density_dict[i] * 10**3 * s_consts.N_A / gas_amu_weight[i])
        for i in density_dict
    }

    return gas_density  # [particles/m^3]


def index_of_refraction(gas_density_dict: dict[str, float]) -> dict[str, float]:
    """
    Calculates index of refraction

    Parameters:
        species in [kg/m^3]

    Returns:
        dict: A dictionary containing
            - dilute: dilute index of refraction
            - dense: dense index of refraction
    """
    pol_consts = constants_tables.polarizability()  # [m^3]
    density = gas_density(gas_density_dict)  # [particles/m3]
    n_const = {}  # [ ]

    for i in gas_density_dict:
        # Pol constants are in m^3
        alpha = polarizability_cgs_to_si(pol_consts[i] * 1e6)
        n_const[i] = alpha * density[i]  # (a_i N_i)

    # add all n_i
    temp = 0.0
    for i in n_const.values():
        temp += i

    n_return = {}
    n_return["dilute"] = 1 + temp / (2 * s_consts.epsilon_0)
    n_temp = temp / (3 * s_consts.epsilon_0)
    n_return["dense"] = ((2 * n_temp + 1) / (1 - n_temp)) ** 0.5

    return n_return


def dielectric_material_const(n_const):
    # n ~ sqrt(e_r)
    dielectric_const_0 = s_consts.epsilon_0  # [F/m]
    dielectric = {}
    dielectric["dilute"] = dielectric_const_0 * n_const["dilute"] ** 2
    dielectric["dense"] = dielectric_const_0 * n_const["dense"] ** 2
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


# Calculate polarizability (uses equation 4 from the paper)
def buldakov_expansion(vibrational_number, rotational_number, molecule):
    # Load constants
    spectroscopy_const = constants_tables.spectroscopy_constants(molecule)
    derivative_const = constants_tables.polarizability_derivatives(molecule)
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
    K_1 = 79  # [K/mbar]
    K_2 = 4800  # [K]

    refractivity = K_2 * vapor_pressure / temperature
    refractivity += pressure
    refractivity *= K_1 / temperature
    refractivity *= 10**-6

    return refractivity + 1


def atmospheric_gladstone_dale(altitude=0.0, gas_composition_dict=None):
    gladstone_const = gladstone_dale()  # [m3/kg]

    if gas_composition_dict is None:
        gas_composition_dict = {}
        gas_composition_dict["N"] = 0.0
        gas_composition_dict["O"] = 0.0
        gas_composition_dict["NO"] = 0.0
        gas_composition_dict["N2"] = 0.79
        gas_composition_dict["O2"] = 0.21

    tmp = 0
    for i in gas_composition_dict.keys():
        tmp += gas_composition_dict[i] * gladstone_const[i]

    return tmp


def gladstone_dale(gas_density_dict=None):  # [kg/m3
    gas_amu_weight = aero.air_atomic_molar_mass()  # [g/mol]
    avogadro_number = s_consts.N_A  # [particles/mol]
    dielectric_const = s_consts.epsilon_0  # [F/m]
    pol_consts = constants_tables.polarizability()  # [m^3]

    # Convert CGS to SI
    pol_consts.update(
        {n: 4 * np.pi * dielectric_const * pol_consts[n] for n in pol_consts.keys()}
    )  # [Fm^2]

    # Calculate Gladstone dale
    gladstone_dale_const = {}
    for i in pol_consts:
        gladstone_dale_const[i] = (
            pol_consts[i]
            / (2 * dielectric_const)
            * (avogadro_number / gas_amu_weight[i])
            * 1e3
        )  # [m^3/kg]

    gladstone_dale_dict = {}
    if not gas_density_dict:
        return gladstone_dale_const  # [m^3/kg]
    else:
        gladstone_dale_dict["gladstone_dale"] = 0.0
        for i in gas_density_dict:
            gladstone_dale_dict[i] = (
                gladstone_dale_const[i] * gas_density_dict[i]
            ) / sum(gas_density_dict.values())
            gladstone_dale_dict["gladstone_dale"] += (
                gladstone_dale_const[i] * gas_density_dict[i]
            )
        gladstone_dale_dict["gladstone_dale"] /= sum(gas_density_dict.values())

        return gladstone_dale_dict  # [m^3/kg]

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


def index_of_refraction_density_temperature(
    temperature_K: float,
    mass_density: float,
    molecule: str = "Air",
    wavelength_nm: float = 633,
) -> dict[str, float]:
    """
    Calculates dilute and dense index of refraction as a
    function of mass density and temperature.
    Uses Kerl approximation for polarizability

    Parameters:
        temperature_K: reference temperature in [K]
        mass density in [kg/m^3]
        molecule: H2, N2, O2, Air(default)
        wavelength_nm: signal's wavelength in [nm], 633(default) [nm]

    Returns:
        dict: A dictionary containing
            - dilute: dilute index of refraction
            - dense: dense index of refraction
    """
    # Calculates polarizability using Kerl
    pol_kerl_air_m3 = kerl_polarizability_temperature(
        temperature_K, "Air", wavelength_nm
    )
    pol_kerl_SI = conversions.polarizability_cgs_to_si(pol_kerl_air_m3 * 1e6)

    if molecule == "Air":
        molar_mass_air = (
            0.78 * molmass.Formula("O2").mass
            + 0.21 * molmass.Formula("N2").mass
            + 0.01 * molmass.Formula("Ar").mass
        )
        molar_density = mass_density * s_consts.N_A / molar_mass_air * 1e3
    else:
        molar_density = conversions.mass_density_to_molar_density(
            mass_density, molecule
        )
    tot_pol_molar = molar_density * pol_kerl_SI
    n_return = {}
    n_return["dilute"] = 1 + tot_pol_molar / (2 * s_consts.epsilon_0)
    n_temp = tot_pol_molar / (3 * s_consts.epsilon_0)
    n_return["dense"] = ((2 * n_temp + 1) / (1 - n_temp)) ** 0.5

    return n_return


def index_of_refraction(mass_density_dict: dict[str, float]) -> dict[str, float]:
    """
    Calculates dilute and dense index of refraction as a
    function of mass  density

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


def dielectric_material_const(index_of_refraction: float) -> float:
    """
    Calculates the dielectric medium's constant

    Parameters:
        index_of_refraction: index of refraction

    Returns:
        material's dielectric constant in [F/m]
    """
    # n ~ sqrt(e_r)
    return s_consts.epsilon_0 * index_of_refraction**2


def optical_path_length(index_of_refraction: float, distance: float) -> float:
    """
    Calculates the optical path length

    Parameters:
        index_of_refraction: index of refraction
        distance: length

    Returns:
        Optical Path Length in units of distance
    """
    if np.shape(index_of_refraction) != np.shape(distance):
        raise ValueError("Index of refraction and distance must have the same length")
    index_avg = 0.5 * (index_of_refraction[:-1] + index_of_refraction[1:])
    return np.mean(np.cumsum(index_avg * np.diff(distance)))


def optical_path_difference_rms(opd: float, avg_ax: int = 0) -> float:
    """
    Calculates the optical path difference RMS

    Parameters:
        opd: Optical Path Difference
        avg_ax: axis where average is performed, 0 (default)

    Returns
        Optical Path Difference Root-Mean-Squared
    """
    # Validate the input array
    if not isinstance(opd, np.ndarray):
        raise ValueError("opd must be a numpy array")

    if avg_ax not in [0, 1, 2, 3]:
        raise ValueError("avg_ax must be one of [0, 1, 2, 3]")
    return np.sqrt(np.mean((opd - np.mean(opd, axis=avg_ax, keepdims=True)) ** 2))


def phase_variance(opd_rms: float, wavelength_nm: float) -> float:
    """
    Calculates phase variance

    Parameters:
        opd_rm: Optical Path Difference RMS in units of [m]
        wavelength_nm: Wavelength of light in units of [nm]

    Returns:
        Phase difference, unit-less
    """
    return (2 * np.pi * opd_rms / (wavelength_nm * 1e-9)) ** 2


def strehl_ratio(phase_variance: float) -> float:
    """
    Calculates the Strehl ratio

    Parameters:
        phase_variance: phase variance

    Returns:
        Strehl ratio
    """

    return np.exp(-phase_variance)


def optical_path_difference(opl: np.array, avg_ax: int = 0) -> float:
    """
    Calculates the optical path difference

    Parameters:
        opl: has to be a numpy array of shape [time, x_axis, y_axis, z_axis]
        avg_ax: axis where average is performed, 0 (default)

    Returns:
        numpy array of the same shape as the optical path length
    """
    if not isinstance(opl, np.ndarray):
        raise ValueError("opl must be a numpy array")

    if avg_ax not in [0, 1, 2, 3]:
        raise ValueError("avg_ax must be one of [0, 1, 2, 3]")
    return opl - np.mean(opl, axis=avg_ax, keepdims=True)


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


def brewster_angle(
    medium_index_of_refraction: float, vacuum_index_of_refraction: float = 1.0
) -> float:
    """
    Calculates the Brewster angle. Note,
    that medium_index_of_refraction should be greater than the
    vacuum_index_of_refraction

    Parameters:
        medium_index_of_refraction: medium's index of refraction
        vacuum_index_of_refraction: vacuum's index of refraction, 1.0 (default)

    Returns:
        Brewster angle in [degs]
    """
    return np.rad2deg(
        np.arctan(medium_index_of_refraction / vacuum_index_of_refraction)
    )


def total_internal_reflection_angle(
    medium_index_of_refraction: float, vacuum_index_of_refraction: float = 1.0
) -> float:
    """
    Calculates the critical angle that causes total internal reflection. Note,
    that medium_index_of_refraction should be greater than the
    vacuum_index_of_refraction

    Parameters:
        medium_index_of_refraction: medium's index of refraction
        vacuum_index_of_refraction: vacuum's index of refraction, 1.0 (default)

    Returns:
        Critical angle in [degs]
    """
    return np.rad2deg(
        np.arcsin(vacuum_index_of_refraction / medium_index_of_refraction)
    )


def normal_incidence_reflectance(
    medium_index_of_refraction: float, vacuum_index_of_refraction: float = 1.0
) -> float:
    """
    Calculates the reflectance at a normal incidence. Note,
    that medium_index_of_refraction should be greater than the
    vacuum_index_of_refraction

    Parameters:
        medium_index_of_refraction: medium's index of refraction
        vacuum_index_of_refraction: vacuum's index of refraction, 1.0 (default)

    Returns:
       Reflectance at normal incidence in [ ]
    """
    return (
        (medium_index_of_refraction - vacuum_index_of_refraction)
        / (medium_index_of_refraction + vacuum_index_of_refraction)
    ) ** 2


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

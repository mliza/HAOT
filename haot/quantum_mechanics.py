"""
    Date:   11/26/2024
    Author: Martin E. Liza
    File:   quantum_mechanics.py
    Def:    Contains Quantum Mechanics functions.
"""

import molmass
import numpy as np
import scipy.constants as s_consts
from haot import constants as constants_tables
from haot import conversions


def zero_point_energy(molecule: str) -> float:
    """
    Calculates zero-point energy (ZPE) using spectroscopy constants for
    diatomic molecules

    Parameters:
        molecule: NO+, N2+, O2+, NO, N2, O2, H2

    Returns:
        zero point energy in [cm^-1]

    Examples:
        >> zero_point_energy('O2')

    Reference:
        Experimental Vibrational Zero-Point Energies Diatomic Molecules
        (https://doi.org/10.1063/1.2436891)
    """
    if molecule not in ["NO+", "N2+", "O2+", "NO", "N2", "O2", "H2"]:
        raise ValueError("This function only supports NO+, N2+, O2+, NO, N2, O2, H2")
    spectroscopy_const = constants_tables.spectroscopy_constants(molecule)

    scope_var = spectroscopy_const["alpha_e"]
    scope_var *= spectroscopy_const["omega_e"]
    scope_var /= spectroscopy_const["B_e"]

    zpe = spectroscopy_const["omega_e"] / 2
    zpe -= spectroscopy_const["omega_xe"] / 2
    zpe += spectroscopy_const["omega_ye"] / 8
    zpe += spectroscopy_const["B_e"] / 4
    zpe += scope_var / 12
    zpe += scope_var**2 / (144 * spectroscopy_const["B_e"])
    return zpe  # [1/cm]


def vibrational_partition_function(
    vibrational_number: int, temperature_K: float, molecule: str
) -> float:
    """
    Calculates the vibrational partition function base in the harmonic
    terms only for diatomic molecules

    Parameters:
        vibrational_number: vibrational quantum number (has to be positive)
        temperature_K: reference temperature in [K]
        molecule: NO+, N2+, O2+, NO, N2, O2

    Returns:
        vibrational partition function

    Examples:
        >> vibrational_partition_function(2, 350.0, 'N2')
    """
    z_vib = 0.0
    for v in range(vibrational_number + 1):
        z_vib += boltzmann_factor(
            temperature_K,
            molecule,
            vibrational_number=v,
            rotational_number=None,
            born_opp_flag=False,
        )
    return z_vib


def rotational_partition_function(
    rotational_number: int, temperature_K: float, molecule: str
) -> float:
    """
    Calculates the rotational partition function base in the harmonic
    terms only for diatomic molecules

    Parameters:
        rotational_number: rotational quantum number (has to be positive)
        temperature_K: reference temperature in [K]
        molecule: NO+, N2+, O2+, NO, N2, O2

    Returns:
        rotational partition function

    Examples:
        >> rotational_partition_function(2, 350.0, 'N2')
    """
    z_rot = 0.0
    for j in range(rotational_number + 1):
        z_rot += boltzmann_factor(
            temperature_K,
            molecule,
            vibrational_number=None,
            rotational_number=j,
            born_opp_flag=False,
        )
    return z_rot


def born_oppenheimer_partition_function(
    vibrational_number: int, rotational_number: int, temperature_K: float, molecule: str
) -> float:
    """
    Calculates the partition function using the Born-Oppenheimer approximation

    Parameters:
        vibrational_number: vibrational quantum number (has to be positive)
        rotational_number: rotational quantum number (has to be positive)
        temperature_K: reference temperature in [K]
        molecule: NO+, N2+, O2+, NO, N2, O2

    Returns:
        partition function using the Born-Oppenheimer approximations

    Examples:
        >> born_oppenheimer_partition_function(2, 4, 500.0, 'O2')
    """
    z_bo = 0.0
    for j in range(rotational_number + 1):
        for v in range(vibrational_number + 1):
            z_bo += boltzmann_factor(
                temperature_K,
                molecule,
                vibrational_number=v,
                rotational_number=j,
                born_opp_flag=True,
            )
    return z_bo


def potential_dunham_coef_012(molecule: str) -> tuple[float, float, float]:
    """
    Calculates the 0th, 1st, and 2nd Dunham potential coefficients for diatomic
    molecules

    Parameters:
        molecule: NO+, N2+, O2+, NO, N2, O2

    Returns:
       Dunham potential coefficients 0, 1, and 2

    Examples:
        >> [a_0, a_1, a_2] = potential_dunham_coef_012('N2')

    Reference:
        Dunham potential energy coefficients of the hydrogen halides and carbon
        monoxide (https://doi.org/10.1016/0022-2852(76)90323-4)

        Anharmonic Potential Constants and Their Dependence upon Bond Length
        (https://doi.org/10.1063/1.1731952)
    """
    spectroscopy_const = constants_tables.spectroscopy_constants(molecule)
    a_0 = spectroscopy_const["omega_e"] ** 2 / (4 * spectroscopy_const["B_e"])
    a_1 = -(
        spectroscopy_const["alpha_e"]
        * spectroscopy_const["omega_e"]
        / (6 * spectroscopy_const["B_e"] ** 2)
        + 1
    )
    a_2 = (5 / 4) * a_1**2 - (2 / 3) * (
        spectroscopy_const["omega_xe"] / spectroscopy_const["B_e"]
    )
    return (a_0, a_1, a_2)


def potential_dunham_coeff_m(a_1: float, a_2: float, m: int) -> float:
    """
    Calculates higher order Dunham potential coefficients

    Parameters:
        a_1: Dunham potential coefficient
        a_2: Dunham potential coefficient
        m:   desired Dunham potential coefficient

    Returns:
       Dunham potential coefficient at m

    Examples:
        >> a_3 = potential_dunham_coef_m(a_1, a_2, 3)

    Reference:
        A recursion formula for the coefficients of Dunham potential (https://doi.org/10.1016/j.theochem.2003.12.003)
    """

    tmp = (12 / a_1) ** (m - 2)
    tmp *= 2 ** (m + 1) - 1
    tmp *= (a_2 / 7) ** (m - 1)
    for i in range(m - 2):
        tmp *= 1 / (m + 2 - i)

    return tmp


def boltzmann_factor(
    temperature_K: float,
    molecule: str,
    vibrational_number: int = None,
    rotational_number: int = None,
    born_opp_flag: bool = False,
) -> float:
    """
    Calculates the Boltzmann factor at a given vibrational quantum number,
    and/or rotational quantum number. If the born_opp_flag is provided, it will
    calculate the Boltzmann factor using the Born-Oppenheimer approximation

    Parameters:
        temperature_K: reference temperature in [K]
        molecule: NO+, N2+, O2+, NO, N2, O2
        vibrational_number: vibrational quantum number (has to be positive)
        rotational_number: rotational quantum number (has to be positive)
        born_opp_flag: uses the Born-Oppenheimer approximation, False (default)

    Returns:
        Boltzmann factor at given temperature, rotational and/or vibrational quantum number

    Examples:
        >> boltzmann_factor(500.0, 'O2', 3)

        >> boltzmann_factor(500.0, 'O2', 3, 2)

        >> boltzmann_factor(500.0, 'O2', 3, 1, True)
    """
    # Initialize energy terms, degeneracy and thermal beta
    energy_vib_k = 0
    energy_rot_k = 0
    degeneracy_rotation = 1
    thermal_beta = 1 / (s_consts.k * temperature_K)

    # Calculates Energy levels
    if not born_opp_flag:
        if vibrational_number is not None:
            energy_vib_k = vibrational_energy_k(vibrational_number, molecule)
        if rotational_number is not None:
            energy_rot_k = rotational_energy_k(rotational_number, molecule)
            degeneracy_rotation = 2 * rotational_number + 1
        tot_energy = conversions.wavenumber_to_joules(energy_vib_k + energy_rot_k)
    else:
        degeneracy_rotation = 2 * rotational_number + 1
        tot_energy = conversions.wavenumber_to_joules(
            born_oppenheimer_approximation(
                vibrational_number, rotational_number, molecule
            )
        )
    return degeneracy_rotation * np.exp(-tot_energy * thermal_beta)


def boltzmann_distribution(
    temperature_K: float,
    molecule: str,
    vibrational_number: int = None,
    rotational_number: int = None,
    born_opp_flag: bool = False,
) -> np.ndarray:
    """
    Calculates the Boltzmann distribution function at a given vibrational quantum number,
    and/or rotational quantum number. If the born_opp_flag is provided, it will
    calculate the Boltzmann factor using the Born-Oppenheimer approximation

    Parameters:
        temperature_K: reference temperature in [K]
        molecule: NO+, N2+, O2+, NO, N2, O2
        vibrational_number: vibrational quantum number (has to be positive)
        rotational_number: rotational quantum number (has to be positive)
        born_opp_flag: uses the Born-Oppenheimer approximation, False (default)

    Returns:
        Boltzmann distribution [vibrational_number, rotational_number]

    Examples:
        >> boltzmann_distribution(500.0, 'O2', 3, 2, true)
    """
    # Calculates partition functions if vibrational or rotational numbers are provided
    if not born_opp_flag:
        z_rot = 1
        z_vib = 1
        if vibrational_number is not None:
            z_vib = vibrational_partition_function(
                vibrational_number, temperature_K, molecule
            )
        if rotational_number is not None:
            z_rot = rotational_partition_function(
                rotational_number, temperature_K, molecule
            )
        z_tot = z_rot * z_vib
    else:
        z_tot = born_oppenheimer_partition_function(
            vibrational_number, rotational_number, temperature_K, molecule
        )

    # Create the distribution array based on the inputs provided
    if vibrational_number and rotational_number:
        tmp = np.zeros([vibrational_number + 1, rotational_number + 1])
        for j in range(rotational_number + 1):
            for v in range(vibrational_number + 1):
                tmp[v][j] = boltzmann_factor(
                    temperature_K=temperature_K,
                    molecule=molecule,
                    vibrational_number=v,
                    rotational_number=j,
                    born_opp_flag=born_opp_flag,
                )
    elif vibrational_number:
        tmp = np.zeros(vibrational_number + 1)
        for v in range(vibrational_number + 1):
            tmp[v] = boltzmann_factor(
                temperature_K=temperature_K,
                molecule=molecule,
                vibrational_number=v,
                rotational_number=None,
                born_opp_flag=False,
            )
    elif rotational_number:
        tmp = np.zeros(rotational_number + 1)
        for j in range(rotational_number + 1):
            tmp[j] = boltzmann_factor(
                temperature_K=temperature_K,
                molecule=molecule,
                vibrational_number=None,
                rotational_number=j,
                born_opp_flag=False,
            )
    return tmp / z_tot


def born_oppenheimer_approximation(
    vibrational_number: int, rotational_number: int, molecule: str
) -> float:
    """
    Calculates the energy at a rotational and vibrational quantum number,
    using the Born-Oppenheimer approximation

    Parameters:
        molecule: NO+, N2+, O2+, NO, N2, O2
        vibrational_number: vibrational quantum number (has to be positive),
        rotational_number: rotational quantum number (has to be positive)

    Returns:
        Vibro-rotational energy at given vibrational and rotational quantum
        numbers in [cm^-1]

    Examples:
        >> born_oppenheimer_approximation(2,3,'N2')
    """
    spectroscopy_constants = constants_tables.spectroscopy_constants(molecule)

    vib_levels = vibrational_number + 1 / 2
    rot_levels = rotational_number * (rotational_number + 1)

    # Harmonic vibration and rotation terms
    harmonic = spectroscopy_constants["omega_e"] * vib_levels
    harmonic += spectroscopy_constants["B_e"] * rot_levels

    # Anharmonic vibration and rotation terms
    anharmonic = spectroscopy_constants["omega_xe"] * vib_levels**2
    anharmonic += spectroscopy_constants["D_e"] * rot_levels**2

    # Interaction between vibration and rotation modes
    interaction = spectroscopy_constants["alpha_e"] * vib_levels * rot_levels

    return harmonic - anharmonic - interaction  # [cm^1]


def vibrational_energy_k(vibrational_number: int, molecule: str) -> float:
    """
    Calculates the vibrational energy at a given vibrational quantum number,
    using for the harmonic terms

    Parameters:
        vibrational_number: vibrational quantum number (has to be positive)
        molecule: NO+, N2+, O2+, NO, N2, O2

    Returns:
        harmonic terms of the vibrational energy in [cm^-1]

    Examples:
        >> vibrational_energy_k(2, 'N2')
    """
    spectroscopy_constants = constants_tables.spectroscopy_constants(molecule)
    # Calculates the vibrational energy in units of wave number
    vib_levels = vibrational_number + 1 / 2
    return spectroscopy_constants["omega_e"] * vib_levels  # [cm^-1]


def rotational_energy_k(rotational_number: int, molecule: str) -> float:
    """
    Calculates the rotational energy at a given rotational quantum number,
    using for the harmonic terms

    Parameters:
        rotational_number: rotational quantum number (has to be positive)
        molecule: NO+, N2+, O2+, NO, N2, O2

    Returns:
        harmonic terms of the rotational energy in [cm^-1]

    Examples:
        >> rotational_energy_k(2, 'N2')
    """
    spectroscopy_constants = constants_tables.spectroscopy_constants(molecule)
    # Calculates the rotational energy in units of wave number
    rot_levels = rotational_number * (rotational_number + 1)
    return spectroscopy_constants["B_e"] * rot_levels  # [cm^-1]


def reduced_mass_kg(molecule_1: str, molecule_2: str) -> float:
    """
    Calculates the reduced mass in [kg]

    Parameters:
        molecule_1: name of molecule one
        molecule_2: name of molecule two

    Returns:
        reduced mass in [kg]

    Examples:
        >> reduced_mass_kg('N', 'O')
    """
    m_1 = molmass.Formula(molecule_1).mass
    m_2 = molmass.Formula(molecule_2).mass
    mu = m_1 * m_2 / (m_1 + m_2)

    return conversions.molar_mass_to_kilogram(mu)


def molecular_spring_constant(molecule: str) -> float:
    """
    Calculates the molecular spring constant in [N/m]

    Parameters:
        molecule: NO+, N2+, O2+, NO, N2, O2

    Returns:
        spring force constant in [N/m]

    Examples:
        >> molecular_spring_constant('N2')
    """
    # Split masses
    m_1 = molecule[0]
    m_2 = m_1 if molecule[1] == str(2) else molecule[1]
    spectroscopy_const = constants_tables.spectroscopy_constants(molecule)
    mass_kg = reduced_mass_kg(m_1, m_2)
    freq_rps = conversions.wavenumber_to_angular_frequency(
        spectroscopy_const["omega_e"]
    )

    return mass_kg * freq_rps**2


# TODO: Missing Translational Energy
def tranlational_energy(principal_number_x, principal_number_y, principal_number_z):
    print("TODO: Missing implementation of this function")

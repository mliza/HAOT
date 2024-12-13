"""
    Date:   12/11/2024
    Author: Martin E. Liza
    File:   conversions.py
    Def:    Contains conversion functions.
"""

import numpy as np
import scipy.constants as s_consts


def polarizability_cgs_to_si(polarizability_cgs: float) -> float:
    """
    Converts volumetric polarizability (CGS) to atomic polarizability (SI)

    Parameters:
        polarizability_cgs: volumetric polarizability in [cm^3]

    Returns:
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


def wavenumber_to_electronvolt(wavenumber_cm: float) -> float:
    """
    Converts wavenumber to electron-volt

    Parameters:
        wavenumber_cm: energy in [cm^-1]

    Returns:
        energy in [eV]
    """
    return wavenumber_to_joules(wavenumber_cm) / s_consts.eV


def wavenumber_to_joules(wavenumber_cm: float) -> float:
    """
    Converts wavenumber to Joules

    Parameters:
        wavenumber_cm: energy in [cm^-1]

    Returns:
        energy in [J]
    """
    return wavenumber_cm * s_consts.c * 100 * s_consts.h


def molar_mass_to_kilogram(molar_mass_gmol: float) -> float:
    """
    Converts molar mass [g/mol] to mass [kg]

    Parameters:
        molar_mass_gmol: mass in [g/mol]

    Returns:
        mass in [kg]
    """
    return molar_mass_gmol * 1e-3 / s_consts.N_A
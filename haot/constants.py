"""
    Date:   10/26/2024
    Author: Martin E. Liza
    File:   constants.py
    Def:    This file contains constants used across the package
"""


def smith_atmospheric_constants() -> tuple[float, float]:
    """
    Constants used to calculate the atmospheric index of refraction

    Returns:
        K_1 in [K]
        K_2 in [K/mbar]

    Reference:
        The constants in the equation for atmospheric refractive index at radio frequencies (https://ieeexplore.ieee.org/document/4051437)
    """
    return [79.0, 4800.0]


def sutherland_constants(gas: str) -> dict[str, float]:
    """
    Constants used for Sutherland's law of viscosity and thermal conductivity


    Parameters:
        gas: Air (default), Argon, N2, O2

    Returns:
        dict: A dictionary containing
            - temperature_ref: reference temperature in [K]
            - viscosity_ref: reference viscosity in [Ns/m2]
            - sutherland_visc: Sutherland viscosity constant in [K]
            - conductivity_ref: reference thermal conductivity in [W/mK]
            - sutherland_cond: Sutherland thermal conductivity in [K]

    Reference:
        Viscous Fluid Flow, International Edition, 4th (White F., ISBN 978 1 260
        59786). Table 1.1 and Table 1.2
    """

    dict_out = {}
    if gas == "Air":
        dict_out["temperature_ref"] = 273.0  # [K]
        dict_out["viscosity_ref"] = 1.716e-5  # [kg/ms]
        dict_out["sutherland_visc"] = 111.0  # [K]
        dict_out["conductivity_ref"] = 0.0241  # [W/mK]
        dict_out["sutherland_cond"] = 194.0  # [K]

    if gas == "Argon":
        dict_out["temperature_ref"] = 273.0  # [K]
        dict_out["viscosity_ref"] = 2.125e-5  # [kg/ms]
        dict_out["sutherland_visc"] = 114.0  # [K]
        dict_out["conductivity_ref"] = 0.0163  # [W/mK]
        dict_out["sutherland_cond"] = 170.0  # [K]

    if gas == "N2":
        dict_out["temperature_ref"] = 273.0  # [K]
        dict_out["viscosity_ref"] = 1.663e-5  # [kg/ms]
        dict_out["sutherland_visc"] = 107.0  # [K]
        dict_out["conductivity_ref"] = 0.0242  # [W/mK]
        dict_out["sutherland_cond"] = 150.0  # [K]

    if gas == "O2":
        dict_out["temperature_ref"] = 273.0  # [K]
        dict_out["viscosity_ref"] = 1.919e-5  # [kg/ms]
        dict_out["sutherland_visc"] = 139.0  # [K]
        dict_out["conductivity_ref"] = 0.0244  # [W/mK]
        dict_out["sutherland_cond"] = 240.0  # [K]

    return dict_out


def karl_2003() -> dict[str, float]:
    """
    Gladstone-Dale constants used in Karl's experiments

    Returns:
        dict: A dictionary containing
            - species in units of [m^3/kg]

    Reference:
        High Enthalpy Cylinder Flow in HEG, A Basis for CFD Validation (https://arc.aiaa.org/doi/pdf/10.2514/6.2003-4252)
    """

    dict_out = {
        "N": 0.301e-3,
        "O": 0.182e-3,
        "NO": 0.221e-3,
        "N2": 0.238e-3,
        "O2": 0.190e-3,
    }  # [m^3/kg]
    return dict_out


def polarizability() -> dict[str, float]:
    """
    Polarizabilities at laboratory conditions

    Returns:
        dict: A dictionary of polarizabilities
            - Volumetric species polarizability in CGS units of [cm^3]

    Reference:
        Handbook of Chemistry and Physics, 95th edition
        (https://doi.org/10.1201/b17118)

        Optical spectroscopy of high L n=10 Rydberg states of nitrogen
        (https://doi.org/10.1103/PhysRevA.54.314)

        A numerical study of coupled Hartree-Fock theory for open-shell systems
        (https://doi.org/10.1080/00268977500102811)

        Ab initio calculations of the properties of NO+ in its ground
        electronic state X 1Sigma+
        (https://doi.org/10.1016/0009-2614(93)89356-M)

        Analysis of the 8f, 9f, and 10f, v=1 Rydberg states of N2
        (https://doi.org/10.1103/PhysRevA.44.3007)
    """
    dict_out = {
        "N+": 0.559e-24,
        "O+": 0.345e-24,
        "NO+": 1.021e-24,
        "N2+": 2.386e-24,
        "O2+": 0.238e-24,
        "N": 1.100e-24,
        "O": 0.802e-24,
        "NO": 1.700e-24,
        "N2": 1.7403e-24,
        "O2": 1.5689e-24,
    }  # [cm^3]
    return dict_out


def buldakov_polarizability_derivatives_2016(molecule: str) -> dict[str, float]:
    """
    Polarizability derivative constants

    Parameters:
        molecule: H2, N2, O2

    Reference:
        Temperature dependence of polarizability of diatomic homonuclear
        molecules (https://doi.org/10.1134/BF03355985)

        Handbook of Chemistry and Physics, 95th edition
        (https://doi.org/10.1201/b17118)

        An accurate calculation of the polarizability of the hydrogen molecule
        and its dependence on rotation, vibration and isotopic substitution
        (https://doi.org/10.1080/00268978000103191)

        Theoretical study of the effects of vibrational‐rotational interactions
        on the Raman spectrum of N2 (https://doi.org/10.1063/1.445482)

        Frequency-dependent polarizability of O2 and van der Waals coefficients
        of dimers containing O2 (https://doi.org/10.1063/1.467256)
    """

    dict_out = {}
    if molecule == "H2":
        dict_out["zeroth"] = 0.7849e-30
        dict_out["first"] = 0.90e-30
        dict_out["second"] = 0.49e-30
        dict_out["third"] = -0.85e-30

    if molecule == "N2":
        dict_out["zeroth"] = 1.7801e-30
        dict_out["first"] = 1.86e-30
        dict_out["second"] = 1.2e-30
        dict_out["third"] = -4.6e-30

    if molecule == "O2":
        dict_out["zeroth"] = 1.6180e-30
        dict_out["first"] = 1.76e-30
        dict_out["second"] = 3.4e-30
        dict_out["third"] = -23.7e-30

    return dict_out  # [m^3]


def kerl_interpolation(molecule: str) -> dict[str, float]:
    """
    Constants used in Kerl's extrapolation method

    Parameters:
        molecules: H2, N2, O2, Air

    Returns:
        dict: A dictionary containing
            - groundPolarizability: polarizability at ground level in [m^3]
            - groundFrequency: frequency at ground level in [Hz]
            - b: extrapolation constant in [1/K]
            - c: extrapolation constant in [1/K^2]

    Reference:
        Polarizability a(w,T,rho) of Small Molecules in the Gas Phase
        (https://doi.org/10.1002/bbpc.19920960517)
    """

    dict_out = {}
    if molecule == "H2":
        dict_out["groundPolarizability"] = 0.80320e-30  # [m^3]
        dict_out["groundFrequency"] = 2.1399e16  # [1/s]
        dict_out["b"] = 5.87e-6  # [1/K]
        dict_out["c"] = 7.544e-9  # [1/K^2]
    if molecule == "N2":
        dict_out["groundPolarizability"] = 1.7406e-30  # [m^3]
        dict_out["groundFrequency"] = 2.6049e16  # [1/s]
        dict_out["b"] = 1.8e-6  # [1/K]
        dict_out["c"] = 0.0  # [1/K^2]
    if molecule == "O2":
        dict_out["groundPolarizability"] = 1.5658e-30  # [m^3]
        dict_out["groundFrequency"] = 2.1801e16  # [1/s]
        dict_out["b"] = -2.369e-6  # [1/K]
        dict_out["c"] = 8.687e-9  # [1/K^2]
    if molecule == "Air":
        dict_out["groundPolarizability"] = 1.6970e-30  # [m^3]
        dict_out["groundFrequency"] = 2.47044e16  # [1/s]
        dict_out["b"] = 10.6e-6  # [1/K]
        dict_out["c"] = 7.909e-9  # [1/K^2]

    return dict_out


def spectroscopy_constants(molecule: str) -> dict[str, float]:
    """
    Returns spectroscopy constants

    Parameters:
        molecule: NO+, N2+, O2+, NO, N2, O2, H2

    Returns:
        dict: A dictionary containing spectroscopy constants
            - omega_e: vibrational constant - first term in [cm^-1]
            - omega_xe: vibrational constant – second term in [cm^-1]
            - omega_ye: vibrational constant – third term in [cm^-1]
            - B_e: equilibrium rotational constant in [cm^-1]
            - alpha_e: rotational constant in [cm^-1]
            - D_e: centrifugal distortion in [cm^-1]
            - r_e: internuclear distance in [m]

    Reference:
        NO+ https://webbook.nist.gov/cgi/cbook.cgi?Name=NO%2B&Units=SI&cDI=on

        N2+ https://webbook.nist.gov/cgi/cbook.cgi?Name=N2%2B&Units=SI&cDI=on

        O2+ https://webbook.nist.gov/cgi/cbook.cgi?Name=O2%2B&Units=SI&cDI=on

        NO https://webbook.nist.gov/cgi/cbook.cgi?Name=NO&Units=SI&cDI=on

        N2 https://webbook.nist.gov/cgi/cbook.cgi?Name=N2&Units=SI&cDI=on

        O2 https://webbook.nist.gov/cgi/cbook.cgi?Name=O2&Units=SI&cDI=on

        H2 https://webbook.nist.gov/cgi/cbook.cgi?Name=H2&Units=SI&cDI=on
    """
    dict_out = {}

    if molecule == "NO+":
        dict_out["omega_e"] = 2376.72
        dict_out["omega_xe"] = 16.255
        dict_out["omega_ye"] = -0.01562
        dict_out["B_e"] = 1.997195
        dict_out["alpha_e"] = 0.018790
        dict_out["D_e"] = 6.64e-6
        dict_out["r_e"] = 1.06322e-10

    if molecule == "N2+":
        dict_out["omega_e"] = 2207.0115
        dict_out["omega_xe"] = 16.0616
        dict_out["omega_ye"] = -0.04289
        dict_out["B_e"] = 1.93176
        dict_out["alpha_e"] = 0.0181
        dict_out["D_e"] = 6.10e-6
        dict_out["r_e"] = 1.11642e-10

    if molecule == "O2+":
        dict_out["omega_e"] = 1905.892
        dict_out["omega_xe"] = 16.489
        dict_out["omega_ye"] = 0.02057
        dict_out["B_e"] = 1.689824
        dict_out["alpha_e"] = 0.019363
        dict_out["D_e"] = 5.32e-6
        dict_out["r_e"] = 1.1164e-10

    if molecule == "NO":
        dict_out["omega_e"] = 1904.1346
        dict_out["omega_xe"] = 14.08836
        dict_out["omega_ye"] = 0.01005
        dict_out["B_e"] = 1.704885
        dict_out["alpha_e"] = 0.0175416
        dict_out["D_e"] = 0.54e-6
        dict_out["r_e"] = 1.15077e-10

    if molecule == "N2":
        dict_out["omega_e"] = 2358.57
        dict_out["omega_xe"] = 14.324
        dict_out["omega_ye"] = -2.26e-3
        dict_out["B_e"] = 1.998241
        dict_out["alpha_e"] = 0.017318
        dict_out["D_e"] = 5.76e-6
        dict_out["r_e"] = 1.2126e-10

    if molecule == "O2":
        dict_out["omega_e"] = 1580.161
        dict_out["omega_xe"] = 11.95127
        dict_out["omega_ye"] = 0.0458489
        dict_out["B_e"] = 1.44562
        dict_out["alpha_e"] = 0.0159305
        dict_out["D_e"] = 4.839e-6
        dict_out["r_e"] = 1.20752e-10

    if molecule == "H2":
        dict_out["omega_e"] = 4401.21
        dict_out["omega_xe"] = 121.33
        dict_out["omega_ye"] = 0.0
        dict_out["B_e"] = 60.853
        dict_out["alpha_e"] = 3.062
        dict_out["D_e"] = 0.0471
        dict_out["r_e"] = 0.74144

    return dict_out

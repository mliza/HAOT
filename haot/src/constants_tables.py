"""
    Date:   10/26/2024
    Author: Martin E. Liza
    File:   constants_tables.py
    Def:    This file contains spectroscopy data used across the package
"""


def sutherland_constants(gas="Air"):
    """
    Constants used for Sutherland's law of viscosity and thermal conductivity

    Reference:
        https://doc.comsol.com/5.5/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.27.html

    Parameters:
        gas (str): Air (default), Argon, N2, O2

    Returns:
        dictionary (dict):
            temperature_ref (float): reference temperature [K]
            viscosity_ref (float): reference viscosity [Ns/m2]
            sutherland_visc (float): Sutherland viscosity constant [K]
            conductivity_ref (float): reference thermal conductivity [W/mK]
            sutherland_cond (float): Sutherland thermal conductivity [K]

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


def karl_2003():
    """
    Constants used in the Karl method

    Reference:
        https://arc.aiaa.org/doi/pdf/10.2514/6.2003-4252

    Returns:
       dict_out (dict):
    """

    dict_out = {
        "N": 0.301e-3,
        "O": 0.182e-3,
        "NO": 0.221e-3,
        "N2": 0.238e-3,
        "O2": 0.190e-3,
    }  # [m3/kg]
    return dict_out


# Polarizability
def polarizability():
    dict_out = {
        "N+": 0.559e-30,
        "O+": 0.345e-30,
        "NO+": 1.021e-30,
        "N2+": 2.386e-30,
        "O2+": 0.238e-30,
        "N": 1.100e-30,
        "O": 0.802e-30,
        "NO": 1.700e-30,
        "N2": 1.7403e-30,
        "O2": 1.5689e-30,
    }  # [m^3]
    return dict_out


def polarizability_derivatives(molecule="N2"):
    # https://link.springer.com/content/pdf/10.1134/BF03355985.pdf
    """
    Constants from Buldakov Paper:
    "Temperature Dependence of Polarizability
    of Diatomic Homonuclear Molecules."

    References:
        - Each constant's reference is provided in the respective comments or documentation.
    """
    dict_out = {}
    if molecule == "H2":
        dict_out["zeroth"] = (
            0.7849e-30  # https://edisciplinas.usp.br/pluginfile.php/4557662/mod_resource/content/1/CRC%20Handbook%20of%20Chemistry%20and%20Physics%2095th%20Edition.pdf
        )
        dict_out["first"] = (
            0.90e-30  # https://www.tandfonline.com/doi/abs/10.1080/00268978000103191
        )
        dict_out["second"] = (
            0.49e-30  # https://www.tandfonline.com/doi/abs/10.1080/00268978000103191
        )
        dict_out["third"] = (
            -0.85e-30
        )  # https://www.tandfonline.com/doi/abs/10.1080/00268978000103191

    if molecule == "N2":
        dict_out["zeroth"] = (
            1.7801e-30  # https://edisciplinas.usp.br/pluginfile.php/4557662/mod_resource/content/1/CRC%20Handbook%20of%20Chemistry%20and%20Physics%2095th%20Edition.pdf
        )
        dict_out["first"] = (
            1.86e-30  # M. A. Buldakov, B. V. Korolev, I. I. Matrosov, and T. N. Popova, Opt. Spektrosk. 63, 775 (1987) MISSING REFERENCE
        )
        dict_out["second"] = (
            1.2e-30  # M. A. Buldakov, B. V. Korolev, I. I. Matrosov, and T. N. Popova, Opt. Spektrosk. 63, 775 (1987) MISSING REFERENCE
        )
        dict_out["third"] = (
            -4.6e-30
        )  # https://pubs.aip.org/aip/jcp/article-abstract/78/9/5287/777120/Theoretical-study-of-the-effects-of-vibrational?redirectedFrom=fulltext

    if molecule == "O2":
        dict_out["zeroth"] = (
            1.6180e-30  # https://edisciplinas.usp.br/pluginfile.php/4557662/mod_resource/content/1/CRC%20Handbook%20of%20Chemistry%20and%20Physics%2095th%20Edition.pdf
        )
        dict_out["first"] = (
            1.76e-30  # M. A. Buldakov, B. V. Korolev, I. I. Matrosov, and T. N. Popova, Opt. Spektrosk. 63, 775 (1987) MISSING REFERENCE
        )
        dict_out["second"] = (
            3.4e-30  # M. A. Buldakov, B. V. Korolev, I. I. Matrosov, and T. N. Popova, Opt. Spektrosk. 63, 775 (1987) MISSING REFERENCE
        )
        dict_out["third"] = (
            -23.7e-30
        )  # https://pubs.aip.org/aip/jcp/article-abstract/100/2/1297/482621/Frequency-dependent-polarizabilities-of-O2-and-van?redirectedFrom=fulltext

    return dict_out  # [m^3]


def kerl_interpolation(molecule="N2"):
    # https://onlinelibrary.wiley.com/doi/10.1002/bbpc.19920960517
    # Check reference in paper
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
        dict_out["c"] = 0.0
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


def spectroscopy_constants(molecule="N2"):
    # https://doi.org/10.1063/1.2436891
    # minimun electronic energy T_e [cm^-1]
    # harmonic frequency (omega_e)  [cm^-1]
    # first anharmonic correction (omega_e x_e) [cm^-1]
    # (omega_e y_e)
    # equilibrium rotational constant (B_e)
    # anharmonic correction to the rotational constant (alpha_e) [cm^-1]
    # centrifugal distortion constant (D_e) [cm^-1]
    # binding energy (D_o) [eV]
    # equilibrium internuclear distance(R_e) [\r{A}]
    # ionization potential (IP) [eV]

    dict_out = {}

    # https://webbook.nist.gov/cgi/cbook.cgi?Name=NO%2B&Units=SI&cDI=on
    if molecule == "NO+":
        dict_out["omega_e"] = 2376.72
        dict_out["omega_xe"] = 16.255
        dict_out["omega_ye"] = -0.01562
        dict_out["B_e"] = 1.997195
        dict_out["alpha_e"] = 0.018790
        dict_out["D_e"] = 6.64e-6
        dict_out["r_e"] = 1.06322e-10  # [m]

    # https://webbook.nist.gov/cgi/cbook.cgi?Name=N2%2B&Units=SI&cDI=on
    if molecule == "N2+":
        dict_out["omega_e"] = 2207.0115
        dict_out["omega_xe"] = 16.0616
        dict_out["omega_ye"] = -0.04289
        dict_out["B_e"] = 1.93176
        dict_out["alpha_e"] = 0.0181
        dict_out["D_e"] = 6.10e-6
        dict_out["r_e"] = 1.11642e-10

    # https://webbook.nist.gov/cgi/cbook.cgi?Name=O2%2B&Units=SI&cDI=on
    if molecule == "O2+":
        dict_out["omega_e"] = 1905.892
        dict_out["omega_xe"] = 16.489
        dict_out["omega_ye"] = 0.02057
        dict_out["B_e"] = 1.689824
        dict_out["alpha_e"] = 0.019363
        dict_out["D_e"] = 5.32e-6
        dict_out["r_e"] = 1.1164e-10

    # https://webbook.nist.gov/cgi/cbook.cgi?Name=NO&Units=SI&cDI=on
    if molecule == "NO":
        dict_out["omega_e"] = 1904.1346
        dict_out["omega_xe"] = 14.08836
        dict_out["omega_ye"] = 0.01005
        dict_out["B_e"] = 1.704885
        dict_out["alpha_e"] = 0.0175416
        dict_out["D_e"] = 0.54e-6
        dict_out["r_e"] = 1.15077e-10

    # https://webbook.nist.gov/cgi/cbook.cgi?Name=N2&Units=SI&cDI=on
    if molecule == "N2":
        dict_out["omega_e"] = 2358.57
        dict_out["omega_xe"] = 14.324
        dict_out["omega_ye"] = -0.00226
        dict_out["B_e"] = 1.998241
        dict_out["alpha_e"] = 0.017318
        dict_out["D_e"] = 5.9e-6
        dict_out["r_e"] = 1.2126e-10

    # https://webbook.nist.gov/cgi/cbook.cgi?Name=O2&Units=SI&cDI=on
    if molecule == "O2":
        dict_out["omega_e"] = 1580.161
        dict_out["omega_xe"] = 11.95127
        dict_out["omega_ye"] = 0.0458489
        dict_out["B_e"] = 1.44562
        dict_out["alpha_e"] = 0.0159305
        dict_out["D_e"] = 4.839e-6
        dict_out["r_e"] = 1.20752e-10

    return dict_out  # wavenumber units  #[cm^-1]

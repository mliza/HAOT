'''
    Date:   10/26/2024
    Author: Martin E. Liza
    File:   constant_tables.py
    Def:    This file contains spectroscopy data used across the package
'''

# Gladstone-Dale constants 
def karl_2003(): #follows SU2 MutationPP format 
# https://arc.aiaa.org/doi/pdf/10.2514/6.2003-4252
    dict_out = { 'N'  : 0.301E-3,
                 'O'  : 0.182E-3,
                 'NO' : 0.221E-3,
                 'N2' : 0.238E-3,
                 'O2' : 0.190E-3 } #[m3/kg]
    return dict_out

# Polarizability
def polarizability(): #follows SU2 MutationPP format
    dict_out = { #'e-'  : 0.0, 
                 'N+'   : 0.559E-30,
                 'O+'   : 0.345E-30,
                 'NO+'  : 1.021E-30,
                 'N2+'  : 2.386E-30,
                 'O2+'  : 0.238E-30,
                 'N'    : 1.100E-30,
                 'O'    : 0.802E-30,
                 'NO'   : 1.700E-30,
                 'N2'   : 1.7403E-30,
                 'O2'   : 1.5689E-30 } # [m^3]  
    return dict_out 

    
def polarizability_derivatives(molecule='N2'):
    # https://link.springer.com/content/pdf/10.1134/BF03355985.pdf
    """
    Constants below come from Buldakov Paper: 'Temperature Dependence of Polarizability
of Diatomic Homonuclear Molecules'. Reference for each constant are provided
below
    """
    dict_out = { }
    if molecule == 'H2':
        dict_out['zeroth'] = 0.7849E-30 #https://edisciplinas.usp.br/pluginfile.php/4557662/mod_resource/content/1/CRC%20Handbook%20of%20Chemistry%20and%20Physics%2095th%20Edition.pdf
        dict_out['first'] = 0.90E-30 #https://www.tandfonline.com/doi/abs/10.1080/00268978000103191
        dict_out['second'] = 0.49E-30 #https://www.tandfonline.com/doi/abs/10.1080/00268978000103191 
        dict_out['third'] = -0.85E-30 #https://www.tandfonline.com/doi/abs/10.1080/00268978000103191

    if molecule == 'N2':
        dict_out['zeroth'] = 1.7801E-30 #https://edisciplinas.usp.br/pluginfile.php/4557662/mod_resource/content/1/CRC%20Handbook%20of%20Chemistry%20and%20Physics%2095th%20Edition.pdf
        dict_out['first'] = 1.86E-30 #M. A. Buldakov, B. V. Korolev, I. I. Matrosov, and T. N. Popova, Opt. Spektrosk. 63, 775 (1987) MISSING REFERENCE 
        dict_out['second'] = 1.2E-30 #M. A. Buldakov, B. V. Korolev, I. I. Matrosov, and T. N. Popova, Opt. Spektrosk. 63, 775 (1987) MISSING REFERENCE 
        dict_out['third'] = -4.6E-30 #https://pubs.aip.org/aip/jcp/article-abstract/78/9/5287/777120/Theoretical-study-of-the-effects-of-vibrational?redirectedFrom=fulltext

    if molecule == 'O2':
        dict_out['zeroth'] = 1.6180E-30 #https://edisciplinas.usp.br/pluginfile.php/4557662/mod_resource/content/1/CRC%20Handbook%20of%20Chemistry%20and%20Physics%2095th%20Edition.pdf
        dict_out['first'] = 1.76E-30 #M. A. Buldakov, B. V. Korolev, I. I. Matrosov, and T. N. Popova, Opt. Spektrosk. 63, 775 (1987) MISSING REFERENCE 
        dict_out['second'] = 3.4E-30 #M. A. Buldakov, B. V. Korolev, I. I. Matrosov, and T. N. Popova, Opt. Spektrosk. 63, 775 (1987) MISSING REFERENCE 
        dict_out['third'] = -23.7E-30 #https://pubs.aip.org/aip/jcp/article-abstract/100/2/1297/482621/Frequency-dependent-polarizabilities-of-O2-and-van?redirectedFrom=fulltext

    return dict_out #[m^3]

def kerl_interpolation(molecule='N2'):
    #https://onlinelibrary.wiley.com/doi/10.1002/bbpc.19920960517
    # Check reference in paper
    dict_out = { }
    if molecule == 'H2':
        dict_out['groundPolarizability'] = 0.80320E-30 #[m^3]
        dict_out['groundFrequency'] = 2.1399E16 # [1/s]
        dict_out['b'] = 5.87E-6 # [1/K]
        dict_out['c'] = 7.544E-9 # [1/K^2]

    if molecule == 'N2':
        dict_out['groundPolarizability'] = 1.7406E-30 #[m^3]
        dict_out['groundFrequency'] = 2.6049E16 # [1/s]
        dict_out['b'] = 1.8E-6 # [1/K]
        dict_out['c'] = 0.0 
    if molecule == 'O2':
        dict_out['groundPolarizability'] = 1.5658E-30 #[m^3]
        dict_out['groundFrequency'] = 2.1801E16 # [1/s]
        dict_out['b'] = -2.369E-6 # [1/K]
        dict_out['c'] = 8.687E-9# [1/K^2]
    if molecule == 'Air':
        dict_out['groundPolarizability'] = 1.6970E-30 #[m^3]
        dict_out['groundFrequency'] = 2.47044E16 # [1/s]
        dict_out['b'] = 10.6E-6 # [1/K]
        dict_out['c'] = 7.909E-9# [1/K^2]

    return dict_out

def spectroscopy_constants(molecule='N2'):
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

    dict_out = { }

    #https://webbook.nist.gov/cgi/cbook.cgi?Name=NO%2B&Units=SI&cDI=on
    if molecule == 'NO+':
        dict_out['omega_e'] = 2376.72
        dict_out['omega_xe'] = 16.255
        dict_out['omega_ye'] = -0.01562
        dict_out['B_e'] = 1.997195
        dict_out['alpha_e'] = 0.018790
        dict_out['D_e'] = 6.64E-6
        dict_out['r_e'] = 1.06322E-10 #[m] 

    #https://webbook.nist.gov/cgi/cbook.cgi?Name=N2%2B&Units=SI&cDI=on
    if molecule == 'N2+':
        dict_out['omega_e'] = 2207.0115
        dict_out['omega_xe'] = 16.0616
        dict_out['omega_ye'] = -0.04289
        dict_out['B_e'] = 1.93176
        dict_out['alpha_e'] = 0.0181
        dict_out['D_e'] = 6.10E-6
        dict_out['r_e'] = 1.11642E-10 

    #https://webbook.nist.gov/cgi/cbook.cgi?Name=O2%2B&Units=SI&cDI=on
    if molecule == 'O2+':
        dict_out['omega_e'] = 1905.892
        dict_out['omega_xe'] = 16.489
        dict_out['omega_ye'] = 0.02057
        dict_out['B_e'] = 1.689824
        dict_out['alpha_e'] = 0.019363
        dict_out['D_e'] = 5.32E-6
        dict_out['r_e'] = 1.1164E-10 

    #https://webbook.nist.gov/cgi/cbook.cgi?Name=NO&Units=SI&cDI=on
    if molecule == 'NO':
        dict_out['omega_e'] = 1904.1346
        dict_out['omega_xe'] = 14.08836
        dict_out['omega_ye'] = 0.01005
        dict_out['B_e'] = 1.704885
        dict_out['alpha_e'] = 0.0175416
        dict_out['D_e'] = 0.54E-6 
        dict_out['r_e'] = 1.15077E-10 

    #https://webbook.nist.gov/cgi/cbook.cgi?Name=N2&Units=SI&cDI=on
    if molecule == 'N2':
        dict_out['omega_e'] = 2358.57
        dict_out['omega_xe'] = 14.324
        dict_out['omega_ye'] = -0.00226
        dict_out['B_e'] = 1.998241
        dict_out['alpha_e'] = 0.017318
        dict_out['D_e'] = 5.9E-6
        dict_out['r_e'] = 1.2126E-10 

    #https://webbook.nist.gov/cgi/cbook.cgi?Name=O2&Units=SI&cDI=on
    if molecule == 'O2':
        dict_out['omega_e'] = 1580.161
        dict_out['omega_xe'] = 11.95127
        dict_out['omega_ye'] = 0.0458489
        dict_out['B_e'] = 1.44562
        dict_out['alpha_e'] = 0.0159305
        dict_out['D_e'] = 4.839E-6 
        dict_out['r_e'] = 1.20752E-10 

    return dict_out #wavenumber units  #[cm^-1]

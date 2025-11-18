"""
Date:   11/11/2024
Author: Martin E. Liza
File:   coordinates.py
Def:    Contains coordinates functions.
"""

import numpy as np


def _earth_data() -> dict:
    """
    Returns Earth geometric parameters.
    Source: https://en.wikipedia.org/wiki/Earth_radius
    """
    semi_minor_earth_radius_m = 6356752.3142
    semi_major_earth_radius_m = 6378137.0
    eccentricity = np.sqrt(
        1 - (semi_minor_earth_radius_m / semi_major_earth_radius_m) ** 2
    )
    dict_out = {
        "semi_minor_earth_radius_m": semi_minor_earth_radius_m,
        "semi_major_earth_radius_m": semi_major_earth_radius_m,
        "earth_eccentricity": eccentricity,
    }

    return dict_out


def haversine_distance(
    latitude_1: float,
    longitude_1: float,
    latitude_2: float,
    longitude_2: float,
    sphere_radius: float = 6371e3,
) -> float:
    """
    Calculates the Haversine distance between two points

    Parameters:
        latitude_1: latitude of point 1 in [degs]
        longitude_1: longitude of point 1 in [degs]
        latitude_2: latitude of point 2 in [degs]
        longitude_2: longitude of point 2 in [degs]
        sphere_radius: sphere radius, 6371E3 (default, earth radius in [m])

    Returns:
        Haversine distance in units of sphere radius
    """

    # Lat and Lon distance in rads and interest lat and lon
    lat_distance = np.deg2rad(latitude_2 - latitude_1)
    lon_distance = np.deg2rad(longitude_2 - longitude_1)
    lat_1 = np.deg2rad(latitude_1)
    lat_2 = np.deg2rad(latitude_2)

    tmp = (
        1
        - np.cos(lat_distance)
        + np.cos(lat_1) * np.cos(lat_2) * (1 - np.cos(lon_distance))
    ) / 2
    return 2 * sphere_radius * np.asin(np.sqrt(tmp))


def lla_to_ecef(*args) -> np.ndarray:
    """
    Converts LLA position to ECEF position

    Parameters:
        *args:
            Either of the following forms:
                (1) lat_deg, lon_deg, alt_m : float
                    Latitude [degrees], longitude [degrees], and altitude [m].
                (2) [lat_deg, lon_deg, alt_m] : array_like
                    Iterable (list, tuple, or NumPy array) containing latitude [degrees],
                    longitude [degrees], and altitude [m].
    Returns:
        ecef_x, ecef_y, ecef_z: float ECEF coordinates [m]
    """
    if len(args) == 1:
        lat_degs, lon_deg, alt_m = args[0]
    elif len(args) == 3:
        lat_deg, lon_deg, alt_m = args
    else:
        raise ValueError("Pass either (lat, lon, alt) or a vector [lat, lon, alt].")

    c_lat = np.cos(np.deg2rad(lat_deg))
    c_lon = np.cos(np.deg2rad(lon_deg))
    s_lat = np.sin(np.deg2rad(lat_deg))
    s_lon = np.sin(np.deg2rad(lon_deg))

    earth_dict = _earth_data()
    radius_of_curvature = earth_dict["semi_major_earth_radius_m"] / np.sqrt(
        1.0 - (earth_dict["earth_eccentricity"] * s_lat) ** 2
    )

    ecef_x = (radius_of_curvature + alt_m) * c_lat * c_lon
    ecef_y = (radius_of_curvature + alt_m) * c_lat * s_lon
    ecef_z = (
        (1.0 - earth_dict["earth_eccentricity"] ** 2) * radius_of_curvature + alt_m
    ) * s_lat

    return np.array([ecef_x, ecef_y, ecef_z])


def ecef_to_lla(*args) -> np.ndarray:
    """
    Converts ECEF position to LLA position

    Source: https://danceswithcode.net/engineeringnotes/geodetic_to_ecef/geodetic_to_ecef.html

    Parameters:
        *args:
            Either of the following forms:
                (1) ecef_x, ecef_y, ecef_z : float
                     ECEF Cartesian coordinates [m].
                (2) [ecef_x, ecef_y, ecef_z] : array_like
                    Iterable (list, tuple, or NumPy array) containing
                    ECEF coordinates [x, y, z] in meters.

    Returns:
        lla : numpy.ndarray of shape (3,)
            Geodetic coordinates [latitude (deg), longitude (deg), altitude (m)].
    """
    if len(args) == 1:
        ecef_x, ecef_y, ecef_z = args[0]
    elif len(args) == 3:
        ecef_x, ecef_y, ecef_z = args
    else:
        raise ValueError("Pass either (lat, lon, alt) or a vector [lat, lon, alt].")

    earth_dict = _earth_data()
    eccentricity_2 = earth_dict["earth_eccentricity"] ** 2
    a_1 = earth_dict["semi_major_earth_radius_m"] * eccentricity_2
    a_2 = a_1**2
    a_3 = 0.5 * a_1 * eccentricity_2
    a_4 = 2.5 * a_2
    a_5 = a_1 + a_3
    a_6 = 1 - eccentricity_2
    radius = np.linalg.norm([ecef_x, ecef_y, ecef_z])
    u = a_2 / radius
    v = a_3 - (a_4 / radius)
    s_2 = (ecef_z / radius) ** 2
    c_2 = (ecef_x**2 + ecef_y**2) / radius**2
    longitude_deg = np.degrees(np.arctan(ecef_y / ecef_x))

    if c_2 > 0.3:
        s = (np.abs(ecef_z) / radius) * (1.0 + c_2 * (a_1 + u + s_2 * v) / radius)
        latitude_rad = np.arcsin(s)
        c = np.sqrt(1.0 - s**2)
    else:
        c = np.sqrt(c_2) * (1 - s_2 * (a_5 - u - c_2 * v) / radius)
        latitude_rad = np.arccos(c)
        s = np.sqrt(1.0 - c**2)

    g = 1.0 - eccentricity_2 * s**2
    rg = earth_dict["semi_major_earth_radius_m"] / np.sqrt(g)
    rf = a_6 * rg
    u = np.sqrt(c_2) * radius - (rg * c)
    v = np.abs(ecef_z) - (rf * s)
    f = (c * u) + (s * v)
    m = (c * v) - (s * u)
    p = m / (rf / g + f)
    latitude_deg = np.degrees(latitude_rad + p)
    altitude_m = f + (0.5 * m * p)
    if ecef_z < 0:
        latitude_deg *= -1
    return np.array([latitude_deg, longitude_deg, altitude_m])


def roll_dcm(roll_deg: float) -> np.ndarray:
    """
    Creates a roll DCM

    Parameters:
        roll_deg: roll angle in [degs]

    Returns: roll DCM
    """
    roll_rad = np.deg2rad(roll_deg)
    dcm = np.zeros((3, 3))
    dcm[0][0] = 1.0
    dcm[1][1] = np.cos(roll_rad)
    dcm[1][2] = np.sin(roll_rad)
    dcm[2][1] = -np.sin(roll_rad)
    dcm[2][2] = np.cos(roll_rad)

    return dcm


def pitch_dcm(pitch_deg: float) -> np.ndarray:
    """
    Creates a pitch DCM

    Parameters:
        pitch_deg: pitch angle in [degs]

    Returns: pitch DCM
    """
    pitch_rad = np.deg2rad(pitch_deg)
    dcm = np.zeros((3, 3))
    dcm[0][0] = np.cos(pitch_rad)
    dcm[0][2] = -np.sin(pitch_rad)
    dcm[1][1] = 1.0
    dcm[2][0] = np.sin(pitch_rad)
    dcm[2][2] = np.cos(pitch_rad)

    return dcm


def yaw_dcm(yaw_deg: float) -> np.ndarray:
    """
    Creates a yaw DCM

    Parameters:
        yaw_deg: yaw angle in [degs]

    Returns: yaw DCM
    """
    yaw_rad = np.deg2rad(yaw_deg)
    dcm = np.zeros((3, 3))
    dcm[0][0] = np.cos(yaw_rad)
    dcm[0][1] = np.sin(yaw_rad)
    dcm[1][0] = -np.sin(yaw_rad)
    dcm[1][1] = np.cos(yaw_rad)
    dcm[2][2] = 1.0

    return dcm


def euler_321_dcm(*args) -> np.ndarray:
    """
    Parameters:
        *args:
            Either of the following forms:
                (1) roll_deg, pitch_deg, yaw_deg : float
                (2) [roll_deg, pitch_deg, yaw_deg] : array_like
    Returns:
        standard 321 (Yaw->Pitch->Roll) DCM
        ecef_x, ecef_y, ecef_z: float ECEF coordinates [m]
    """
    if len(args) == 1:
        roll_deg, pitch_deg, yaw_deg = args[0]
    elif len(args) == 3:
        roll_deg, pitch_deg, yaw_deg = args
    else:
        raise ValueError(
            "Pass either (roll, pitch, yaw) or a vector [roll, pitch, yaw]."
        )

    dcm = np.zeros((3, 3))
    c_roll = np.cos(np.deg2rad(roll_deg))
    s_roll = np.sin(np.deg2rad(roll_deg))
    c_pitch = np.cos(np.deg2rad(pitch_deg))
    s_pitch = np.sin(np.deg2rad(pitch_deg))
    c_yaw = np.cos(np.deg2rad(yaw_deg))
    s_yaw = np.sin(np.deg2rad(yaw_deg))

    dcm[0][0] = c_pitch * c_yaw
    dcm[0][1] = c_pitch * s_yaw
    dcm[0][2] = -s_pitch
    dcm[1][0] = s_roll * s_pitch * c_yaw - c_roll * s_yaw
    dcm[1][1] = s_roll * s_pitch * s_yaw + c_roll * c_yaw
    dcm[1][2] = s_roll * c_pitch
    dcm[2][0] = c_roll * s_pitch * c_yaw + s_roll * s_yaw
    dcm[2][1] = c_roll * s_pitch * s_yaw - s_roll * c_yaw
    dcm[2][2] = c_roll * c_pitch

    return dcm


def euler_angles_from_dcm(dcm: np.ndarray) -> np.ndarray:
    """
    Parameters:
        dcm: directional cosine matrix
    Returns:
        euler angles, [roll_deg, pitch_deg, yaw_deg] in [degrees]
    """

    roll_deg = np.rad2deg(np.arctan(dcm[1][2] / dcm[2][2]))
    pitch_deg = np.rad2deg(np.arcsin(-dcm[0][2]))
    yaw_deg = np.rad2deg(np.arctan(dcm[0][1] / dcm[0][0]))

    return np.array([roll_deg, pitch_deg, yaw_deg])

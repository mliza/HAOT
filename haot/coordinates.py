import numpy as np


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

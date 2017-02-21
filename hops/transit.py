__all__ = ['transit']


import numpy as np

from transit_flux_drop import *
from exoplanet_orbit import *


def transit(limb_darkening_coefficients, rp_over_rs,
            period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array,
            method='claret', precision=3):

    position_vector = exoplanet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return transit_flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance,
                             method=method, precision=precision)

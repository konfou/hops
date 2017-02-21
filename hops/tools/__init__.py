import numpy as np

from scipy.optimize import curve_fit

import matplotlib as plt

import ephem

import pyfits as pf


def cartesian_to_polar(x_position, y_position, x_ref_position, y_ref_position):

    x_position, y_position = float(x_position), float(y_position)
    x_ref_position, y_ref_position = float(x_ref_position), float(y_ref_position)

    radius = np.sqrt((x_position - x_ref_position) ** 2 + (y_position - y_ref_position) ** 2)

    if (x_position - x_ref_position) > 0:
        if (y_position - y_ref_position) >= 0:
            angle = np.arctan((y_position - y_ref_position) / (x_position - x_ref_position))
        else:
            angle = 2.0 * np.pi + np.arctan((y_position - y_ref_position) / (x_position - x_ref_position))
    elif (x_position - x_ref_position) < 0:
        angle = np.arctan((y_position - y_ref_position) / (x_position - x_ref_position)) + np.pi
    else:
        if (y_position - y_ref_position) >= 0:
            angle = np.pi / 2
        else:
            angle = -np.pi / 2

    return radius, angle


def jd_to_hjd(ra_target, dec_target, julian_date, degrees=False):

    if degrees:
        ra_target *= np.pi / 180
        dec_target *= np.pi / 180
    else:
        k = ephem.Equatorial(ra_target, dec_target)
        ra_target = k.ra
        dec_target = k.dec

    sun = ephem.Sun()
    sun.compute(ephem.date(julian_date - 2415020))
    ra_sun, dec_sun = float(sun.ra), float(sun.dec)

    a = 149597870700.0 / ephem.c
    b = np.sin(dec_target) * np.sin(dec_sun)
    c = np.cos(dec_target) * np.cos(dec_sun) * np.cos(ra_target - ra_sun)

    heliocentric_julian_date = julian_date - (a * (b + c)) / (24.0 * 60.0 * 60.0)

    return heliocentric_julian_date


def gaussian(x_array, model_norm, model_floor, model_mean, model_std):

    return model_floor + ((model_norm / abs(model_std) / np.sqrt(2.0 * np.pi)) *
                          np.exp(- 0.5 * (model_mean - x_array) * (model_mean - x_array) / (model_std * model_std)))


def fit_gaussian(data_x_array, data_y_array):

    mean = data_x_array[np.argmax(data_y_array)]
    std = np.abs(data_x_array[np.argmin(np.abs(data_y_array - np.max(data_y_array) / 2))] - mean)
    norm = (np.max(data_y_array) - np.min(data_y_array)) / gaussian(mean, 1.0, 0.0, mean, std)
    floor = np.min(data_y_array)

    try:

        [norm, floor, mean, std], covariance = \
            curve_fit(gaussian, data_x_array, data_y_array, p0=[norm, floor, mean, std])

        return norm, floor, mean, std

    except RuntimeError:

        print 'fit_gaussian: could not find a Gaussian'

        return np.nan, np.nan, np.nan, np.nan


def distribution1d(data_array, step=None, binning=None):

    data_array = np.array(data_array).flatten()
    data_array = np.sort(data_array)

    if binning:
        start = np.mod(data_array.size, binning)
        data_array = np.mean(np.reshape(data_array[start:], (data_array.size / binning, binning)), 1)

    if not step:
        step = np.sqrt(np.median((data_array - np.median(data_array)) ** 2)) / 5.0

    min_value = min(data_array)
    max_value = max(data_array)
    bins_number = int((max_value - min_value) / step) + 2

    bins = min_value + step / 2. + np.arange(bins_number) * step
    bins_i = bins - step / 2
    bins_f = bins + step / 2
    counts = np.bincount(np.int_((data_array - min_value) / step))
    counts = np.insert(counts, len(counts), np.zeros(bins_number - len(counts)))

    return bins_i, bins_f, counts


def fit_distribution1d_gaussian(data_array, step=None, binning=None):

    bins_i, bins_f, counts = distribution1d(data_array, step=step, binning=binning)

    norm, floor, mean, std = fit_gaussian(0.5 * (bins_i + bins_f), counts)

    return norm, floor, mean, std


def fits_like(fits):
    copy_fits = [pf.PrimaryHDU(header=fits[0].header, data=fits[0].data)]
    for j in fits[1:]:
        copy_fits.append(pf.ImageHDU(header=j.header, data=j.data))
    return pf.HDUList(copy_fits)


def fit_2d_gauss(data_array, predicted_x_mean=None, predicted_y_mean=None, search_window=None):

    found = False

    def function_2d_gauss((x_array, y_array), model_norm, model_floor,
                          model_x_mean, model_y_mean, model_x_std, model_y_std, model_theta):

        a = (np.cos(model_theta) ** 2) / (2 * model_x_std ** 2) + (np.sin(model_theta) ** 2) / (2 * model_y_std ** 2)
        b = -(np.sin(2 * model_theta)) / (4 * model_x_std ** 2) + (np.sin(2 * model_theta)) / (4 * model_y_std ** 2)
        c = (np.sin(model_theta) ** 2) / (2 * model_x_std ** 2) + (np.cos(model_theta) ** 2) / (2 * model_y_std ** 2)

        return (model_floor + model_norm * np.exp(- (a * ((x_array - model_x_mean) ** 2)
                                                  + 2.0 * b * (x_array - model_x_mean) * (y_array - model_y_mean)
                                                  + c * ((y_array - model_y_mean) ** 2)))).flatten()

    if not predicted_x_mean:
        predicted_x_mean = np.where(data_array == np.max(data_array))[1][0]

    if not predicted_y_mean:
        predicted_y_mean = np.where(data_array == np.max(data_array))[0][0]

    if not search_window:
        search_window = max(data_array.shape)

    norm, floor, x_mean, y_mean, x_std, y_std = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    while not found and search_window < max(data_array.shape):

        y_min = int(max(predicted_y_mean - search_window, 0))
        y_max = int(min(predicted_y_mean + search_window, len(data_array)))
        x_min = int(max(predicted_x_mean - search_window, 0))
        x_max = int(min(predicted_x_mean + search_window, len(data_array[0])))

        cropped_data_array = data_array[y_min:y_max, x_min:x_max]
        cropped_x_data_array = np.arange(x_min, x_max)
        cropped_y_data_array = np.arange(y_min, y_max)

        norm = data_array[int(predicted_y_mean)][int(predicted_x_mean)] - np.min(cropped_data_array)
        floor = np.min(cropped_data_array)
        x_mean = cropped_x_data_array[np.where(cropped_data_array == np.max(cropped_data_array))[1][0]]
        y_mean = cropped_y_data_array[np.where(cropped_data_array == np.max(cropped_data_array))[0][0]]
        x_std = max(1.0, np.abs(cropped_x_data_array[np.argmin(
            np.abs(np.sum(cropped_data_array, 0) - np.max(np.sum(cropped_data_array, 0)) / 2))] - x_mean))
        y_std = x_std
        theta = 0.0

        try:
            [norm, floor, x_mean, y_mean, x_std, y_std, theta], covariance = \
                curve_fit(function_2d_gauss,
                          np.meshgrid(cropped_x_data_array, cropped_y_data_array),
                          cropped_data_array.flatten(),
                          p0=[norm, floor, x_mean, y_mean, x_std, y_std, theta])

            x_mean += 0.5
            y_mean += 0.5
            x_std = abs(x_std)
            y_std = abs(y_std)

            if norm > 3 * np.sqrt(covariance[0][0]):
                found = True

            else:
                search_window *= 2.0

        except RuntimeError:

            search_window *= 2.0

    if found:
        return norm, floor, x_mean, y_mean, x_std, y_std

    else:
        print 'fit_2d_gauss: could not find a 2D Gaussian'
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


def find_centroids(data_array, x_low=0, x_upper=None, y_low=0, y_upper=None,
                   x_centre=None, y_centre=None,
                   mean=None, std=None, std_limit=3.0, burn_limit=None, star_std=2,
                   flux_order=False):

    if not x_upper:
        x_upper = data_array.shape[1]

    if not y_upper:
        y_upper = data_array.shape[0]

    if not x_centre:
        x_centre = data_array.shape[1] / 2

    if not y_centre:
        y_centre = data_array.shape[0] / 2

    if not mean or not std:

        fit_norm, fit_floor, fit_mean, fit_std = fit_distribution1d_gaussian(data_array)

        if np.isnan(fit_norm):
            fit_mean = np.mean(data_array)
            fit_std = np.std(data_array)

        if not mean:
            mean = fit_mean

        if not std:
            std = fit_std

    if not burn_limit:
        burn_limit = np.max(data_array)

    data_array = np.full_like(data_array[y_low:y_upper, x_low:x_upper], data_array[y_low:y_upper, x_low:x_upper])

    noise_limit = mean + std_limit * std

    test = []

    for i in range(-star_std, star_std + 1):
        for j in range(-star_std, star_std + 1):
            test.append(np.roll(np.roll(data_array, i, 0), j, 1))

    stars = np.where((data_array < burn_limit) & (np.max(test, 0) == data_array) & (np.min(test, 0) > noise_limit))
    stars = [np.sqrt((stars[1] + x_low - x_centre) ** 2 + (stars[0] + y_low - y_centre) ** 2),
             stars[1] + x_low, stars[0] + y_low, np.sum(test, 0)[stars]]
    stars = np.swapaxes(stars, 0, 1)

    del data_array
    del test

    if not flux_order:
        stars = stars[stars[:, 0].argsort()]
        return stars

    else:
        stars = stars[stars[:, 3].argsort()][::-1]
        return stars

"""GEM Fig. 6b noise distribution with the single-cone efficiencies of Table 1."""

import numpy as np


# Filled marker centres in a Poppler rendering of PDF p. 24 of GEM TN-93-373:
# pdftoppm -f 24 -l 24 -singlefile -scale-to 1550 -png REPORT.pdf figure
# Axis zero, slope and baseline tilt are read from the printed ticks.
PIXEL_X = [
    252,272,291,311,331,351,370,390,410,430,450,470,488,508,528,548,567,588,
    606,626,646,665,685,705,724,744,764,783,802,823,842,862,882,902,922,942,
    961,980,1000,1020,1039,1058,
]
PIXEL_Y = [
    1025,1025,1025,1025,1024,1024,1024,1025,1022,1022,1018,1019,1016,1000,
    968,934,900,729,846,882,901,952,970,981,994,1000,1005,1009,1011,1015,
    1021,1022,1021,1020,1021,1026,1022,1025,1025,1025,1024,1027,
]
MEAN_GEV = 1.48
MARGINS_GEV = np.array([3., 4., 5.])
PUBLISHED_PASS = np.array([.890, .924, .946])


def distribution():
    x, y = np.asarray(PIXEL_X), np.asarray(PIXEL_Y)
    centers = (x - 578.5) / 32.85
    heights = np.maximum(0., 1025. + (x - 250.) * 3. / 820. - y)
    edges = np.r_[centers[0] - .3, .5 * (centers[:-1] + centers[1:]), centers[-1] + .3]
    cdf = np.r_[0., np.cumsum(heights) / heights.sum()]
    return edges - MEAN_GEV, cdf


EDGES, CDF = distribution()
RAW_PASS = np.interp(MARGINS_GEV, EDGES, CDF)


def pass_probability(activity_gev, radius, *, scale_power=1., margin=5., calibrated=True):
    """Probability that activity + mean-subtracted cone noise passes the margin."""
    scale = (np.asarray(radius) / .45) ** scale_power
    raw = np.interp((margin - np.asarray(activity_gev)) / scale, EDGES, CDF)
    if not calibrated:
        return raw
    # Allow for graph resolution and tails outside the displayed range. Rescale
    # the probability in each interval to the tabulated efficiency at its ends.
    return np.interp(raw, np.r_[0., RAW_PASS, 1.], np.r_[0., PUBLISHED_PASS, 1.])

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from math import *

# Returns the centroided x and y coordinates of a region of a fits image (usually a region around the object of interest)
# Parameters are approx x and y center of object, size of region, and fits image file name
def findCentroid(x, y, dimension, fits_image):
    image = fits.getdata(fits_image)
    dim_change = int(np.floor(dimension / 2.))
    sliced_image = image[(y - dim_change):(y + dim_change + 1), (x - dim_change):(x + dim_change + 1)]
    radius = dimension / 2.
    r = dimension
    x_initial = x
    y_initial = y

    # Determine which pixels are only partially in the aperture
    for i in range(r):
        for j in range(r):
            frac = 0
            for a in np.arange((i),(i+1.),(0.01)):
                for b in np.arange((j),(j+1.),(0.01)):
                    if sqrt((a - radius)**2 + (b - radius)**2) <= radius:
                        frac += 1
            fraction = frac / 10000.
            sliced_image[i,j] = (fraction * sliced_image[i,j])
    x_columns = np.sum(sliced_image, 0)
    y_rows = np.sum(sliced_image, 1)

    # Find "x bar" (x coordinate of centroid)
    x_summ = 0
    x_summ_denom = 0
    for i in range(len(x_columns)):
        x_summ += (x_columns[i] * i)
        x_summ_denom += x_columns[i]
    x_bar = x_summ / x_summ_denom
    x_bar_fit = x_bar + (x_initial - np.floor(dimension/2))

    # Find "y bar" (y coordinate of centroid)
    y_summ = 0
    for i in range(len(y_rows)):
        y_summ += y_rows[i] * i
    y_summ_denom = 0
    for y in y_rows:
        y_summ_denom += y
    y_bar = y_summ / y_summ_denom
    y_bar_fit = y_bar + (y_initial - np.floor(dimension/2))
    centroid = [x_bar_fit, y_bar_fit]

    # Find x standard deviation
    x_std = 0
    for i in range(len(x_columns)):
        x_std += ((((x_bar - i)**2)*x_columns[i]) / (x_summ_denom * (x_summ_denom - 1)))
    final_x_std = sqrt(x_std)

    # Find y standard deviation
    y_std = 0
    for i in range(len(y_rows)):
        y_std += ((((y_bar - i)**2)*y_rows[i]) / (y_summ_denom * (y_summ_denom - 1)))
    final_y_std = sqrt(y_std)

    centroid = [x_bar_fit, y_bar_fit]
    uncertainty = [final_x_std, final_y_std]

    return centroid
print(findCentroid(465, 593, 9, "1929 SH.00000025.ENTERED COORDINATES.REDUCED.FIT")) # Put in fits file name

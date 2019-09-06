# See circular centroiding for more a accurate centroids.
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import math

def find_centroid_of_image(fits_image):

    image = fits.getdata(fits_image)


    x_initial = int(input("Enter x: "))
    y_initial = int(input("Enter y "))
    x_dimension = 3#int(input("Enter x dimension: "))
    y_dimension = 3#int(input("Enter x dimension: "))
    dim_change_x = int(np.floor(x_dimension / 2))
    dim_change_y = int(np.floor(y_dimension / 2))
    sliced_image = image[(y_initial - dim_change_y):(y_initial + dim_change_y + 1), (x_initial - dim_change_x):(x_initial + dim_change_x + 1)]
    #Use centroiding method from part a but with standard deviations(uncertainties)
    print(sliced_image)
    x_columns = np.sum(sliced_image, 0)
    y_rows = np.sum(sliced_image, 1)
    #Find x bar (x coordinate of centroid)
    x_summ = 0
    x_summ_denom = 0
    for i in range(len(x_columns)):
        x_summ += (x_columns[i] * i)
        x_summ_denom += x_columns[i]
    x_bar = x_summ / x_summ_denom

    x_bar_fit = x_bar + (x_initial - np.floor(x_dimension/2))
    #Find x standard deviation
    x_std = 0
    for i in range(len(x_columns)):
        x_std += ((((x_bar - i)**2)*x_columns[i]) / (x_summ_denom * (x_summ_denom - 1)))
    final_x_std = math.sqrt(x_std)

    #Find y bar (y coordinate of centroid)
    y_summ = 0
    for i in range(len(y_rows)):
        y_summ += y_rows[i] * i
    y_summ_denom = 0
    for y in y_rows:
        y_summ_denom += y
    y_bar = y_summ / y_summ_denom
    y_bar_fit = y_bar + (y_initial - np.floor(y_dimension/2))

    #Find y standard deviation
    y_std = 0
    for i in range(len(y_rows)):
        y_std += ((((y_bar - i)**2)*y_rows[i]) / (y_summ_denom * (y_summ_denom - 1)))
    final_y_std = math.sqrt(y_std)

    centroid = [x_bar_fit, y_bar_fit]
    uncertainty = [final_x_std, final_y_std]

    return centroid, uncertainty
print(find_centroid_of_image("aptest.FIT"))

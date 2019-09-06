# Photometry that does not account for fractions of a pixel inside/outside of aperature.
# It insteads either accepts them or rejects those borderline pixels.
from astropy.io import fits
import numpy as np
from math import sqrt, pi, log10


image = fits.getdata("aptest.fit")
x_ = int(input("Enter x: "))
y_ = int(input("Enter y: "))
acceptance = input("Accept border pixels? ")

def findPixelCount(x, y, radius):

    dim_change_x = radius
    dim_change_y = radius
    sliced_image = image[(y - dim_change_y):(y + dim_change_y + 1), (x - dim_change_x):(x + dim_change_x + 1)]

    r = radius*2 +1
    star_plus_sky = 0
    counter = 0
    #Determine which pixels are only partially in the aperture
    for i in range(r):
        for j in range(r):
            if acceptance == "y":
                if (sqrt((i - 0.5 - radius)**2 + (j + 0.5 - radius)**2) <= radius) or (sqrt((i - 0.5 - radius)**2 + (j - 0.5 - radius)**2) <= radius) or (sqrt((i + 0.5 - radius)**2 + (j - 0.5 - radius)**2) <= radius) or (sqrt((i + 0.5 - radius)**2 + (j + 0.5 - radius)**2) <= radius):
                    star_plus_sky += sliced_image[i,j]
                    counter += 1
            else:
                if (sqrt((i - 0.5 - radius)**2 + (j + 0.5 - radius)**2) <= radius) and (sqrt((i - 0.5 - radius)**2 + (j - 0.5 - radius)**2) <= radius) and (sqrt((i + 0.5 - radius)**2 + (j - 0.5 - radius)**2) <= radius) and (sqrt((i + 0.5 - radius)**2 + (j + 0.5 - radius)**2) <= radius):
                    star_plus_sky += sliced_image[i,j]
                    counter += 1
    return star_plus_sky, counter


star_sky = findPixelCount(x_, y_, 5)[0]
print(star_sky)
ann_in =  findPixelCount(x_, y_, 8)[0]
ann_out = findPixelCount(x_, y_, 13)[0]
ann = ann_out - ann_in
n_aperture = findPixelCount(x_, y_, 5)[1]
print(n_aperture)
n_ann_out = findPixelCount(x_, y_, 13)[1]
n_ann_in = findPixelCount(x_, y_, 8)[1]
n_ann = n_ann_out - n_ann_in
avg_sky = ann / n_ann
signal = star_sky - (avg_sky * n_aperture)


def findSNR(S, n_ap, n_an, sky, dark, gain, read):
    S = S*gain
    sky = sky * gain
    numerator = sqrt(S)
    denom_term_1 = 1
    denom_term_2 = n_ap*(1+(n_ap/n_an))
    p_squared = read**2 + (gain / sqrt(12.))**2
    denom_term_3 = (sky + dark + p_squared) / S
    denominator = sqrt(denom_term_1 + denom_term_2 * denom_term_3)
    SNR = numerator / denominator
    return SNR

SNR = findSNR(16811.5, n_aperture, n_ann, avg_sky, 10, 2, 15)
print("Signal = {0}".format(signal))
print("SNR = {0}".format(SNR))
print("Noise error margin = {0}".format(signal / SNR))
print("inst mag = {0}".format(-2.5*log10(signal)))
mag_uncertainty = 1.0857 / SNR
print("Mag uncertainty = {0}".format(mag_uncertainty))

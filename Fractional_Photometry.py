# In a given fits image, this program finds finds the apparent magnitude of the asteroid.
# Reference stars must be available, with their image positions and magnitudes known and given.
from astropy.io import fits
import numpy as np
from math import sqrt, pi, log10
import matplotlib.pyplot as plt
import multiprocessing as mp

num_cores = mp.cpu_count()
pool = mp.Pool(num_cores)
image = fits.getdata("1929 SH.00000002.ENTERED COORDINATES.REDUCED.FIT")


x_coords =  [842,844,906,842,601,580]



y_coords =   [503,617,656,816,743,765]



actual_mags =   [14.448,14.457,12.571,12.932,11.510,13.620]



ast_x = 617
ast_y = 591
rad = 8 # Aperture size
rad_ann_in = 10 # Inner annulus size
rad_ann_out = 13 # Outer annulus size


def findPixelCount(x, y, radius):
    dim_change_x = radius
    dim_change_y = radius
    sliced_image = image[(y - dim_change_y):(y + dim_change_y + 1), (x - dim_change_x):(x + dim_change_x + 1)]
    r = radius*2 + 1
    # Determine which pixels are only partially in the aperture
    tot_pix_intensity = 0.0
    for i in range(r):
        for j in range(r):
            frac = 0
            for a in np.arange((i),(i+1.),(0.01)):
                for b in np.arange((j),(j+1.),(0.01)):
                    if sqrt((a - radius)**2 + (b - radius)**2) <= radius:
                        frac += 1
            fraction = frac / 10000.
            tot_pix_intensity +=  (fraction * sliced_image[i,j])
    pix_ct = tot_pix_intensity
    return pix_ct

signals = []
x_axis = []
inst_mags = []
# For the reference stars
for i in range(len(x_coords)):
    star_sky = findPixelCount(x_coords[i], y_coords[i], rad)
    ann_in = findPixelCount(x_coords[i], y_coords[i], rad_ann_in)
    ann_out = findPixelCount(x_coords[i], y_coords[i], rad_ann_out)
    ann = ann_out - ann_in
    n_aperture = pi * (rad**2)
    n_ann = (pi * (rad_ann_out**2)) - (pi * (rad_ann_in**2))
    avg_sky = ann / n_ann
    signal = star_sky - (avg_sky * n_aperture)
    signals.append(signal)
    x_axis_term = log10(signal)
    x_axis.append(x_axis_term)
    inst_mag = -2.5*log10(signal)
    inst_mags.append(inst_mag)
# For the asteroid specifically
star_sky = findPixelCount(ast_x, ast_y, rad)
ann_in = findPixelCount(ast_x, ast_y, rad_ann_in)
ann_out = findPixelCount(ast_x, ast_y, rad_ann_out)
ann = ann_out - ann_in
n_aperture = pi * (rad**2)
n_ann = (pi * (rad_ann_out**2)) - (pi * (rad_ann_in**2))
avg_sky = ann / n_ann
signal = star_sky - (avg_sky * n_aperture)
inst_mag = -2.5*log10(signal)
print("inst mag for asteroid = {0}".format(inst_mag))

# Find asteroid's inst mag uncertainty
def findSNR(S, n_ap, n_an, sky, dark, read, gain):
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
SNR = findSNR(signal, n_aperture, n_ann, avg_sky, 10, 15, 2)
mag_uncertainty = 1.0857 / SNR
print("Mag uncertainty for asteroid = {0}".format(mag_uncertainty))

# Plotting
(m, b), cov = np.polyfit(x_axis, actual_mags, 1, cov=True)
y_axis = []
for i in x_axis:
    y_axis.append(m*i + b)
app_mag = inst_mag + b
print("Constant = {0}".format(b))
print("Covariance: {0}".format(cov))
print("Asteroid apparent mag = {0}".format(app_mag))

# Chi squared
calculated_star_mags = [] # Calculated, observed magnitudes of ref stars
for i in range(len(inst_mags)):
    calculated_star_mags.append(inst_mags[i] + b)
chi_squared = 0
for i in range(len(inst_mags)):
    chi_squared += ((calculated_star_mags[i] - actual_mags[i]) **2 ) / actual_mags[i]
print("χ2 = {0}".format(chi_squared))
print("App mag uncertainty = {0}".format(sqrt( (mag_uncertainty**2) + (cov[1,1]**2) )))
plt.xlabel("log(Signal)", fontname="Times New Roman")
plt.ylabel("Apparent magnitude", fontname="Times New Roman")
plt.scatter(x_axis, actual_mags, color="#000000")
plt.scatter((inst_mag/-2.5), m*(inst_mag/-2.5) + b, color="C1")
plt.text((inst_mag/-2.5) + 0.1, m*(inst_mag/-2.5) + b, "1627 Ivar", fontname="Times New Roman", color="red")
plt.xticks(fontname="Times New Roman")
plt.yticks(fontname="Times New Roman")
plt.plot(x_axis, y_axis, color="#000000")
plt.title("Observation 1: Photometry - Image 00000047: Apparent Magnitude vs. log(Signal)", fontname="Times New Roman" )
plt.text(3.8,11.5, "slope = {0}\nC = {1}\nχ2 = {2} ".format(m,b, chi_squared), fontname="Times New Roman")
plt.show()

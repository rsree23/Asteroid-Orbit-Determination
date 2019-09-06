'''
This is a Least Squares Plate Reduction (LSPR) program. It does "astrometry."
Basically, given the (x,y) positions of some reference stars (in the same image), and their right ascensions
and declinations, the program will output the right ascension and declination of the asteroid (the unknown)
if provided the centroided (x,y) position of the asteroid on the fits image.
Need to find the six plate constants: b1, b2, a11, a12, a21, a22
Need to read x centroid, y centroid, alpha, and delta from file for each of the 12 reference stars
Store each in their own lists (will be four lists)
Access those lists for elements in the matrices
'''
import numpy as np
from math import sqrt, sin, cos, tan, asin, acos, atan, radians, degrees

# Get numbers from file
LSPRfile = open("Obs2.txt", "r")
file_list = LSPRfile.read()
file_list = file_list.split()
value_list = []
for value in file_list:
    value_list.append(eval(value))
value_list = np.asarray(value_list)
# Find each of the four lists from the one file list
x_list = []
y_list = []
alpha_list = []
delta_list = []
k = 0
N = 12
for i in range(len(value_list)):
    if i != 0 and i % 4 ==0:
        k += 1
    if i == 4 * k:
        x_list.append(value_list[i])
    elif i == 1 + 4*k:
        y_list.append(value_list[i])
    elif i == 2 + 4*k:
        alpha_list.append(value_list[i])
    elif i == 3 + 4*k:
        delta_list.append(value_list[i])

# Create the sum alpha 3x1 matrix
alpha_sum = 0
alpha_x_sum = 0
alpha_y_sum = 0
for a in range(len(alpha_list)):
    alpha_sum += alpha_list[a]
    alpha_x_sum += alpha_list[a] * x_list[a]
    alpha_y_sum += alpha_list[a] * y_list[a]
alpha_sum_mat = np.array([ [alpha_sum], [alpha_x_sum], [alpha_y_sum] ])

# Create the sum delta 3x1 matrix
delta_sum = 0
delta_x_sum = 0
delta_y_sum = 0
for d in range(len(delta_list)):
    delta_sum += delta_list[d]
    delta_x_sum += delta_list[d] * x_list[d]
    delta_y_sum += delta_list[d] * y_list[d]
delta_sum_mat = np.array([ [delta_sum], [delta_x_sum], [delta_y_sum] ])

# Create the A matrix and its inverse, both 3x3
x_sum = 0
y_sum = 0
x_sq_sum = 0
y_sq_sum = 0
xy_sum = 0
for i in range(len(x_list)):
    x_sum += x_list[i]
    y_sum += y_list[i]
    x_sq_sum += x_list[i] ** 2
    y_sq_sum += y_list[i] ** 2
    xy_sum += x_list[i] * y_list[i]
A_matrix = np.array([ [N, x_sum, y_sum], [x_sum, x_sq_sum, xy_sum], [y_sum, xy_sum, y_sq_sum] ])
A_inverse = np.linalg.inv(A_matrix)


# Find b1, a11, a12
alpha_product = np.matmul(A_inverse, alpha_sum_mat)
b1 = alpha_product[0,0]
a11 = alpha_product[1,0]
a12 = alpha_product[2,0]
# Find b2, a21, a22
delta_product = np.matmul(A_inverse, delta_sum_mat)
b2 = delta_product[0,0]
a21 = delta_product[1,0]
a22 = delta_product[2,0]

# Functions to find fit alpha and fit delta for given coordinates
def findAlpha(x, y):
    alph = b1 + a11 * x + a12 * y
    return alph
def findDelta(x, y):
    delt = b2 + a21 * x + a22 * y
    return delt


# User inputs (x,y) centroid for their desired object
x_unknown = float(input("Enter x centroid of unknown object: "))
y_unknown = float(input("Enter y centroid of unknown object: "))

# Find alpha and delta
alpha_deg = b1 + a11 * x_unknown + a12 * y_unknown
delta = b2 + a21 * x_unknown + a22 * y_unknown



# Convert alpha to hh:mm:ss.ss
def hourConversion(des):
    des_hrs = int(des)
    des_minutes = int(60 * (des - des_hrs))
    des_seconds_1 = (des_hrs + (des_minutes / 60))
    des_seconds = 3600 * (des - des_seconds_1)
    des_true = "{0}:{1}:{2}".format(des_hrs, des_minutes, round(des_seconds, 4))
    return des_true
alpha_true = hourConversion(alpha_deg)

# Convert delta to hh:mm:ss.ss
delta_true = hourConversion(delta)

# Find alpha uncertainty
chi_a = 0
for i in range(len(alpha_list)):
    chi_a += (alpha_list[i] - findAlpha(x_list[i], y_list[i])) ** 2
alpha_uncertainty = 3600 * sqrt(chi_a / (N - 3))

# Find delta uncertainty
chi_d = 0
for i in range(len(delta_list)):
    chi_d += (delta_list[i] - findDelta(x_list[i], y_list[i])) ** 2
delta_uncertainty = 3600 * sqrt(chi_d / (N - 3))


# Print the results
print("Plate constants: ")
print("b1 = {0} deg".format(b1))
print("b2 = {0} deg".format(b2))
print("a11 = {0} deg/pix".format(a11))
print("a12 = {0} deg/pix".format(a12))
print("a21 = {0} deg/pix".format(a21))
print("a22 = {0} deg/pix".format(a22))

print("Uncertainty: ")
print("RA: {0} arcsec".format(round(alpha_uncertainty, 2)))
print("Dec: {0} arcsec".format(round(delta_uncertainty, 2)))

print("Astrometry for (x,y) = (", x_unknown, y_unknown, ")")
print("RA = {0}".format(alpha_true))
print("Dec = {0}".format(delta_true))





# Optional part
# Convert to radians!!!
# Xi = ξ
# Eta = η

doflats = input("Do field flattening? ")

if doflats == "yes" or "y" or "Yes":



    def findSumOfList(a):
        summ = 0
        for i in a:
            summ += i
        return summ

    summ_alpha = radians(findSumOfList(alpha_list))
    A = summ_alpha / N
    summ_delta = radians(findSumOfList(delta_list))
    D = summ_delta / N

    L = (3911e-3) / (24e-6)

    H = 0
    xisum = 0
    etasum = 0
    xisum_x = 0
    xisum_y = 0
    etasum_x = 0
    etasum_y = 0
    for i in range(N):
        H = sin(radians(delta_list[i])) * sin(D) + cos(radians(delta_list[i])) * cos(D) * cos(radians(alpha_list[i]) - A)
        xi = ((cos(radians(delta_list[i])) * sin(radians(alpha_list[i]) - A)) / H) - (x_list[i] / L)
        eta = ((sin(radians(delta_list[i])) * cos(D) - cos(radians(delta_list[i])) * sin(D) * cos(radians(alpha_list[i]) - A)) / H) - (y_list[i] / L)
        xisum += xi
        etasum += eta
        xisum_x += xi * x_list[i]
        etasum_x += eta * x_list[i]
        xisum_y += xi * y_list[i]
        etasum_y += eta * y_list[i]



    xi_sum_mat = np.array([[xisum], [xisum_x], [xisum_y]])
    eta_sum_mat = np.array([[etasum], [etasum_x], [etasum_y]])

    # Find b1p, a11p, a12p
    xi_product = np.matmul(A_inverse, xi_sum_mat)
    b1p = xi_product[0,0]
    a11p = xi_product[1,0]
    a12p = xi_product[2,0]
    # Find b2p, a21p, a22p
    eta_product = np.matmul(A_inverse, eta_sum_mat)
    b2p = eta_product[0,0]
    a21p = eta_product[1,0]
    a22p = eta_product[2,0]



    def findXi(x, y):
        the_xi = b1p + a11p * x + a12p * y + (x / L)
        the_eta = b2p + a21p * x + a22p * y + (y / L)

        triangle = cos(D) - the_eta * sin(D)
        gamma = sqrt((triangle ** 2) + (the_xi ** 2))

        tanalpha = the_xi / triangle
        final_true_alpha = degrees(atan(tanalpha) + A)
        return final_true_alpha


    def findEta(x, y):
        the_xi = b1p + a11p * x + a12p * y + (x / L)
        the_eta = b2p + a21p * x + a22p * y + (y / L)

        triangle = cos(D) - the_eta * sin(D)
        gamma = sqrt((triangle ** 2) + (the_xi ** 2))

        tandelta = (sin(D) + the_eta * cos(D)) / gamma
        final_true_delta = degrees(atan(tandelta))
        return final_true_delta


    chia = 0
    for i in range((len(alpha_list))):
        chia += ((alpha_list[i] - (findXi(x_list[i], y_list[i]))) ** 2)
    alph = 3600 * sqrt(chia / (N-3))
    print(chia)
    chid = 0
    for i in range(N):
        chid += ((delta_list[i] - (findEta(x_list[i], y_list[i]))) ** 2)
    delt = 3600 * sqrt(chid / (N-3))
    print(chid)


    real_alpha = hourConversion(findXi(x_unknown, y_unknown))
    real_delta = hourConversion(findEta(x_unknown, y_unknown))


    print("-------------------------")
    print("Flat corrections...")
    print("Plate constants: ")
    print("b1 = {0}".format(b1p))
    print("b2 = {0}".format(b2p))
    print("a11 = {0}".format(a11p))
    print("a12 = {0}".format(a12p))
    print("a21 = {0}".format(a21p))
    print("a22 = {0}".format(a22p))

    print("Uncertainty: ")
    print("RA: {0} arcsec".format(round(alph, 2)))
    print("Dec: {0} arcsec".format(round(delt, 2)))

    print("Astrometry for (x,y) = (", x_unknown, y_unknown, ")")
    print("RA = {0}".format(real_alpha))
    print("Dec = {0}".format(real_delta))
else:
    print(" ")

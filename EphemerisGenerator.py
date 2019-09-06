# This is, in a way, a reverse of Orbit_Determination.py:
# Given six orbital elements, it outputs the right ascension and declination of the asteroid.
# Note, M will need to be changed into E (see commented-out funtions below)
import numpy as np
from math import radians, sin, cos, tan, acos, asin, atan, atan2, degrees, sqrt, pi
from functions_OD import *
a = 1.8534126044942794
e = 0.39378153749170036
I = radians(8.406512659898615)
O = radians(133.17295119651467)
w = radians(167.76478487779573)
E = 5.877758741500021
##M_naught = radians(353.5214322928663) #For July 22, 2018
##t_naught = convertToJulianDays(2018, 7, 22, 6, 0, 0) #July 22, 2018
##t_final = convertToJulianDays(2018, 6, 23, 5, 40, 17.399) #Date of observation
#Calculate mean angular velocity, n
k = 0.01720209895 #in AU^1/2 per day, k = sqrt(mu)
n = k / sqrt(a ** 3)
epsilon = radians(23.4352)
#Function finds the mean anomaly, M, given an M initial, epoch, new time, and mean motion value
##def findMeanAnomaly(m_initial, t_initial, new_t, mean_motion): #times in Julian days
##	mean_anomaly = m_initial - mean_motion * (t_initial - t_final)
##	return mean_anomaly
##M = findMeanAnomaly(M_naught, t_naught, t_final, n)
##
##
###Function to find eccentric anomaly, mean_a = M, M = E - e*sin(E)
##def findE(M):
##	E = M
##	calculated_m = (E - e * sin(E))
##	while abs(M - calculated_m) > 1e-6:
##		calculated_m = (E - e * sin(E))
##		fx = M - E + e*sin(E)
##		fxprime = e*cos(E) - 1
##		E = E - (fx/fxprime)
##	return E
##E = findE(M)

#Find physics coordinates x, y, z and put into a 3x1 matrix
x = a * (cos(E) - e)
y = a * sqrt(1- (e **2)) * sin(E)
z = 0
physics = np.array([[x], [y], [z]])

#Initialize transform matrices... will be by the negative of the angles to go from physics to ecliptic, so regular transform matrices
# -w transform about Z axis
w_matrix = np.array([[cos(w), -sin(w), 0], [sin(w), cos(w), 0], [0,0,1]])
# -I transform about that x axis
I_matrix = np.array([[1,0,0], [0,cos(I), -sin(I)], [0, sin(I), cos(I)]])
# -O transform about that z axis
O_matrix = np.array([[cos(O), -sin(O), 0], [sin(O), cos(O), 0], [0,0,1]])

first_transform = np.matmul(w_matrix, physics)
second_transform = np.matmul(I_matrix, first_transform)
third_transform = np.matmul(O_matrix, second_transform)

#Ecliptic vector, r, from sun to asteroid
ecliptic = third_transform


#Earth to sun vector, R, from JPL
big_R = np.array([[-2.626483184093047E-02], [9.322034083425700E-01], [4.040774355056518E-01]])
#Range vector, from earth to asteroid
rho_ecl = big_R + ecliptic

#Convert rho from ecliptic to equatorial, x rotation by -epsilon
equ_rotation = np.array([[1, 0, 0], [0, cos(epsilon), -sin(epsilon)], [0, sin(epsilon), cos(epsilon)]])
rho_equ = np.matmul(equ_rotation, rho_ecl)

#Find alpha and delta from rho in equatorial coordinates
def magnitude(vec):
	unsquared = 0
	for x in range(3):
		unsquared += vec[x, 0] ** 2
	mag = sqrt(unsquared)
	return mag


rho_equ_mag = magnitude(rho_equ)
def findDelta(vec):
	sin_d = sin(vec[2,0] / rho_equ_mag)
	d = asin(sin_d)
	return d
delta = degrees(findDelta(rho_equ))

def findAlpha(vec):
        cos_a = vec[0,0] / (rho_equ_mag * cos(delta))
        sin_a = vec[1,0] / (rho_equ_mag * cos(delta))
        a = atan2(sin_a, cos_a)
        if a < 0:
                a = a + 2*pi
        return a
alpha = degrees(findAlpha(rho_equ))
alpha_hrs = alpha / 15.

#Format the decimals in hh:mm:ss
def formatting(des):
    des_hrs = int(des)
    des_minutes = int(60 * (des - des_hrs))
    des_seconds_1 = (des_hrs + (des_minutes / 60))
    des_seconds = 3600 * (des - des_seconds_1)
    des_true = "{0}:{1}:{2}".format(des_hrs, abs(des_minutes), abs(round(des_seconds, 2)))
    return des_true

print("RA: {0}".format(alpha_hrs))
print("Dec: {0}".format(delta))

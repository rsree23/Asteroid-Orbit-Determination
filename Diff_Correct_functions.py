# This file has the differential correction functions that are called in Orbit_Determination.py
import numpy as np
from math import sqrt, pi, degrees, radians, sin, cos, tan, asin, acos, atan, atan2
from functions_OD import *
k = 0.01720209895
mu = 1.
epsilon = radians(23.4352)

# Finds the magnitude of a vector (analogous to np.linalg.norm)
def magnitude(vec):
    unsquared = 0
    for x in range(3):
        unsquared += (vec[x, 0] ** 2)
    mag = sqrt(unsquared)
    return mag

# Format the decimals in hh:mm:ss
def formatting(des):
    des_hrs = int(des)
    des_minutes = int(60 * (des - des_hrs))
    des_seconds_1 = (des_hrs + (des_minutes / 60))
    des_seconds = 3600 * (des - des_seconds_1)
    des_true = "{0}:{1}:{2}".format(des_hrs, abs(des_minutes), abs(round(des_seconds, 2)))
    return des_true


# Functions below used for differential correction
def findAlphaAndDelt(r2_x, r2_y, r2_z, r2_dot_x, r2_dot_y, r2_dot_z, t_naught, tee, t_want, bigRx, bigRy, bigRz):
    r2 = np.array([[r2_x], [r2_y], [r2_z]])
    r2_dot = np.array([[r2_dot_x], [r2_dot_y], [r2_dot_z]])

    #Convert from equatorial to ecliptic
    ecl_transform = np.array([[1, 0, 0], [0, cos(epsilon), sin(epsilon)], [0, -sin(epsilon), cos(epsilon)]])
    r2 = np.matmul(ecl_transform, r2)
    r2_dot = np.matmul(ecl_transform, r2_dot)
    r2 = np.array([r2[0,0], r2[1,0], r2[2,0]])
    r2_dot = np.array([r2_dot[0,0], r2_dot[1,0], r2_dot[2,0]])

    #Finds the magnitude of the vectors
    r2_mag = np.linalg.norm(r2)
    r2_dot_mag = np.linalg.norm(r2_dot)

    #Find a, the semi-major axis
    a = ((2. / r2_mag) - (r2_dot_mag ** 2.)) ** -1.

    #Find e, the eccentricity
    #h = cross product of r and rdot, angular momentum per unit mass
    h = np.cross(r2, r2_dot)
    h_mag = np.linalg.norm(h)
    e = sqrt(1. - ((h_mag) ** 2. / a))

    #Find I (Inclination) and O (Capital Omega, longitude of the ascending node)
    cos_I = h[2] / h_mag
    I = angleCheck(acos(cos_I))
    sin_O = h[0] / (sin(I) * h_mag)
    cos_O = -h[1] / (sin(I) * h_mag)
    O = angleCheck(atan2(sin_O, cos_O))

    #Find w (argument of the perehelion)
    #First find w + f, then f, then w
    sin_wf = r2[2] / (r2_mag * sin(I))
    cos_wf = ((r2[0] / r2_mag) + cos(I) * sin_wf * sin(O)) / cos(O)
    w_plus_f = angleCheck(atan2(sin_wf, cos_wf))
    cos_f =  (((a * (1. - e ** 2.)) / r2_mag) - 1.) / e
    sin_f = (np.dot(r2, r2_dot) / (e * r2_mag)) * sqrt(a * (1.-e**2.))
    f = angleCheck(atan2(sin_f, cos_f))
    w = angleCheck(w_plus_f - f)

    #Find E2, then M2, then M0 (mean anomoly)
    cos_E2 = (1. - (r2_mag / a)) / e
    E2 = acos(cos_E2)
    #E2 and f need to be in the same half plane / quadrant
    if f > pi and f < 2. * pi:
        E2 = 2. * pi - E2
    M2 = E2 - e * sin(E2) #Finds mean anomoly for this observation night
    n = k / (a**1.5)
    M0 = M2 + n*(t_naught-tee) #Finds mean anomoly on July 22, 2018 at 6:00 UT

########Ephem generation##############
    M_naught = M0

    #Calculate mean angular velocity, n
    n = k / sqrt(a ** 3)

    #Function finds the mean anomaly, M, given an M initial, epoch, new time, and mean motion value
    M = n * (t_want - t_naught) + M_naught


    #Function to find eccentric anomaly, mean_a = M, M = E - e*sin(E)
    E = M
    calculated_m = (E - e * sin(E))
    while abs(M - calculated_m) > 1e-6:
        calculated_m = (E - e * sin(E))
        fx = M - E + e*sin(E)
        fxprime = e*cos(E) - 1
        E = E - (fx/fxprime)

    #Find physics coordinates x, y, z and put into a 3x1 matrix
    x = a * (cos(E) - e)
    y = a * sqrt(1- (e **2)) * sin(E)
    z = 0
    physics = np.array([[x], [y], [z]])

    #Initialize transform matrices... will be by the negative of the angles to go from physics to ecliptic, so regular transform matrices
    # -w transform about Z axis
    w_matrix = np.array([[cos(w), -sin(w), 0], [sin(w), cos(w), 0], [0,0,1]])
    I_matrix = np.array([[1,0,0], [0,cos(I), -sin(I)], [0, sin(I), cos(I)]])
    O_matrix = np.array([[cos(O), -sin(O), 0], [sin(O), cos(O), 0], [0,0,1]])

    first_transform = np.matmul(w_matrix, physics)
    second_transform = np.matmul(I_matrix, first_transform)
    third_transform = np.matmul(O_matrix, second_transform)

    #ECLIPTIC vector, r, from sun to asteroid
    ecliptic = third_transform
    ecliptic = np.array([ecliptic[0,0], ecliptic[1,0], ecliptic[2,0]])

    #Earth to sun vector, R, from JPL, ECLIPTIC
    big_R = np.array([bigRx, bigRy, bigRz])
    #Range vector, from earth to asteroid
    rho_ecl = big_R + ecliptic

    rho_ecl = np.array([ [rho_ecl[0]], [rho_ecl[1]], [rho_ecl[2]]  ])
    #Convert rho from ecliptic to equatorial, x rotation by -epsilon
    equ_rotation = np.array([[1, 0, 0], [0, cos(epsilon), -sin(epsilon)], [0, sin(epsilon), cos(epsilon)]])
    rho_equ = np.matmul(equ_rotation, rho_ecl)

    #Find alpha and delta from rho in equatorial coordinates
    rho_equ_mag = magnitude(rho_equ)

    sin_d = sin(rho_equ[2,0] / rho_equ_mag)
    d = asin(sin_d)

    delta = d

    cos_a = rho_equ[0,0] / (rho_equ_mag * cos(delta))
    sin_a = rho_equ[1,0] / (rho_equ_mag * cos(delta))
    a = atan2(sin_a, cos_a)
    if a < 0:
        a = a + 2*pi
    alpha = a

    return alpha, delta

def differentiallyCorrect(r, rdot, t_nought, time_index, Rx, Ry, Rz, ra_hours, ra_mins, ra_secs, dec_degrees, dec_mins, dec_secs):
    num_obs = 5 * 2 ##Number of observation nights times two to account for alpha and delta
    time_list = [convertToJulianDays(2018,6,23,5,40,17.399), convertToJulianDays(2018,6,28,3,46,3.530), convertToJulianDays(2018,7,2,3,33,35.710), convertToJulianDays(2018,7,10,6,0,11.299), convertToJulianDays(2018,7,18,5,58,8.390)]
    change = 1e-4
    te = time_list[time_index]
    #Transform big R from equatorial to ecliptic
    Rx1 = np.array([[Rx[0]], [Ry[0]], [Rz[0]]])
    Rx2 = np.array([[Rx[1]], [Ry[1]], [Rz[1]]])
    Rx3 = np.array([[Rx[2]], [Ry[2]], [Rz[2]]])
    Rx4 = np.array([[Rx[3]], [Ry[3]], [Rz[3]]])
    Rx5 = np.array([[Rx[4]], [Ry[4]], [Rz[4]]])
    ecl_transform = np.array([[1, 0, 0], [0, cos(epsilon), sin(epsilon)], [0, -sin(epsilon), cos(epsilon)]])
    Rx1 = np.matmul(ecl_transform, Rx1)
    Rx2 = np.matmul(ecl_transform, Rx2)
    Rx3 = np.matmul(ecl_transform, Rx3)
    Rx4 = np.matmul(ecl_transform, Rx4)
    Rx5 = np.matmul(ecl_transform, Rx5)
    Rxs = np.array([float(Rx1[0]), float(Rx2[0]), float(Rx3[0]), float(Rx4[0]), float(Rx5[0])])
    Rys = np.array([float(Rx1[1]), float(Rx2[1]), float(Rx3[1]), float(Rx4[1]), float(Rx5[1])])
    Rzs = np.array([float(Rx1[2]), float(Rx2[2]), float(Rx3[2]), float(Rx4[2]), float(Rx5[2])])


    #Observed RA and Dec - Calculated RA and Dec (the triangle alpha term)
    obs_alpha_list = []
    obs_delta_list = []
    for i in range(5):
        obs_alpha_list.append(convertRAToRadians(ra_hours[i], ra_mins[i], ra_secs[i]))
        obs_delta_list.append(convertDecToRadians(dec_degrees[i], dec_mins[i], dec_secs[i]))
    observed_matrix = np.array([  [obs_alpha_list[0], obs_delta_list[0]], [obs_alpha_list[1], obs_delta_list[1]], [obs_alpha_list[2], obs_delta_list[2]], [obs_alpha_list[3], obs_delta_list[3]], [obs_alpha_list[4], obs_delta_list[4]]   ])

    calc_alpha_list = []
    calc_delta_list = []
    for i in range(5):
        calc_alpha_list.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        calc_delta_list.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    calculated_matrix = np.array([[calc_alpha_list[0], calc_delta_list[0]], [calc_alpha_list[1], calc_delta_list[1]], [calc_alpha_list[2], calc_delta_list[2]], [calc_alpha_list[3], calc_delta_list[3]], [calc_alpha_list[4], calc_delta_list[4]]   ])


    triangle_alpha = np.sum(observed_matrix - calculated_matrix)
    RMS_initial = sqrt(np.sum((observed_matrix - calculated_matrix) **2 ) / (num_obs - 6))

    #Find the partial derivatives as matrices
    #Find the elements for the matrices
    alphs_xpos = []
    alphs_ypos = []
    alphs_zpos = []
    alphs_xdotpos = []
    alphs_ydotpos = []
    alphs_zdotpos = []
    delts_xpos = []
    delts_ypos = []
    delts_zpos = []
    delts_xdotpos = []
    delts_ydotpos = []
    delts_zdotpos = []

    alphs_xneg = []
    alphs_yneg = []
    alphs_zneg = []
    alphs_xdotneg = []
    alphs_ydotneg = []
    alphs_zdotneg = []
    delts_xneg = []
    delts_yneg = []
    delts_zneg = []
    delts_xdotneg = []
    delts_ydotneg = []
    delts_zdotneg = []


    for i in range(5):
        alphs_xpos.append(findAlphaAndDelt(r[0] + change, r[1], r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_xpos.append(findAlphaAndDelt(r[0] + change, r[1], r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    for i in range(5):
        alphs_xneg.append(findAlphaAndDelt(r[0] - change, r[1], r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_xneg.append(findAlphaAndDelt(r[0] - change, r[1], r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    for i in range(5):
        alphs_ypos.append(findAlphaAndDelt(r[0], r[1] + change, r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_ypos.append(findAlphaAndDelt(r[0], r[1] + change, r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    for i in range(5):
        alphs_yneg.append(findAlphaAndDelt(r[0], r[1] - change, r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_yneg.append(findAlphaAndDelt(r[0], r[1] - change, r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    for i in range(5):
        alphs_zpos.append(findAlphaAndDelt(r[0], r[1], r[2] + change, rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_zpos.append(findAlphaAndDelt(r[0], r[1], r[2] + change, rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    for i in range(5):
        alphs_zneg.append(findAlphaAndDelt(r[0], r[1], r[2] - change, rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_zneg.append(findAlphaAndDelt(r[0], r[1], r[2] - change, rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])



    for i in range(5):
        alphs_xdotpos.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0] + change, rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_xdotpos.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0] + change, rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    for i in range(5):
        alphs_xdotneg.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0] - change, rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_xdotneg.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0] - change, rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    for i in range(5):
        alphs_ydotpos.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1] + change, rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_ydotpos.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1] + change, rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    for i in range(5):
        alphs_ydotneg.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1] - change, rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_ydotneg.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1] - change, rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    for i in range(5):
        alphs_zdotpos.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1], rdot[2] + change, t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_zdotpos.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1], rdot[2] + change, t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    for i in range(5):
        alphs_zdotneg.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1], rdot[2] - change, t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        delts_zdotneg.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1], rdot[2] - change, t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    #Find d/dx
    x_plus_triangle = np.array([  [alphs_xpos[0], delts_xpos[0]], [alphs_xpos[1], delts_xpos[1]], [alphs_xpos[2], delts_xpos[2]], [alphs_xpos[3], delts_xpos[3]], [alphs_xpos[4], delts_xpos[4]]   ])
    x_minus_triangle = np.array([  [alphs_xneg[0], delts_xneg[0]], [alphs_xneg[1], delts_xneg[1]], [alphs_xneg[2], delts_xneg[2]], [alphs_xneg[3], delts_xneg[3]], [alphs_xneg[4], delts_xneg[4]]   ])
    d_dx = np.sum((x_plus_triangle - x_minus_triangle) / (2*change))

    #Find d/dy
    y_plus_triangle = np.array([  [alphs_ypos[0], delts_ypos[0]], [alphs_ypos[1], delts_ypos[1]], [alphs_ypos[2], delts_ypos[2]], [alphs_ypos[3], delts_ypos[3]], [alphs_ypos[4], delts_ypos[4]]   ])
    y_minus_triangle = np.array([  [alphs_yneg[0], delts_yneg[0]], [alphs_yneg[1], delts_yneg[1]], [alphs_yneg[2], delts_yneg[2]], [alphs_yneg[3], delts_yneg[3]], [alphs_yneg[4], delts_yneg[4]]   ])
    d_dy = np.sum((y_plus_triangle - y_minus_triangle) / (2*change))

    #Find d/dz
    z_plus_triangle = np.array([  [alphs_zpos[0], delts_zpos[0]], [alphs_zpos[1], delts_zpos[1]], [alphs_zpos[2], delts_zpos[2]], [alphs_zpos[3], delts_zpos[3]], [alphs_zpos[4], delts_zpos[4]]   ])
    z_minus_triangle = np.array([  [alphs_zneg[0], delts_zneg[0]], [alphs_zneg[1], delts_zneg[1]], [alphs_zneg[2], delts_zneg[2]], [alphs_zneg[3], delts_zneg[3]], [alphs_zneg[4], delts_zneg[4]]   ])
    d_dz = np.sum((z_plus_triangle - z_minus_triangle) / (2*change))


    #Find d/dxdot
    xdot_plus_triangle = np.array([  [alphs_xdotpos[0], delts_xdotpos[0]], [alphs_xdotpos[1], delts_xdotpos[1]], [alphs_xdotpos[2], delts_xdotpos[2]], [alphs_xdotpos[3], delts_xdotpos[3]], [alphs_xdotpos[4], delts_xdotpos[4]]   ])
    xdot_minus_triangle = np.array([  [alphs_xdotneg[0], delts_xdotneg[0]], [alphs_xdotneg[1], delts_xdotneg[1]], [alphs_xdotneg[2], delts_xdotneg[2]], [alphs_xdotneg[3], delts_xdotneg[3]], [alphs_xdotneg[4], delts_xdotneg[4]]   ])
    d_dxdot = np.sum((xdot_plus_triangle - xdot_minus_triangle) / (2*change))

    #Find d/dydot
    ydot_plus_triangle = np.array([  [alphs_ydotpos[0], delts_ydotpos[0]], [alphs_ydotpos[1], delts_ydotpos[1]], [alphs_ydotpos[2], delts_ydotpos[2]], [alphs_ydotpos[3], delts_ydotpos[3]], [alphs_ydotpos[4], delts_ydotpos[4]]   ])
    ydot_minus_triangle = np.array([  [alphs_ydotneg[0], delts_ydotneg[0]], [alphs_ydotneg[1], delts_ydotneg[1]], [alphs_ydotneg[2], delts_ydotneg[2]], [alphs_ydotneg[3], delts_ydotneg[3]], [alphs_ydotneg[4], delts_ydotneg[4]]   ])
    d_dydot = np.sum((ydot_plus_triangle - ydot_minus_triangle) / (2*change))


    #Find d/dzdot
    zdot_plus_triangle = np.array([  [alphs_zdotpos[0], delts_zdotpos[0]], [alphs_zdotpos[1], delts_zdotpos[1]], [alphs_zdotpos[2], delts_zdotpos[2]], [alphs_zdotpos[3], delts_zdotpos[3]], [alphs_zdotpos[4], delts_zdotpos[4]]   ])
    zdot_minus_triangle = np.array([  [alphs_zdotneg[0], delts_zdotneg[0]], [alphs_zdotneg[1], delts_zdotneg[1]], [alphs_zdotneg[2], delts_zdotneg[2]], [alphs_zdotneg[3], delts_zdotneg[3]], [alphs_zdotneg[4], delts_zdotneg[4]]   ])
    d_dzdot = np.sum((zdot_plus_triangle - zdot_minus_triangle) / (2*change))

    #Find the a matrix
    a_matrix = np.array([ [triangle_alpha * d_dx], [triangle_alpha * d_dy], [triangle_alpha * d_dz], [triangle_alpha * d_dxdot], [triangle_alpha * d_dydot], [triangle_alpha * d_dzdot] ])

    #Find the J matrix

    J_matrix = np.array([ [d_dx ** 2, d_dx * d_dy, d_dx * d_dz, d_dx * d_dxdot, d_dx * d_dydot, d_dx * d_dzdot],
                          [d_dx * d_dy, d_dy ** 2, d_dy * d_dz, d_dy * d_dxdot, d_dy * d_dydot, d_dy * d_dzdot],
                          [d_dx * d_dz, d_dy * d_dz, d_dz ** 2, d_dz * d_dxdot, d_dz * d_dydot, d_dz * d_dzdot],
                          [d_dx * d_dxdot, d_dy * d_dxdot, d_dz * d_dxdot, d_dxdot ** 2, d_dxdot * d_dydot, d_dxdot * d_dzdot],
                          [d_dx * d_dydot, d_dy * d_dydot, d_dz * d_dydot, d_dxdot * d_dydot, d_dydot ** 2, d_dydot * d_dzdot],
                          [d_dx * d_dzdot, d_dy * d_dzdot, d_dz * d_dzdot, d_dxdot * d_dzdot, d_dydot * d_dzdot, d_dzdot ** 2] ])

    #Find the x matrix
    x_matrix = np.linalg.lstsq(J_matrix, a_matrix, rcond=True) #Used instead of inverse J since determinant of J = 0

    change_x = x_matrix[0][0,0]
    change_y = x_matrix[0][1,0]
    change_z = x_matrix[0][2,0]
    change_xdot = x_matrix[0][3,0]
    change_ydot = x_matrix[0][4,0]
    change_zdot = x_matrix[0][5,0]

    r[0] += change_x
    r[1] += change_y
    r[2] += change_z
    rdot[0] += change_xdot
    rdot[1] += change_ydot
    rdot[2] += change_zdot


    new_orbital_elements = findOrbitalElements(r, rdot, t_nought, te)
    calc_alpha_list = []
    calc_delta_list = []
    for i in range(5):
        calc_alpha_list.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[0])
        calc_delta_list.append(findAlphaAndDelt(r[0], r[1], r[2], rdot[0], rdot[1], rdot[2], t_nought, te, time_list[i], Rxs[i], Rys[i], Rzs[i])[1])
    calculated_matrix = np.array([[calc_alpha_list[0], calc_delta_list[0]], [calc_alpha_list[1], calc_delta_list[1]], [calc_alpha_list[2], calc_delta_list[2]], [calc_alpha_list[3], calc_delta_list[3]], [calc_alpha_list[4], calc_delta_list[4]]   ])


    triangle_alpha = np.sum(observed_matrix - calculated_matrix)
    RMS_final = sqrt(np.sum((observed_matrix - calculated_matrix) **2 ) / (num_obs - 6))



    return new_orbital_elements, RMS_initial, RMS_final, r, rdot

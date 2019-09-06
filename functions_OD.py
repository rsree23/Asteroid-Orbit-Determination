# These functions are called in Orbit_Determination.py.
import numpy as np
from math import sqrt, pi, degrees, radians, sin, cos, tan, asin, acos, atan, atan2

epsilon = radians(23.4352)
k = 0.01720209895

#Converts civil date/time to Julain days
def convertToJulianDays(Y, M, D, H, mins, sec):
    J_nought = 367. * Y - int((7. * (Y + int((M + 9.) / 12.))) / 4.) + int((275. * M) / 9) + D + 1721013.5
    T = (H) + (mins / 60.) + (sec / 3600)
    JD = J_nought + (T / 24.)
    return JD

#Converts hh mm ss right ascension to decimal radians
def convertRAToRadians(RA_hrs, RA_mins, RA_secs):
    RA_in_hours = RA_hrs + (RA_mins / 60.) + (RA_secs / 3600.)
    RA = radians(RA_in_hours * 15)
    return RA

#Converts dd mm ss declination to decimal radians
def convertDecToRadians(dec_hours, dec_mins, dec_secs):
    dec = radians(dec_hours + (dec_mins / 60.) + (dec_secs / 3600.))
    return dec

#Will account for negative angles
def angleCheck(angle):
    if angle < 0:
        angle = angle + 2.*pi
    return angle

#Baby OD
def findOrbitalElements(r2_new, r2_dot_new, t0, t2):
    r2 = r2_new
    r2_dot = r2_dot_new
    r2_x = r2[0]
    r2_y = r2[1]
    r2_z = r2[2]
    r2_dot_x = r2_dot[0]
    r2_dot_y = r2_dot[1]
    r2_dot_z = r2_dot[2]
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
    M0 = M2 + n*(t0-t2) #Finds mean anomoly on July 22, 2018 at 6:00 UT

    return a, e, I, O, w, M0, r2, r2_dot

# This is the main orbit determination program. Given right ascensions and declinations (see input file),
# This will provide the six orbital elements that describe the orbit and location of the asteroid in space.
# This also differentially corrects those elements to, basically, make them more precise.
import numpy as np
from math import sqrt, pi, degrees, radians, sin, cos, tan, asin, acos, atan, atan2
from functions_OD import *
from Diff_Correct_functions import *

# Get numbers from file
ODfile = open("OD_Input_File.txt", "r")
file_list = ODfile.read()
file_list = file_list.split()
value_list = []
for value in file_list:
    value_list.append(eval(value))
value_list = np.asarray(value_list) # Every single value from the input file

mids = [1,2,3] # Middle observation indices
fir = 0 # Index of first observation
las = 4 # Index of last observation

for mid in mids:
    # Find each of the four lists from the one file list
    years = []
    months = []
    days = []
    hours = []
    mins =[]
    secs = []
    ra_hours = []
    ra_mins = []
    ra_secs = []
    dec_degrees = []
    dec_mins = []
    dec_secs = []
    Rx = []
    Ry = []
    Rz = []
    k = 0 # Different k (just a counter), is assigned to proper value below
    for i in range(len(value_list)):
        if i != 0 and i % 15 ==0:
            k += 1
        if i == 15 * k:
            years.append(value_list[i])
        elif i == 1 + 15*k:
            months.append(value_list[i])
        elif i == 2 + 15*k:
            days.append(value_list[i])
        elif i == 3 + 15*k:
            hours.append(value_list[i])
        elif i == 4 + 15*k:
            mins.append(value_list[i])
        elif i == 5 + 15*k:
            secs.append(value_list[i])
        elif i == 6 + 15*k:
            ra_hours.append(value_list[i])
        elif i == 7 + 15*k:
            ra_mins.append(value_list[i])
        elif i == 8 + 15*k:
            ra_secs.append(value_list[i])
        elif i == 9 + 15*k:
            dec_degrees.append(value_list[i])
        elif i == 10 + 15*k:
            dec_mins.append(value_list[i])
        elif i == 11 + 15*k:
            dec_secs.append(value_list[i])
        elif i == 12 + 15*k:
            Rx.append(value_list[i])
        elif i == 13 + 15*k:
            Ry.append(value_list[i])
        elif i == 14 + 15*k:
            Rz.append(value_list[i])

    # Constants used throughout
    k = 0.01720209895 #proper k
    speed_o_light = 173.264498209 #in AU/day
    mu = 1. #Gaussian
    epsilon = radians(23.4352)

    # Times defined (Non-Gaussian and Gaussian)
    t0 = convertToJulianDays(2018, 7, 22, 6, 0, 0)
    t1 = convertToJulianDays(years[fir], months[fir], days[fir], hours[fir], mins[fir], secs[fir])
    t2 = convertToJulianDays(years[mid], months[mid], days[mid], hours[mid], mins[mid], secs[mid])
    t3 = convertToJulianDays(years[las], months[las], days[las], hours[las], mins[las], secs[las])
    tau1 = k * (t1-t2)
    tau3 = k * (t3-t2)
    tau = tau3 - tau1

    # Define Right Ascensions and Declinations from file
    alpha = [convertRAToRadians(ra_hours[fir], ra_mins[fir], ra_secs[fir]), convertRAToRadians(ra_hours[mid], ra_mins[mid], ra_secs[mid]), convertRAToRadians(ra_hours[las], ra_mins[las], ra_secs[las])]
    delta = [convertDecToRadians(dec_degrees[fir], dec_mins[fir], dec_secs[fir]), convertDecToRadians(dec_degrees[mid], dec_mins[mid], dec_secs[mid]), convertDecToRadians(dec_degrees[las], dec_mins[las], dec_secs[las])]
    alpha1 = alpha[0]
    alpha2 = alpha[1]
    alpha3 = alpha[2]
    delta1 = delta[0]
    delta2 = delta[1]
    delta3 = delta[2]

    # Define rho hat unit vectors
    rho_hat_1 = np.array([cos(alpha1)*cos(delta1), sin(alpha1)*cos(delta1), sin(delta1)])
    rho_hat_2 = np.array([cos(alpha2)*cos(delta2), sin(alpha2)*cos(delta2), sin(delta2)])
    rho_hat_3 = np.array([cos(alpha3)*cos(delta3), sin(alpha3)*cos(delta3), sin(delta3)])
    rho_hat = np.array([rho_hat_1, rho_hat_2, rho_hat_3])

    # Define Sun vector (R)
    big_R_1 = np.array([Rx[fir], Ry[fir], Rz[fir]])
    big_R_2 = np.array([Rx[mid], Ry[mid], Rz[mid]])
    big_R_3 = np.array([Rx[las], Ry[las], Rz[las]])
    big_R = np.array([big_R_1, big_R_2, big_R_3])

    # D constants - need for scalar equations of range
    D1 = []
    D2 = []
    D3 = []
    for i in range(3):
        D1i = np.dot(np.cross(big_R[i], rho_hat_2), rho_hat_3)
        D1.append(D1i)
    for i in range(3):
        D2i = np.dot(np.cross(rho_hat_1, big_R[i]), rho_hat_3)
        D2.append(D2i)
    for i in range(3):
        D3i = np.dot(rho_hat_1, np.cross(rho_hat_2, big_R[i]))
        D3.append(D3i)
    # Explicitly defined
    D0 = np.dot(rho_hat[0], np.cross(rho_hat[1], rho_hat[2]))
    D11 = D1[0]
    D12 = D1[1]
    D13 = D1[2]
    D21 = D2[0]
    D22 = D2[1]
    D23 = D2[2]
    D31 = D3[0]
    D32 = D3[1]
    D33 = D3[2]

    # Solve 8th order polynomial to get an initial value for r2
    # Define constants
    A1 = tau3 / tau
    A3 = -tau1 / tau
    B1 = (A1 * (tau**2 - tau3**2)) / 6.
    B3 = (A3 * (tau**2 - tau1**2)) / 6.
    A = (A1 * D21 - D22 + A3 * D23) / -D0
    B = (B1 * D21 + B3 * D23) / -D0
    E = -2. * np.dot(rho_hat_2, big_R_2)
    F = np.linalg.norm(big_R_2) ** 2.
    a = -1. *(A**2. + A*E + F)
    b = -1. * mu * (2. * A * B + B * E)
    c = -1. * (mu**2.) * (B**2.)

    # 8th order polynomial: (r2**8) + (a * r2**6) + (b * r2**3) + (c) = 0
    polynomial_coefficients = [c, 0, 0, b, 0, 0, a, 0, 1]
    polynomial_roots = np.polynomial.polynomial.polyroots(polynomial_coefficients)

    # Set r2_initial equal to the real root(s) from the list above
    r2_mag = 0
    poly_reals = []
    for i in polynomial_roots:
        if np.isreal(i) == True and np.real(i)>0:
            real_part = np.real(i)
            poly_reals.append(real_part)

    for root in poly_reals:
        r2_mag = root

        # Truncated Taylor series to get f and g
        f1 = 1 - (mu / (2 * (r2_mag**3))) * (tau1 ** 2)
        f3 = 1 - (mu / (2 * (r2_mag**3))) * (tau3 ** 2)
        g1 = tau1 - (mu / (6 * (r2_mag**3))) * (tau1 ** 3)
        g3 = tau3 - (mu / (6 * (r2_mag**3))) * (tau3 ** 3)

        # Can now define c values
        c1 = g3 / (f1*g3 - g1*f3)
        c2 = -1
        c3 = -g1 / (f1*g3 - g1*f3)
        c_list = [c1, c2, c3]

        # Find rho scalars
        rho1 = (c1*D11 + c2*D12 + c3*D13) / (c1*D0)
        rho2 = (c1*D21 + c2*D22 + c3*D23) / (c2*D0)
        rho3 = (c1*D31 + c2*D32 + c3*D33) / (c3*D0)

        # r vectors
        r1 = np.asarray(rho1 * rho_hat_1 - big_R_1)
        r2 = np.asarray(rho2 * rho_hat_2 - big_R_2)
        r3 = np.asarray(rho3 * rho_hat_3 - big_R_3)

        # r2 dot vector
        d1 = (-f3 / (f1*g3 - f3*g1))
        d3 = (f1 / (f1*g3 - f3*g1))
        r2_dot = np.asarray( d1* r1 +  d3* r3)

        # All above to get INITIAL values for r2 and r2 dot vectors
        # Now the beginning of the iterations...

        r2_new = r2
        r2_dot_new = r2_dot
        r2_mag = 0

        while abs(r2_mag - np.linalg.norm(r2_new)) > 1e-10:
            # Updating vectors before iteration
            r2_initial = r2_new
            r2_dot_initial = r2_dot_new
            r2_mag = np.linalg.norm(r2_initial)

            # Light travel time corrections before every iteration
            tnew1 = t1 - (rho1 / speed_o_light)
            tnew2 = t2 - (rho2 / speed_o_light)
            tnew3 = t3 - (rho3 / speed_o_light)
            tau1 = k * (tnew1-tnew2)
            tau3 = k * (tnew3-tnew2)

            # Define constants used in 4 term Taylor series...
            u = mu / (r2_mag**3)
            z = np.dot(r2_initial, r2_dot_initial) / (r2_mag**2)
            q = (np.dot(r2_dot_initial, r2_dot_initial) / (r2_mag**2)) - u



            # f and g with extended Taylor series
            f1 = 1 - 0.5*u*(tau1**2) + (1/2)*u*z*(tau1**3) + (1/24)*(3*u*q - 15*u*(z**2) + u**2)*(tau1**4)
            f3 = 1 - 0.5*u*(tau3**2) + (1/2)*u*z*(tau3**3) + (1/24)*(3*u*q - 15*u*(z**2) + u**2)*(tau3**4)
            g1 = tau1 - (1/6)*u*(tau1**3) + 0.25*u*z*(tau1**4)
            g3 = tau3 - (1/6)*u*(tau3**3) + 0.25*u*z*(tau3**4)

            # Can now define c values
            c1 = g3 / (f1*g3 - g1*f3)
            c2 = -1
            c3 = -g1 / (f1*g3 - g1*f3)

            # Find rho scalars
            rho1 = (c1*D11 + c2*D12 + c3*D13) / (c1*D0)
            rho2 = (c1*D21 + c2*D22 + c3*D23) / (c2*D0)
            rho3 = (c1*D31 + c2*D32 + c3*D33) / (c3*D0)

            # r vectors
            r1 = rho1 * rho_hat_1 - big_R_1
            r2 = rho2 * rho_hat_2 - big_R_2
            r3 = rho3 * rho_hat_3 - big_R_3

            r1_new = r1
            r2_new = r2
            r3_new = r3
            # r2 dot vector
            r2_dot_new = (-f3 / (f1*g3 - f3*g1)) * r1_new + (f1 / (f1*g3 - f3*g1)) * r3_new


        # The Orbital Elements!!!!
        orbital_elements = findOrbitalElements(r2_new, r2_dot_new, t0, t2)
        a = orbital_elements[0]
        e = orbital_elements[1]
        I = orbital_elements[2]
        O = orbital_elements[3]
        w = orbital_elements[4]
        M = orbital_elements[5]
        print("Middle Observation = {0} using Scalar Equation of Lagrange root = {1}".format(mid + 1, root))
        print("Position and velocity vectors, in AU and AU/day:")
        print("Position vector: {0}\nVelocity vector: {1}".format(r2_new, r2_dot_new * k)) #Since velocity was in AU/Gaussian day, multiply by k to get back into AU/day
        print("Range to asteroid, in AU: {0}".format(rho2))
        print("Orbital elements")
        print("a = {0} AU".format(a))
        print("e = {0}".format(e))
        print("i = {0} degrees".format(degrees(I)))
        print("Ω = {0} degrees".format(degrees(O)))
        print("ω = {0} degrees".format(degrees(w)))
        print("M on July 22, 2018 6:00 UT = {0} degrees".format(degrees(M)))

        # Differential correction option
        if rho2 > 0:
            user_input = input("Differentially correct orbital elements with Middle Observation = {0} (y/n)?\n".format(mid + 1))
            if user_input == "y":
                correction_elements = differentiallyCorrect(r2_new, r2_dot_new, t0, mid, Rx, Ry, Rz, ra_hours, ra_mins, ra_secs, dec_degrees, dec_mins, dec_secs)
                new_orbital_elements = correction_elements[0]
                a = new_orbital_elements[0]
                e = new_orbital_elements[1]
                I = new_orbital_elements[2]
                O = new_orbital_elements[3]
                w = new_orbital_elements[4]
                M = new_orbital_elements[5]
                firstRMS = correction_elements[1]
                lastRMS = correction_elements[2]
                position_vec = new_orbital_elements[6]
                velocity_vec = new_orbital_elements[7]
                print("Differentially Corrected Orbital elements...")
                print("a = {0} AU".format(a))
                print("e = {0}".format(e))
                print("i = {0} degrees".format(degrees(I)))
                print("Ω = {0} degrees".format(degrees(O)))
                print("ω = {0} degrees".format(degrees(w)))
                print("M on July 22, 2018 6:00 UT = {0} degrees\n".format(degrees(M)))
                print("Initial RMS = {0}".format(firstRMS))
                print("New RMS = {0}".format(lastRMS))
                print("Position = {0}, {1}".format(position_vec, correction_elements[3]))
                print("Velocity = {0}, {1}\n\n\n------------------------------\n\n\n".format(velocity_vec * k, correction_elements[4]))
        else:
            print("Range was negative, these orbital elements are not suitable, cannot be differentially corrected.\n\n\n-----------------------\n\n\n")

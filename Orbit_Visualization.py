# Uses vpython to visualize an orbit in space. Has not been implemented with asteroid's orbital elemets used
# other files.
from vpython import *
from math import cos, sin, sqrt, pi
import numpy as np


a = 2.773017979589484
e = 0.1750074901308245
M = radians(336.0050001501443)
Oprime = radians(108.032597191534)
iprime = radians(16.34548466739393)
wprime = radians(74.95130563682554)

def solvekep(M):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M
    while abs(Mguess - Mtrue) > 1e-004: #absolute value, so greater than instead of less than
        Mguess = Eguess - e*sin(Eguess)
        f = Mtrue - (Eguess - e*sin(Eguess))
        fprime = e*cos(Eguess) - 1
        Eguess = Eguess - (f / fprime) #code did computation on one line with syntax errors
    return Eguess

sqrtmu = 0.01720209895
mu = sqrtmu**2
time = 0
dt = .05
period = sqrt((4*(pi**2)*(a**3))/mu) #added parentheses
r1ecliptic = vector(0, 0, 0)
Mtrue = 2*pi/period*(time) + M
Etrue = solvekep(Mtrue)

# Following function will use numpy matrix multiplication to rotate initial coordinates
# This is as opposed to what the sample did - brute compenent calculations on a single line
# Referenced the sample ephemeris generation guide
def matrixTransform():
    cartesian_x = a*cos(Etrue) - a*e
    cartesian_y = a*sqrt(1 - (e ** 2)) * sin(Etrue)
    cartesian_z = 0
    cartesian_vec = np.array([[cartesian_x], [cartesian_y], [cartesian_z]])

    w_rot_mat = np.array([[cos(wprime),-sin(wprime),0], [sin(wprime),cos(wprime),0], [0,0,1]])
    i_rot_mat = np.array([[1,0,0], [0,cos(iprime),-sin(iprime)], [0,sin(iprime),cos(iprime)]])
    O_rot_mat = np.array([[cos(Oprime),-sin(Oprime),0], [sin(Oprime),cos(Oprime),0], [0,0,1]])

    first_transform = np.matmul(w_rot_mat,cartesian_vec)
    second_transform = np.matmul(i_rot_mat,first_transform)
    third_transform = np.matmul(O_rot_mat, second_transform)

    position = third_transform
    r1ecliptic.x = position[0,0]
    r1ecliptic.y = position[1,0]
    r1ecliptic.z = position[2,0]
    return r1ecliptic

initial_ecliptic = matrixTransform()
asteroid = sphere(pos=vector(initial_ecliptic)*150, radius=15, color=color.red)
asteroid.trail = curve(color=color.white)
sun = sphere(pos=vector(0,0,0), radius=50, color=color.yellow)
# When using pos in vpython, need to specify that it is a vector using vector()
while (1==1):
    rate(200)
    time = time + 1
    Mtrue = 2*pi/period*(time) + M
    Etrue = solvekep(Mtrue)
    current_ecliptic = matrixTransform()
    # Function necessary due to repetition; accounts for initial +current ecliptics
    asteroid.pos = vector(current_ecliptic)*150 #Used vector()
    asteroid.trail.append(pos=asteroid.pos)

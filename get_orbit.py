#!/usr/local/bin/python
#       orbit1.py
#  Integrates the Equations of Motion ( using Euler's method)
#  solves the gravitational motion of the earth in cartesian coordinates
#
#  --------------------------------------------------------------
#   vector coordinates
#   r_earth = ( x, y)
#
#                ^
#                |
#                |
#                |
#                |                 X (x,y)
#                |
#    ------------+------------------------->
#                |
#                |
#                |
#                |
#
#
import string,sys,math
import aphelion
import solve_trajectory
import constants
print 'Planetary motion programme'

#  Mearth := mass of the earth /10^(24) kg
#  length unit = 10^6 km

Mearth = 5.97       # 10^ 24 kg

#
#  Earth's Initial Conditions: x0, vx0, y0, vy0
#

x0 = 1.471e2  # 10^6km     perihelion distance

vy0 = 2.6168  #10^6 km/day   Speed at the perihelion

aphelion.solveAphelion( x0, vy0)

solve_trajectory.getTrajectory(Mearth, x0, vy0)




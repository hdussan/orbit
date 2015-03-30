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
print 'Planetary motion programme'

#  Mearth := mass of the earth /10^(24) kg
#  GMsun := G \times Msun
#  length unit = 10^6 km

Mearth = 5.97       # 10^ 24 kg
GMsun   = 9.9e2     # (10^6 km)^3 / day^2

#
#  Earth's Initial Conditions: x0, vx0, y0, vy0
#

x0 = 1.471e2  # 10^6km     perihelion distance

vy0 = 2.6168  #10^6 km/day   Speed at the perihelion

aphelion.solveAphelion(GMsun, x0, vy0)

solve_trajectory.getTrajectory(Mearth, GMsun, x0, vy0)




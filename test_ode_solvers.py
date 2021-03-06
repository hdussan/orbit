#!/usr/local/bin/python
#
# Test for methods to solve differential Equation
#
"""
  Only one test until now. Using harmonic oscillator.
  (x  , y) : coordinates of the particle
  (vx , vy): velocity of the particle
"""
import math, sys, string
import numpy as np
import ode_solvers
from ode_solvers import *

def rhsFunction(x, vx, t):
  return  -x

x0 = 1.
vx0 = 0.
t0 = 0.
tf = 10.
oscillator = ordinaryDiffEq1D(rhsFunction, x0, vx0, t0, tf, 40)
tgrid = ode_solvers.makeTimeGrid(t0, tf, 40)

x,vx = oscillator.euler(tgrid)
x2,vx2 = oscillator.rungeKutta2order(tgrid)
x3,vx3 = oscillator.rungeKutta4order(tgrid)
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
xAnalitic = np.cos(tgrid)
vxAnalitic = -np.sin(tgrid)
figure = plt.figure()
plt.plot(tgrid, x, ':')
plt.plot(tgrid, x2, '+')
plt.plot(tgrid, x3, 'o')
plt.plot(tgrid, xAnalitic, '-')
plt.show()



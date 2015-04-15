#!/usr/local/bin/python
#       solveTrajectory.py
#  Solves the equation of motion given the initial
#  velocity and position at the perihelion of the
#  earth. (It uses simple Euler's integration)
#  TODO: Plot results of the trajectory ...
#

import math, sys, string
import numpy as np
import constants
import ode_solvers
from ode_solvers import *
import plot_orbit
def acceleration(x, y, vx, vy, t):
  GMsun = constants.GMsun
  r2 = x * x + y * y
  r3 = math.pow(r2, 3. / 2.)
  ax = -GMsun * x / r3
  ay = -GMsun * y / r3
  return ax,ay

def kineticEnergy(Mass, vx, vy):
  v2 = vx * vx + vy * vy
  return Mass * v2 / 2.0

def potentialEnergy(Mass, x, y):
  GMsun = constants.GMsun
  r = math.sqrt( x * x + y * y)
  return - GMsun * Mass / r

def mechanicalEnergy(Mass, x, y, vx, vy):
  energy0 = kineticEnergy(Mass, vx[0], vy[0]) + potentialEnergy(Mass, x[0], y[0])
  energy = np.array([energy0] * x.size)
  for i in xrange( x.size):
    energy[i] = kineticEnergy(Mass, vx[i], vy[i]) + potentialEnergy(Mass, x[i], y[i])
  return energy



def getTrajectory(Mass, positionAtPerihelion, velocityAtPerihelion):
  GMsun = constants.GMsun
  # initial conditions
  x0 = positionAtPerihelion # 10^6km     perihelion distance
  y0 =  0.0
  vx0 = 0.0
  vy0 = velocityAtPerihelion  #10^6 km/day   Speed at the perihelion

  t0 = 0.0 # days
  tf = 365.0 # days
  numberSteps = 365
  tgrid = ode_solvers.makeTimeGrid(t0, tf, numberSteps)
  planetaryOrbit = ordinaryDiffEq2D(acceleration, x0, y0, vx0, vy0, numberSteps)

  x,y,vx,vy = planetaryOrbit.rungeKutta4order(tgrid)
  energy = mechanicalEnergy(Mass, x, y, vx, vy)

  x2,y2,vx2,vy2 = planetaryOrbit.rungeKutta2order(tgrid)
  energy2 = mechanicalEnergy(Mass, x2, y2, vx2, vy2)

  plot_orbit.compareOrbits(x, y, x2, y2)
  plot_orbit.energyCheckCompare(energy, energy2, tgrid)




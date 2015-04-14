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
  x,y,vx,vy = planetaryOrbit.euler(tgrid)
  energy = mechanicalEnergy(Mass, x, y, vx, vy)

  plot_orbit.pathOrbit(x, y)
  plot_orbit.energyCheck(energy, tgrid)

  """
  from matplotlib import pyplot as plt
  from matplotlib.lines import Line2D
  figure = plt.figure()
  plt.plot(x, y, '+')
  plt.show()
  """
  """
  # initialise kinematic variables for integration
  t_next  = t0
  vx_next = vx0
  vy_next = vy0
  x_next  = x0
  y_next  = y0
  r2 = x_next * x_next + y_next * y_next
  r3 = math.pow(r2, 3. / 2.)

  # initialise energy to keep track of conservation

  kineticEnergy = 0.5 * Mass * (vx_next * vx_next + vy_next * vy_next )
  potentialEnergy= - GMsun * Mass/ math.sqrt(r2)
  Etotal0 = kineticEnergy + potentialEnergy
  Etotal = Etotal0

  counter =0
  dt = 1.0 # step in days
  # arrays for plot
  xlist = []
  ylist = []
  vxlist = []
  vylist = []
  energyList = []
  timeList =[]
  #first data point t = 0 days
  xlist.append(x_next)
  ylist.append(y_next)
  vxlist.append(vx_next)
  vylist.append(vy_next)
  timeList.append(0)
  energyList.append(Etotal0)


  f = open('planetMoves2.dat','w')
  f.write('t [d],   x [10^6 km] ,  y[10^6 km], vx [10^6 km/d], vy[10^6 km/d], Energy\n')

  while counter < 370:

    print >>f, t_next,x_next,y_next,vx_next,vy_next,r3,Etotal
    t_next += dt

    axNext, ayNext = acceleration(x_next, y_next, vx_next, vy_next, t_next)
    vx_next += axNext * dt
    vy_next += ayNext * dt

    x_next += vx_next * dt
    y_next += vy_next * dt

    r2 = x_next * x_next + y_next * y_next
    r3 = math.pow(r2, 3./2.)

    xlist.append(x_next)
    ylist.append(y_next)
    vxlist.append(vx_next)
    vylist.append(vy_next)

    kineticEnergy = 0.5 * Mass * (vx_next * vx_next + vy_next * vy_next )
    potentialEnergy =  - GMsun * Mass / math.sqrt(r2)
    Etotal = kineticEnergy + potentialEnergy
    energyList.append(Etotal)
    timeList.append(t_next)

    counter += 1

  f.close()

  plot_orbit.pathOrbit(xlist, ylist)
  plot_orbit.energyCheck(energyList, timeList)
"""




#!/usr/local/bin/python
#       solveTrajectory.py
#  Solves the equation of motion given the initial
#  velocity and position at the perihelion of the
#  earth. (It uses simple Euler's integration)
#  TODO: Plot results of the trajectory ...
#

import math, sys, string
import numpy as np
import plot_orbit
def getTrajectory(Mass, GMsun, positionAtPerihelion, velocityAtPerihelion):
  # initial conditions
  x0 = positionAtPerihelion # 10^6km     perihelion distance
  y0 =  0.0
  vx0 = 0.0
  vy0 = velocityAtPerihelion  #10^6 km/day   Speed at the perihelion

  t0 = 0.0 # days
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
  Etotal = kineticEnergy + potentialEnergy

  counter =0
  dt = 1.0 # step in days
  # arrays for plot
  xlist = []
  ylist = []
  vxlist = []
  vylist = []

  f = open('planetMoves2.dat','w')
  f.write('t [d],   x [10^6 km] ,  y[10^6 km], vx [10^6 km/d], vy[10^6 km/d], Energy\n')

  while counter < 370:

    print >>f, t_next,x_next,y_next,vx_next,vy_next,r3,Etotal
    t_next += dt

    vx_next -= GMsun * x_next * dt / r3
    vy_next -= GMsun * y_next * dt / r3

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
    counter += 1

  f.close()

  plot_orbit.pathOrbit(xlist, ylist)





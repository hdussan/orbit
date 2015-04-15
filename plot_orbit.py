#!/usr/local/bin/python
#  Plots orbit of the planet about the sun
#  given the position components
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
# Plot orbit
def pathOrbit(x, y):
  plotName = 'orbit_trajectory.png'
  figure = plt.figure()
  plt.axis('scaled')
  plt.xlim(-160, 160)
  plt.ylim(-160, 160)
  plt.title('Trajectory')
  plt.xlabel('x [$10^6$ km]')
  plt.ylabel('y [$10^6$ km]')
  plt.text(4,4, r'Sun')
  plt.grid(True)
  plt.plot(x, y, '-')
  plt.plot(0,0,'o')
  #plt.savefig(plotName)
  plt.show()

def velocities(vx,vy, t):
  plotName = 'velocity_vs_time.png'
  figure = plt.figure()
  plt.xlabel('t [days]')
  plt.ylabel('v [$10^6$ km /day]')
  plt.plot(t, vx, '-')
  plt.plot(t,vy, 'r--')
  plt.show()


def compareOrbits(x1, y1, x2, y2):
  plotName = 'orbit_trajectories.png'
  figure = plt.figure()
  plt.axis('scaled')
  plt.xlim(-160, 160)
  plt.ylim(-160, 160)
  plt.title('Trajectory')
  plt.xlabel('x [$10^6$ km]')
  plt.ylabel('y [$10^6$ km]')
  plt.text(4,4, r'Sun')
  plt.grid(True)
  plt.plot(x1, y1, '-')
  plt.plot(x2, y2, 'r--')
  plt.plot(0,0,'o')
  #plt.savefig(plotName)
  plt.show()


def energyCheck(energyIteration, time):
  energy0 = energyIteration[0]
  plotName = 'energy_deviation.png'
  figure = plt.figure()
  plt.xlim(0,366)
  plt.ylim(-0.05,0.05)
  ax = figure.add_subplot(111)
  ax.set_ylabel('Energy deviation [$10^{42}$ kg km$^2$ / day$^2$]')
  ax.set_xlabel( 't [days]')
  ax.plot(time, energyIteration - energy0, 'x')
  #plt.savefig(plotName)
  plt.show()

# Deviation from the initial value of the total energy
# This is a check on how good the integration of the Equation of Motion is.
def energyCheckCompare(energyIteration1, energyIteration2, time):
  energy0 = energyIteration1[0]
  plotName = 'energy_deviation_comp.png'
  figure = plt.figure()
  plt.xlim(0,366)
  plt.ylim(-0.001,0.001)
  ax = figure.add_subplot(111)
  ax.set_ylabel('Energy deviation [$10^{42}$ kg km$^2$ / day$^2$]')
  ax.set_xlabel( 't [days]')
  ax.plot(time, energyIteration1 - energy0, 'x', label = "Runge-Kutta 2nd order")
  ax.plot(time, energyIteration2 - energy0, '+', label = "Euler")
  #plt.savefig(plotName)
  plt.show()


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
  plt.xlim(-160, 160)
  plt.ylim(-160, 160)
  plt.title('Trajectory')
  plt.xlabel('x [$10^6$ km]')
  plt.ylabel('y [$10^6$ km]')
  plt.text(4,4, r'Sun')
  plt.grid(True)
  plt.plot(x, y, '-')
  plt.plot(0,0,'o')

  #plt.subplot(211)
  #ax = figure.add_subplot(111)
  #ax.set_ylabel('y[10^6 km]')
  #ax.set_xlabel('x[10^6 km]')

  #ax.plot(x, y,'-')
  #ax.plot(0,0,'o')
  #plt.savefig(plotName)
  plt.show()


def energyCheck(energyIteration, time):
  energy0 = energyIteration[0]
  plotName = 'energy_time.png'
  figure = plt.figure()
  plt.xlim(0,366)
  plt.ylim(-0.1,0.1)
  ax = figure.add_subplot(111)
  ax.set_ylabel('Energy deviation [$10^42$ kg km$^2$ / day$^2$]')
  ax.set_xlabel( 't [days]')
  ax.plot(time, energyIteration - energy0, 'x')
  #plt.savefig(plotName)
  plt.show()


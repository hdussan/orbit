#!/usr/local/bin/python
#  Plots orbit of the planet about the sun
#  given the position components
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
def pathOrbit(xlist, ylist):
  # Plot orbit
  x = np.array(xlist)
  y = np.array(ylist)
  figure = plt.figure()
  ax = figure.add_subplot(111)
  ax.set_ylabel('y[10^6 km]')
  ax.set_xlabel('x[10^6 km]')
  ax.plot(x, y,'-')
  plt.show()

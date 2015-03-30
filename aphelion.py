#!/usr/local/bin/python
#       aphelion.py
#
#    Solution of planetary motion by conservation
#    of energy. This method to find physical quantities
#    at the aphelion.
#    At the aphelion the velocity (vAphelion) satisfies:
#      vAphelion^2 + B vAphelion + C = 0
#    with
#      B =  - 2 G Msun / (v1 * r1)
#      C = -( v1^2 - 2 G Msun/ r1 )
#     where
#       v1 = velocityAtPerihelion
#       r1 = positionAtPerihelion
#     Inputs:
#         GMsun (G times the mass of the sun) [(10^6 km)^3 / day^2 ]
#         positionAtPerihelion [10^6 km]
#         velocityAtPerihelion [10^6 km/day]
#
#
import math
import numpy as np
def solveAphelion( GMsun, positionAtPerihelion, velocityAtPerihelion):
  vAphelion = 0.
  orbitPeriod = 0.
  eccentricity = 0.
  semiMinorAxis = 0.
  semiMajorAxis = 0.

  r1 = positionAtPerihelion
  v1 = velocityAtPerihelion
  perihelionConst = - 2. * GMsun / r1
  B = perihelionConst/ v1
  C = -( v1 * v1 + perihelionConst )

  #first check if there are real solutions:
  discriminant = B * B - 4. * C
  if discriminant >= 0. :
    solutionPlus = ( - B  + math.sqrt(discriminant) ) / 2.
    solutionMinus = ( - B  - math.sqrt(discriminant) ) / 2.
    vAphelion = min(solutionPlus, solutionMinus)
    rAphelion = v1 * r1 / vAphelion
    semiMajorAxis = (rAphelion + r1) / 2.
    semiMinorAxis = math.sqrt(rAphelion * r1)
    orbitPeriod = 2. * math.pi * semiMinorAxis * semiMajorAxis / (r1 * v1)
    eccentricity = (rAphelion - r1)/(rAphelion + r1)
  else :
    print "No real orbit for that initial conditions"


  print " velocity at aphelion = ", vAphelion," x 10^6 km/day"
  print " aphelion distance =  ", rAphelion," x 10^6 km"
  print " Semi-major Axis = ", semiMajorAxis," x 10^6 km"
  print " Semi-Minor Axis = ", semiMinorAxis," x 10^6 km"
  print " Period of the Orbit = ", orbitPeriod, " days"
  print " eccentricity =  ", eccentricity




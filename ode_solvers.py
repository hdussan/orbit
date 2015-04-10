#!/usr/local/bin/python
#
#  Traditional methods to solve numerically
#  Ordinary Differential Equation
#
"""
  Start to write methods to solve ode that
  can be re-used in several situations
  Notation is based on the traditional from the
  notation used for classical equations of motion:
  (x  , y) : coordinates of the particle
  (vx , vy): velocity of the particle
"""
import math, sys, string
import numpy as np
class ordinaryDiffEq1D:
  def __init__(self, rhsFunction, initialX, initialVx, initialTime, finalTime, numberSteps):
    self.a = rhsFunction
    self.x0 = initialX
    self.vx0 = initialVx
    self.t0 = initialTime
    self.tf = finalTime
    self.nsteps = numberSteps

  def makeTImeGrid(self):
    t = np.array([t0] * self.nsteps)
    dt = (self.tf - self.t0) / self.nsteps
    t[0] = self.t0
    for i in range(1, self.nsteps):
      t[i] = t[i - 1] + dt
    return t

  def euler(self, t):
    x = np.array([x0] * self.nsteps)
    vx = np.array([vx0] * self.nsteps)
    dt = t[1] - t[0]

    tNext = self.t0
    counter = 1
    for i in xrange( self.nsteps - 1):
      vx[i + 1] = vx[i] + self.a(x[i], vx[i], t[i]) * dt
      x[i + 1] = x[i] + vx[i + 1] * dt
    return x,vx

  def rungeKutta2order(self, t):
    x = np.array([x0] * self.nsteps)
    vx = np.array([vx0] * self.nsteps)
    dt = t[1] - t[0]
    for i in xrange(self.nsteps - 1):
      vk1 = self.a(x[i], vx[i], t[i]) * dt
      xk1 = vx[i] * dt
      vk2 = self.a(x[i] + xk1 / 2., vx[i] + vk1 / 2., t[i] + dt / 2.) * dt
      xk2 = (vx[i] + vk1 / 2.) * dt

      vx[i + 1] = vx[i] + vk2
      x[i + 1] = x[i] + xk2

    return x,vx

  def rungeKutta4order(self, t):
    x = np.array([x0] * self.nsteps)
    vx = np.array([vx0] * self.nsteps)
    dt = t[1] - t[0]
    for i in xrange(self.nsteps - 1):
      vk1 = self.a(x[i], vx[i], t[i]) * dt
      xk1 = vx[i] * dt
      vk2 = self.a(x[i] + xk1 / 2., vx[i] + vk1 / 2., t[i] + dt / 2.) * dt
      xk2 = (vx[i] + vk1 / 2.) * dt
      vk3 = self.a(x[i] + xk2 / 2., vx[i] + vk2 / 2., t[i] + dt / 2.) * dt
      xk3 = (vx[i] + vk2 / 2.) * dt
      vk4 = self.a(x[i] + xk3, vx[i] + vk3, t[i + 1]) * dt
      xk4 = (vx[i] + vk3) * dt

      vx[i + 1] = vx[i] + (vk1 + 2. * vk2 +  2. * vk3 + vk4) / 6.
      x[i + 1] = x[i] + (xk1 + 2. * xk2 +  2. * xk3 + xk4) / 6.

    return x,vx

class ordinaryDiffEq2D:
  def __init__(self, rhsFunctionX, rhFunctionY, initialX, initialY, initialVx, initialVy, initialTime, finalTime, numberSteps):
    self.ax = rhsFunctionX
    self.ay = rhsFunctionY
    self.x0 = initialX
    self.y0 = initialY
    self.vx0 = initialVx
    self.vy0 = initialVy
    self.t0 = initialTime
    self.tf = finalTime
    self.nsteps = numberSteps

  def printInitialConditions(self):
     print ('x(t = t0) = %d' % self.x0)


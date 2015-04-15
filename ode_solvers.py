#!/usr/local/bin/python
#
#  Traditional methods to solve numerically
#  Ordinary Differential Equation
#
"""
  Start to write methods to solve ode that
  can be re-used in several situations
  Notation is based on the traditional one
  used in classical mechanics equations of motion:
  (x  , y) : coordinates of the particle
  (vx , vy): velocity of the particle
"""
import math, sys, string
import numpy as np

def makeTimeGrid(initialTime, finalTime, numberSteps):
    t0 = initialTime
    tf = finalTime
    nsteps = numberSteps
    t = np.array([t0] * nsteps)
    dt = (tf - t0) / nsteps
    t[0] = t0
    for i in range(1, nsteps):
      t[i] = t[i - 1] + dt
    return t

class ordinaryDiffEq1D:
  def __init__(self, rhsFunction, initialX, initialVx, initialTime, finalTime, numberSteps):
    self.a = rhsFunction
    self.x0 = initialX
    self.vx0 = initialVx
    self.t0 = initialTime
    self.tf = finalTime
    self.nsteps = numberSteps

  def euler(self, t):
    x = np.array([self.x0] * self.nsteps)
    vx = np.array([self.vx0] * self.nsteps)
    dt = t[1] - t[0]

    for i in xrange( self.nsteps - 1):
      vx[i + 1] = vx[i] + self.a(x[i], vx[i], t[i]) * dt
      x[i + 1] = x[i] + vx[i + 1] * dt
    return x,vx

  def rungeKutta2order(self, t):
    x = np.array([self.x0] * self.nsteps)
    vx = np.array([self.vx0] * self.nsteps)
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
    x = np.array([self.x0] * self.nsteps)
    vx = np.array([self.vx0] * self.nsteps)
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
  def __init__(self, rhsFunction, initialX, initialY, initialVx, initialVy, numberSteps):
    #self.ax = rhsFunctionX
    #self.ay = rhsFunctionY
    self.acceleration = rhsFunction
    self.x0 = initialX
    self.y0 = initialY
    self.vx0 = initialVx
    self.vy0 = initialVy
    self.nsteps = numberSteps

  def printMembers(self):
    print self.vx0, ", ", self.vy0
    print self.x0,  ", ", self.y0

  def euler(self, t):
    x = np.array([self.x0] * self.nsteps)
    y = np.array([self.y0] * self.nsteps)
    vx = np.array([self.vx0] * self.nsteps)
    vy = np.array([self.vy0] * self.nsteps)
    dt = t[1] - t[0]
    for i in xrange( self.nsteps - 1):
      ax,ay = self.acceleration(x[i], y[i], vx[i], vy[i], t[i])
      vx[i + 1] = vx[i] + ax * dt
      vy[i + 1] = vy[i] + ay * dt
      x[i + 1] = x[i] + vx[i + 1] * dt
      y[i + 1] = y[i] + vy[i + 1] * dt

    return x,y,vx,vy

  def rungeKutta2order(self, t):
    x = np.array([self.x0] * self.nsteps)
    y = np.array([self.y0] * self.nsteps)
    vx = np.array([self.vx0] * self.nsteps)
    vy = np.array([self.vy0] * self.nsteps)
    dt = t[1] - t[0]
    for i in xrange(self.nsteps - 1):
      ax1,ay1 = self.acceleration(x[i], y[i], vx[i], vy[i], t[i])
      vxk1 = ax1 * dt
      xk1  = vx[i] * dt
      vyk1 = ay1 * dt
      yk1  = vy[i] * dt

      ax2,ay2 = self.acceleration(x[i] + xk1 / 2., y[i] + yk1 / 2., vx[i] + vxk1 / 2., vy[i] + vyk1 / 2., t[i + 1])
      vxk2 = ax2 * dt
      xk2  = (vx[i] + vxk1 / 2.) * dt
      vyk2 = ay2 * dt
      yk2 = (vy[i] + vyk1 / 2.) * dt

      vx[i + 1] = vx[i] + vxk2
      x[i + 1]  = x[i] + xk2
      vy[i + 1] = vy[i] + vyk2
      y[i + 1]  = y[i] + yk2

    return x,y,vx,vy

  def rungeKutta4order(self, t):
    x = np.array([self.x0] * self.nsteps)
    y = np.array([self.y0] * self.nsteps)
    vx = np.array([self.vx0] * self.nsteps)
    vy = np.array([self.vy0] * self.nsteps)
    dt = t[1] - t[0]
    for i in xrange(self.nsteps - 1):
      ax1,ay1 = self.acceleration(x[i], y[i], vx[i], vy[i], t[i])
      vxk1 = ax1 * dt
      xk1 = vx[i] * dt
      vyk1 = ay1 * dt
      yk1 = vy[i] * dt

      ax2,ay2 = self.acceleration(x[i] + xk1 / 2., y[i] + yk1 / 2., vx[i] + vxk1 / 2., vy[i] + vyk1 / 2., t[i] + dt / 2.)
      vxk2 = ax2 * dt
      xk2 = (vx[i] + vxk1 / 2.) * dt
      vyk2 = ay2 * dt
      yk2 = (vy[i] + vyk1 / 2.) * dt

      ax3,ay3 = self.acceleration(x[i] + xk2 / 2., y[i] + yk2 / 2., vx[i] + vxk2 / 2., vy[i] + vyk2 / 2., t[i] + dt / 2.)
      vxk3 = ax3 * dt
      xk3 = (vx[i] + vxk2 / 2.) * dt
      vyk3 = ay3 * dt
      yk3 = (vy[i] + vyk2 / 2.) * dt

      ax4,ay4 = self.acceleration(x[i] + xk3, y[i] + yk3, vx[i] + vxk3, vy[i] + vyk3 , t[i + 1])
      vxk4 = ax4 * dt
      xk4 = (vx[i] + vxk3) * dt
      vyk4 = ay4 * dt
      yk4 = (vy[i] + vyk3) * dt

      vx[i + 1] = vx[i] + (vxk1 + 2. * vxk2 +  2. * vxk3 + vxk4) / 6.
      x[i + 1]  = x[i]  + (xk1  + 2. * xk2  +  2. * xk3  + xk4) / 6.
      vy[i + 1] = vy[i] + (vyk1 + 2. * vyk2 +  2. * vyk3 + vyk4) / 6.
      y[i + 1]  = y[i]  + (yk1  + 2. * yk2  +  2. * yk3  + yk4) / 6.

    return x,y,vx,vy




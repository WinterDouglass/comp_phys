import numpy as np
from matplotlib import pyplot
from matplotlib.colors import ColorConverter as cc
import math

class particle2(object):
    
    def __init__(self, mass=1., x=0., y=0., vx=0., vy=0.):
        self.mass = mass
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
       
    def euler(self, fx, fy, dt):
        self.vx = self.vx + fx*dt
        self.vy = self.vy + fy*dt
        self.x = self.x + self.vx*dt
        self.y = self.y + self.vy*dt
        
    def get_force(self):  # from Lennard-Jones Potential
        r = math.sqrt(self.x*self.x+self.y*self.y)
        r8 = r**8
        r14 = r**14
        frr = 24./r8 - 48./r14 # dV/dr * 1/r 
        fx = self.x * frr
        fy = self.y * frr
        return (fx,fy)
        
    def verlet(self, dt):
        (fx,fy) = self.get_force() # before I move to the new position
        self.x += self.vx*dt + 0.5*fx*dt*dt
        self.y += self.vy*dt + 0.5*fy*dt*dt
        self.vx += 0.5*fx*dt
        self.vy += 0.5*fy*dt
        (fx,fy) = self.get_force() # after I move to the new position
        self.vx += 0.5*fx*dt
        self.vy += 0.5*fy*dt

xDelta = 0.01 # step size for x values


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
         
    def getForce(self):  # from Lennard-Jones Potential
        r = math.sqrt(self.x*self.x+self.y*self.y)
        r8 = r**8
        r14 = r**14
        frr = 24./r8 - 48./r14 # dV/dr * 1/r 
        fx = self.x * frr
        fy = self.y * frr
        return (fx,fy)    
        
    def verlet(self, dt):
        (fx,fy) = self.getForce() # before I move to the new position
        self.x += self.vx*dt + 0.5*fx*dt*dt
        self.y += self.vy*dt + 0.5*fy*dt*dt
        self.vx += 0.5*fx*dt
        self.vy += 0.5*fy*dt
        (fx,fy) = self.getForce() # after move to the new position
        self.vx += 0.5*fx*dt
        self.vy += 0.5*fy*dt
        
###################################
#            Parameters
nRuns = 40 # NUMBER OF SIMULATIONS
dt = 0.1 # TIME STEP
tMax = 13.1 # SIMULATION RUN TIME
x0 = -10. # BEGINING X POSITION
E0 = 1. # PARTICLE ENERGY (TO CALCULATE v0x SINCE WE START WITH ONLY KINETIC ENERGY
v0y = 0. # DON'T WANT ANY INITIAL Y VELOCITY
###################################
# Initializations
b = np.zeros(nRuns)
theta = np.zeros(nRuns)
v0x = math.sqrt(2*E0)
nSteps = int(tMax/dt) 
x = np.zeros(nSteps)
y = np.zeros(nSteps)
vx = np.zeros(nSteps) 
vy = np.zeros(nSteps)

for j in range(0,nRuns):
    b[j] = (j+1)*(2.5/nRuns)
    p = particle2(1., x0, b[j], v0x, v0y)
    for i in range(0,nSteps):
        p.verlet(dt)
        x[i] = p.x
        y[i] = p.y
        vx[i] = p.vx
        vy[i] = p.vy

    tanTheta = vy[nSteps-1]/vx[nSteps-1]
    if tanTheta < 0:
        theta[j] = np.pi + np.arctan(tanTheta)
    else:
        theta[j] = np.arctan(tanTheta)
pyplot.figure()
plt = pyplot.plot(b,theta)
pyplot.ylabel('$\Theta$')
pyplot.xlabel('b');
pyplot.show()

diffCrossSection = np.zeros(nRuns-1)
thetaForPlot = np.zeros(nRuns-1)
for i in range(0,nRuns-1):
    thetaForPlot[i] = theta[i]
    diffCrossSection[i] = (b[i]/np.sin(theta[i])*((b[i+1]-b[i])/(theta[i+1]-theta[i])))
pyplot.figure()
plt = pyplot.plot(thetaForPlot, diffCrossSection)
pyplot.xlabel('$\Theta$')
pyplot.ylabel('$d\sigma/d\Omega$');
pyplot.show()
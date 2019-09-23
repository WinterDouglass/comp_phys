
# challenge 1.3
class particle(object):
    
    def __init__(self, mass=1., y=0., v=0.):
        self.mass = mass
        self.y = y
        self.v = v
        
    def euler(self, f, dt):
        self.y = self.y + self.v*dt
        self.v = self.v + f/self.mass*dt
        
    def euler_cromer(self, f, dt):
        self.v = self.v + f/self.mass*dt
        self.y = self.y + self.v*dt

import numpy as np
from matplotlib import pyplot

g = 9.8            # g acceleration
mass = 0.01        # mass of the particle
y0 = 300.          # initial position
v0 = 0.            # initial velocity
R = 6370000	   # Radius of Earth

dt = 0.5           # time step

gforce = g*mass    # weight
results = [y0,v0,0.,y0,v0,0.,0.]

# yDepForce = g*mass/((1 + y/R)**2)
p = particle(mass, y0, v0)
newP = particle(mass, y0, v0)

#y = [y0] # since we do not know the size of the arrays, we define first a python list
#v = [v0] # the append method is more efficient for lists than arrays
#t = [0.]

# I'm not using arrays because we only care about the last point so why take up memory not being used

def height(y0,p,newP):
    while p.y > 0.:
        fy = -gforce
        p.euler(fy, dt)
        results[0] = p.y
        results[1] = p.v
        results[2] = results[2]+dt
    
    while newP.newY > 0.:
        yDepForce = g*mass/((1 + results[3]/R)**2)
        newP.euler(yDepForce, dt)
        results[3] = newP.newY
        results[4] = newP.newV
        results[5] = results[5]+dt
    
    results[6] = abs((results[1] - results[4])/results[1])
    return(results)
results = height(y0,p,newP)
# set some bounds to do a binary search
yH = 100000000. # upper limit
yL = 1. # lower limit
while True:
    if(results[6]<.01):
        yL = y0
        y0 = (y0+yH)/2
        p = particle(mass,y0,v0)
        newP = particle(mass,y0,v0)
        results = height(y0,p,newP)
    elif(results[6]>0.1):
        yH = y0
        y0 = (y0+yL)/2
        p = particle(mass,y0,v0)
        newP = particle(mass,y0,v0)
    elif(results[6]<0.100001 or results[6]>0.999999):
        break

print(y0)    
    
# t_data = np.array(t) # we convert the list into a numpy array for plotting
# y_data = np.array(y)
# v_data = np.array(v)

# pyplot.plot(t_data, v_data, color="#FF0000", ls='-', lw=3)
# pyplot.xlabel('time(s)')
# pyplot.ylabel('velocity(m/s)');

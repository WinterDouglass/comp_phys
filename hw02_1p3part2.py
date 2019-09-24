class particle2(object):
    
    def __init__(self, mass=1., x=0., y=0., vx=0., vy=0.):
        self.mass = mass
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
       
    def euler(self, fx, fy, dt):
        self.vx = self.vx + fx/self.mass*dt
        self.vy = self.vy + fy/self.mass*dt
        self.x = self.x + self.vx*dt
        self.y = self.y + self.vy*dt

import numpy as np
from matplotlib import pyplot
from matplotlib.colors import ColorConverter as cc
import math

g = 9.8            # g acceleration
v0 = 30.           # initial velocity
k = 0.1
dt = 0.1           # time step

colors = ['red','orange','yellow','green','magenta','cyan','blue','purple','black']
xM = [0]
for angle in range(1,9):
    x = [0]                                  # we need to initialize the arrays for each value of the angle
    y = [0]
    vx = [math.cos(angle*0.1*math.pi/2.)*v0] 
    vy = [math.sin(angle*0.1*math.pi/2.)*v0] 
    t = [0.]

    p = particle2(1., 0., 0., vx[0], vy[0])
    while p.y >= 0.:
        fy = -g-k*(p.vy)*math.sqrt((p.vy)**2 + (p.vx)**2)
        fx = -k*p.vx*math.sqrt((p.vy)**2 + (p.vx)**2)
        p.euler(fx, fy, dt)
        x.append(p.x)
        y.append(p.y)
        vx.append(p.vx)
        vy.append(p.vy)
        t.append(t[-1]+dt)
    
    t_data = np.array(t) # we convert the list into a numpy array for plotting
    x_data = np.array(x)
    y_data = np.array(y)
    vx_data = np.array(vx)
    vy_data = np.array(vy)
    xM.append(p.x)

    my_plot = pyplot.plot(x_data, y_data, color=(colors[angle]), ls='-', lw=3, label = str(angle*0.1))
    pyplot.legend()

pyplot.xlabel('position x(m)')
pyplot.ylabel('position y(m)');
pyplot.show()
xMax = 0.
angleMax = 0.
for i in range(1,len(xM)):
    if(xM[i-1] < xM[i]):
        xMax = xM[i]
        angleMax = i*0.1
print('Max distance: '+ str(xMax)+ ' | Angle of max distance: '+  str(angleMax*90) + ' Degrees')
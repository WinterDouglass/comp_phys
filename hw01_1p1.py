T0 = 10.   # initial temperature
Ts = 83.   # temp. of the environment
r = 0.1    # cooling rate
dt = 0.05  # time step
tmax = 60. # maximum time
nsteps = int(tmax/dt)  # number of steps

T = T0
for i in range(1,nsteps+1):
    new_T = T - r*(T-Ts)*dt
    T = new_T
    #print ('{:20.18f}  {:20.18f}  {:20.18f}'.format(i,i*dt, T))

import numpy as np
from matplotlib import pyplot 
my_time = np.zeros(nsteps)
my_temp = np.zeros(nsteps)

T = T0
my_temp[0] = T0
for i in range(1,nsteps):
    T = T - r*(T-Ts)*dt
    my_time[i] = i*dt
    my_temp[i] = T

#pyplot.plot(my_time, my_temp, color='#003366', ls='-', lw=3)
#pyplot.xlabel('time')
#pyplot.ylabel('temperature');

my_time = np.linspace(0.,tmax,nsteps)

#pyplot.plot(my_time, my_temp, color='#003366', ls='-', lw=3)
#pyplot.xlabel('time')
#pyplot.ylabel('temperature');

def euler(y, f, dx):
    """Computes y_new = y + f*dx
    
    Parameters
    ----------
    y  : float
        old value of y_n at x_n
    f  : float
        first derivative f(x,y) evaluated at (x_n,y_n)
    dx : float
        x step
    """
    
    return y + f*dx


T = T0
for i in range(1,nsteps):
    T = euler(T, -r*(T-Ts), dt)
    my_temp[i] = T


euler = lambda y, f, dx: y + f*dx 

dt = 1.
#my_color = ['#003366','#663300','#660033','#330066']
my_color = ['red', 'green', 'blue', 'black']
for j in range(0,4):
    nsteps = int(tmax/dt)    #the arrays will have different size for different time steps
    my_time = np.linspace(dt,tmax,nsteps) 
    my_temp = np.zeros(nsteps)
    T = T0
    for i in range(1,nsteps):
        T = euler(T, -r*(T-Ts), dt)
        my_temp[i] = T
        
    pyplot.plot(my_time, my_temp, color=my_color[j], ls='-', lw=3)
    dt = dt/2.

    
#pyplot.xlabel('time');
#pyplot.ylabel('temperature');
#pyplot.xlim(8,10)
#pyplot.ylim(48,58);


# here is challenge 1.1 for hw01
nPoints = 20 # number of trials to plot
endTime = 10 # the time that we are interested in
delTime = np.zeros(nPoints)
tempE = np.zeros(nPoints)
for i in range(0, nPoints):
    dt = 1./(i+1) #dt ranges from 1/20 to 1
    delTime[i] = dt #list of dt values
    T = T0
    nSteps = int(endTime/dt)
    for j in range(nSteps):
        T = euler(T, -r*(T-Ts), dt) #the last one will be T ans t = 10
    tempE[i] = T #recording the last found temperature at t = 10 into the list

print (temp)
print (delTime)
pyplot.plot(delTime, temp, 'k.')
pyplot.xlabel('delta time');
pyplot.ylabel('temperature');


# this is for challenge 1.2 I just figured I would keep it all together
# the program should only print the results for the homework problem

def rK(T, r):
    return (-r*(T-Ts))
dt = 1.
T10 = Ts + (T0 - Ts)*np.exp(-r*(10))
temp = np.zeros(nPoints)
delTime = np.zeros(nPoints)
for i in range(nPoints):
    dt = 1./(i+1) #dt ranges from 1/20 to 1
    delTime[i] = dt #list of dt values
    nsteps = int(endTime/dt)    #the arrays will have different size for different time steps
    my_time = np.linspace(dt,tmax,nsteps) 
    my_temp = np.zeros(nsteps)
    T = T0
    for j in range(nsteps):
        k1 = rK(T, r)*dt
        k2 = rK(T + k1/2, r)*dt
        k3 = rK(T + k2/2, r)*dt
        k4 = rK(T + k3, r)*dt
        T = T + 1./6. * (k1 + 2*k2 + 2*k3 + k4)
    #print (T - T10)
    temp[i] = T - T10
    
#pyplot.plot(delTime, temp, 'k.')
#pyplot.ylim(-.0000004,0.0000001);


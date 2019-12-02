# -*- coding: utf-8 -*-
"""
Whole Thing Now
"""
import math
import numpy as np
import time
from matplotlib import pyplot
tInit = time.time()
'''
[gridRand] - gives arbitrarily assigned vertices in an nSq by nSq grid 
in a box of side length L
[perimeter] - calculates the perimeter of a cell defined by vertices in an ordered list
[area] - area of a cell defined by an ordered list of vertices
[energy] - energy of a single cell by a function which calls area and perimeter
[energyTotal] - calls gridRand to get all the vertices then arranges them to calculate
energy for each cell and sums to get theh total
[minEnergy] - calls energyTotal, moves a vertex and recall energyTotal until it gets to a
lower energy state, moves back if energy is higher stays moved if energy lower, for all vertices
until it reaches the lowest state
[prettyPicture] - makes a pretty picture obvi
'''
def gridRand(L,nSq):
    nSq2 = nSq**2    
    x1 = np.zeros(nSq2)
    y1 = np.zeros(nSq2)
    n = 0
    x = [[0.0 for i in range(nSq+1)] for j in range(nSq+1)]
    y = [[0.0 for i in range(nSq+1)] for j in range(nSq+1)]    
    for i in range(nSq):
        for j in range(nSq):
            x[i][j] = np.random.random()*(L)/nSq + i*L/nSq
            y[i][j] = np.random.random()*(L)/nSq + j*L/nSq
            x1[n] = x[i][j]
            y1[n] = y[i][j]
            n += 1           
    # top righty is bottom lefty
    x[nSq][nSq] = x[0][0] + L
    y[nSq][nSq] = y[0][0] + L
    for k in range(nSq):
        # set all right side virtual particles to left side reals
        x[nSq][k] = x[0][k] + L
        y[nSq][k] = y[0][k]
        # set all top virtual particles to bottom reals
        x[k][nSq] = x[k][0]
        y[k][nSq] = y[k][0] + L
        
    return(x,y,x1,y1)

def perimeter(x, y):
    P = math.sqrt((x[0]-x[len(x)-1])**2 + (y[0]-y[len(x)-1])**2) # initialize with the side length of the last and first index
    for i in range(len(x)-1): # sum over the sides with adjacent indices
        P += math.sqrt((x[i]-x[i+1])**2 + (y[i]-y[i+1])**2)
    return(P)
    
def area(x, y):
    A = 0.5*(x[len(x)-1]*y[0] - y[len(y)-1]*x[0]) # order matters sum(x[i]y[i+1]-y[i]x[i+1])...(x[n]y[0]-y[n]x[0])
    for i in range(len(y)-1):
        A += 0.5*(x[i]*y[i+1] - y[i]*x[i+1])
    A = abs(A)
    return(A)

def energy(x, y):
    '''
    these variables need to be set according to some physical consideration
    E[i] = ksi*P[i]**2 + gamma*P[i] + beta*(A[i] - A0)
    '''
    ##################### Set these #######################
    ksi = 1.
    gamma = 2.
    beta = 3.
    A0 = 2.
    #######################################################
    A = area(x, y)
    P = perimeter(x, y)
    E = ksi*P**2 + gamma*P + beta*(A - A0)
    return(E)

def energyTotal(nGrid,x,y):
    totalEnergy = 0.0
    
    for i in range(nGrid):
        for j in range(nGrid):
            # start at bottom left vertex and assign an array to the vertices in ccw order
            eks=[x[i][j],x[i][j+1],x[i+1][j+1],x[i+1][j]]
            wai=[y[i][j],y[i][j+1],y[i+1][j+1],y[i+1][j]]
            totalEnergy += energy(eks,wai)
    return(totalEnergy)

def minEnergy(sideLen,nGrid,x,y,x1,y1,stepSize):
    en = energyTotal(nGrid,x,y)
    caseState = True
    while(caseState):
        eBegin = en
        #prettyPicture(x,y,nGrid)
        for i in range(nGrid):
            for j in range(nGrid):
                # move x and check energy stay mov if lower move back if higher
                x[i][j] += stepSize
                x[nGrid][nGrid] = x[0][0] + sideLen
                for k in range(nGrid):
                    x[nGrid][k] = x[0][k] + sideLen
                    x[k][nGrid] = x[k][0]
                newEn = energyTotal(nGrid,x,y)
                if (newEn - en >= 0.0):
                    x[i][j] -= 2*stepSize
                    x[nGrid][nGrid] = x[0][0] + sideLen
                    for k in range(nGrid):
                        x[nGrid][k] = x[0][k] + sideLen
                        x[k][nGrid] = x[k][0]
                    newerEn = energyTotal(nGrid,x,y)
                    if (newerEn - en >= 0.0):
                        x[i][j] += stepSize
                        x[nGrid][nGrid] = x[0][0] + sideLen
                        for k in range(nGrid):
                            x[nGrid][k] = x[0][k] + sideLen
                            x[k][nGrid] = x[k][0]
                    else:
                        en = newerEn
                else:
                    en = newEn
                # same thing with y    
                y[i][j] += stepSize
                y[nGrid][nGrid] = y[0][0] + sideLen
                for k in range(nGrid):
                    y[nGrid][k] = y[0][k]
                    y[k][nGrid] = y[k][0] + sideLen
                newEn = energyTotal(nGrid,x,y)
                if (newEn - en >= 0.0):
                    y[i][j] -= 2*stepSize
                    y[nGrid][nGrid] = y[0][0] + sideLen
                    for k in range(nGrid):
                        y[nGrid][k] = y[0][k]
                        y[k][nGrid] = y[k][0] + sideLen
                    newerEn = energyTotal(nGrid,x,y)
                    if (newerEn - en >= 0.0):
                        y[i][j] += stepSize
                        y[nGrid][nGrid] = y[0][0] + sideLen
                        for k in range(nGrid):
                            y[nGrid][k] = y[0][k]
                            y[k][nGrid] = y[k][0] + sideLen
                    else:
                        en = newerEn
                else:
                    en = newEn
        if (abs(en - eBegin) < 0.01):
            caseState = False
        elif(en > eBegin):
            print('Does not work?')
            caseState = False
    return(en,x,y)

def prettyPicture(x,y,nGrid):
    eks = [[0.0 for i in range(nGrid+1)] for j in range(nGrid+1)]
    wai = [[0.0 for i in range(nGrid+1)] for j in range(nGrid+1)]   
    for i in range(nGrid+1):
        for j in range(nGrid+1):
            eks[i][j] = x[j][i]
            wai[i][j] = y[j][i]
    for i in range(nGrid+1):
        pyplot.plot(x[i],y[i], '-k')
        pyplot.plot(eks[i],wai[i], '-k')
        pyplot.plot(x[i],y[i],'ro')
    
    pyplot.show()
    return()
    
# start the program with these commands
sideLen = 10
nGrid = 5
stepSize = 0.01

x,y,x1,y1 = gridRand(sideLen,nGrid)
prettyPicture(x,y,nGrid)
energyMinimum = minEnergy(sideLen, nGrid, x, y, x1, y1, stepSize)
prettyPicture(x,y,nGrid)
'''
pyplot.figure()
plt = pyplot.plot(x1,y1, 'k*')
#pyplot.title('')
pyplot.ylabel('$y_j$')
pyplot.xlabel('$x_i$');
pyplot.show()
'''

tFin = time.time()
tDiff = tFin - tInit
print('Run Time: ', tDiff)
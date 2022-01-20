#!/usr/bin/python
import json, math
import numpy as np 
from matplotlib import pyplot as plt

def ftcs(inpath = '2.json', outpath = '2.png', plot = True):
    # load data into variables
    with open(inpath) as f:
        data = json.load(f)
    xmax = 2
    tmax = 10

    delx = 0.05
    delt1 = 0.01
    delt2 = 0.001
    delt3 = 0.0001
    
    diff = data["alpha"]
    c0 = data["c0"]
    cMax = data["cxmax"]
    const = (diff * delt1)/(delx**2)
    init = lambda x : math.sin(math.pi*x/xmax) #initial condition (lambda exp.)

    t1 = np.arange(0, (tmax + delt1), delt1)
    t2 = np.arange(0, (tmax + delt2), delt2)
    t3 = np.arange(0, (tmax + delt3), delt3)

    x = np.arange(0, (xmax + delx), delx)
    rows = len(x)
    cols = len(t1)

    c = np.zeros((rows, cols))

    exact = lambda x,y : math.sin(math.pi*x/xmax) * math.exp(-diff*(math.pi**2)*y/(xmax**2))
    a = np.zeros(rows)
    for d in range(len(x)):
        a[d] = exact(x[d],10)

    # initial conditions
    for k in range(rows):
        c[k][0] = init(x[k])

    # fill concentration matrix incl. BCs
    for i in range(rows):
        for j in range(cols):
            if i == 0: #BC1
                c[i][j] = c0
            elif i == (rows - 1): #BC2
                c[i][j] = cMax
            elif j != 0:          
                c[i][j] = c[i][j-1] + const * (c[i-1][j-1] - 2*c[i][j-1] + c[i+1][j-1])
    
    cols2 = len(t2)
    c2 = np.zeros((rows, cols2))
    const2 = (diff * delt2)/(delx**2)
    # initial conditions
    for k in range(rows):
        c2[k][0] = init(x[k])

    # fill concentration matrix incl. BCs
    for i in range(rows):
        for j in range(cols2):
            if i == 0: #BC1
                c2[i][j] = c0
            elif i == (rows - 1): #BC2
                c2[i][j] = cMax
            elif j != 0:          
                c2[i][j] = c2[i][j-1] + const2 * (c2[i-1][j-1] - 2*c2[i][j-1] + c2[i+1][j-1])
    
    cols3 = len(t3)
    c3 = np.zeros((rows, cols3))
    const3 = (diff * delt3)/(delx**2)
    # initial conditions
    for k in range(rows):
        c3[k][0] = init(x[k])

    # fill concentration matrix incl. BCs
    for i in range(rows):
        for j in range(cols3):
            if i == 0: #BC1
                c3[i][j] = c0
            elif i == (rows - 1): #BC2
                c3[i][j] = cMax
            elif j != 0:          
                c3[i][j] = c3[i][j-1] + const3 * (c3[i-1][j-1] - 2*c3[i][j-1] + c3[i+1][j-1])
    
    # generate plot and save as png 
    initialC = c[:,0]
    tmaxC = c[:,-1]
    tmaxC2 = c2[:,-1]
    tmaxC3 = c3[:,-1]
    e1 = np.zeros(len(tmaxC))
    e2 = np.zeros(len(tmaxC))
    e3 = np.zeros(len(tmaxC))

    for b in range(len(tmaxC)):
        e1[b] = a[b] - tmaxC[b]
        e2[b] = a[b] - tmaxC2[b]
        e3[b] = a[b] - tmaxC3[b]
    
    er1 = sum(e1)/(len(e1))
    er2 = sum(e2)/len(e2)
    er3 = sum(e3)/len(e3)
    print((er1+er2+er3)/3)
    
    deltas = np.array([0.01, 0.001, 0.0001])
    error = np.array([er1, er2, er3])


    if plot:
        #plt.plot(x, initialC, 'k', linestyle = 'dashed', label = 'c(x, 0)')

        #plt.plot(deltas, error)
        #plt.plot(x, e1)
        #plt.plot(x, e2)
        #plt.plot(x, a, 'b')
        #plt.plot(x, tmaxC, 'b', label = 'delt = 0.01')
        #plt.plot(x, tmaxC2, 'k', label = 'delt = 0.')
        #plt.plot(x, tmaxC3, 'r', label = 'delt = 0.05')
        plt.xlabel('x')
        plt.ylabel('c(x, t)')
        #plt.legend()
        
        plt.show()
        #plt.savefig(outpath)
    
    return x, tmaxC

def ftcs2 (inpath = '2.json', outpath = 'ftcs_plot.png', plot = True):
    # load data into variables
    with open(inpath) as f:
        data = json.load(f)
    xmax = data["xmax"]
    tmax = data["tmax"]

    delx = data["delx"]
    delt = data["delt"]
    
    diff = data["alpha"]
    c0 = data["c0"]
    cMax = data["cxmax"]
    const = (diff * delt)/(delx**2)
    init = eval(data["ic"]) #initial condition (lambda exp.)

    t = np.arange(0, (tmax + delt), delt)
    x = np.arange(0, (xmax + delx), delx)
    rows = len(x)
    cols = len(t)
    c = np.zeros((rows, cols))

    # initial conditions
    for k in range(rows):
        c[k][0] = init(x[k])

    # fill concentration matrix incl. BCs
    for i in range(rows):
        for j in range(cols):
            if i == 0: #BC1
                c[i][j] = c0
            elif i == (rows - 1): #BC2
                c[i][j] = cMax
            elif j != 0:          
                c[i][j] = c[i][j-1] + const * (c[i-1][j-1] - 2*c[i][j-1] + c[i+1][j-1])
    
    # generate plot and save as png 
    initialC = c[:,0]
    tmaxC = c[:,-1]
    if plot:
        plt.plot(x, initialC, 'k', linestyle = 'dashed', label = 'c(x, 0)')
        plt.plot(x, tmaxC, 'b', label = 'c(x, t=tmax)')
        plt.xlabel('x')
        plt.ylabel('c(x, t)')
        plt.legend()
        plt.show()
        #plt.savefig(outpath)

    return x, tmaxC
ftcs2()
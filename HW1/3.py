# Laura (Kai Sze) Luk (kluk@andrew.cmu.edu)

#!/usr/bin/python
import json, copy, math
import numpy as np 
from matplotlib import pyplot as plt
from scipy import linalg


def btcs(inpath = 'sample_in.json', outpath = 'btcs_plot.png', plot = True):
    # load data into variables
    with open(inpath) as f:
        data = json.load(f)
    xmax = data["xmax"]
    tmax = data["tmax"]

    delx = data["delx"]
    delt = data["delt"]
    
    diff = data["alpha"]
    c0 = 0
    cMax = 0
    const = (diff * delt)/(delx**2)
    init = lambda x: 4*x - x**3 #initial condition (lambda exp.)

    t = np.arange(0, (tmax + delt), delt)
    x = np.arange(0, (xmax + delx), delx)
    rows = len(x)
    cols = len(t)

    #initial conditions
    c_m = np.zeros(rows)
    for k in range(rows):
        c_m[k] = init(x[k])
    
    initialC = copy.deepcopy(c_m)

    # set up matrix A of constants
    a = diff*delt/(delx**2)
    b = (2*diff*delt/(delx**2)) + 1
    A = np.zeros((rows, rows))
    for i in range(1, rows - 1):
        A[i, i-1] = -a 
        A[i, i+1] = -a 
        A[i, i] = b
    A[0, 0] = 1
    A[rows - 1, rows - 1] = 1

    c = np.zeros(rows)
    d = np.zeros(rows)
    for j in range(len(t)):
        for m in range(1, rows - 1):
            d[i] = c_m[i]
        d[0] = c0
        d[rows - 1] = cMax
        c[:] = linalg.solve(A, d)
        c_m[:] = c
    
    if plot:
        plt.plot(x, initialC, 'k', linestyle = 'dashed', label = 'c(x, 0)')
        plt.plot(x, c, 'b', label = 'c(x, t=tmax)')
        plt.xlabel('x')
        plt.ylabel('c(x, t)')
        plt.legend()
        plt.show()
        #plt.savefig(outpath)
    
    return x, c
btcs()
# Laura (Kai Sze) Luk (kluk@andrew.cmu.edu)

import json
import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt

def fem(inpath = 'uniform4.npy', m = 2.0, outpath = 'fem_results.json', plotfile = 'fem_plot.png', plot = True):
    # load node locations
    nodes = np.load(inpath)
    totalNodes = len(nodes)
    
    # construct array of delta x_i's
    delX = np.zeros(totalNodes - 1)
    for i in range(totalNodes - 1):
        delX[i] = nodes[i + 1] - nodes[i]

    # check delta x_i sums to 1
    assert(sum(delX) == 1)

    # build stiffness matrix
    size = totalNodes - 2
    diag = np.zeros(size)
    for i in range(size):
        value = (1/delX[i]) + (1/delX[i + 1])
        diag[i] = value

    offDiag = np.zeros(size - 1)
    for i in range(size - 1):
        offDiag[i] = -1/delX[i + 1]
    # combine diagonal and off-diagonals 
    matrix = np.diag(diag, 0) + np.diag(offDiag, -1) + np.diag(offDiag, 1)

    # build RHS column 
    rightSide = np.zeros((size, 1))
    for i in range(size):
        rightSide[i, 0] = (m / 2) * (delX[i] + delX[i + 1])

    # solve for concentration array
    c = linalg.solve(matrix, rightSide)
    # include boundary conditions
    c = np.insert(c, size, 0, 0)
    c = np.insert(c, 0, 0, 0)

    # plot results
    if plot:
        t = np.linspace(0, 1)
        exact = np.zeros(len(t))
        for i in range(len(t)):
            exact[i] = t[i] - (t[i] ** 2)
        
        plt.plot(t, exact, 'k', linestyle = 'dashed', label = 'exact')
        plt.plot(nodes, c, 'b', label = 'computed')
        plt.xlabel('x')
        plt.ylabel('c(x)')
        plt.legend()
        plt.show()
        plt.savefig(plotfile)
    
    # save results
    results = dict()
    results['x'] = nodes.tolist()
    c = c.tolist()
    results['comp'] = [item for sublist in c for item in sublist]
    with open(outpath, 'w') as fp:
        json.dump(results, fp)

    pass
#fem()





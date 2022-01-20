import numpy as np
import numpy.linalg as la
import json 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import math
from copy import deepcopy

# courtesy of TA Ryan (1st two functions)

def philj(epsilon, sigma, rc, rij):
    try: # assume rij is an array of values
        rij = rij[rij <= rc] # rij should never equal 0
        if not len(rij): # no elements are below cutoff radius
            return np.array([0.])
    except: # rij is a single numeric value
        if not rij <= rc:
            return 0.
    return 4*epsilon*((sigma/rij)**12 - (sigma/rij)**6)

def latticesumpbc(coords, ncells, rc, a, epsilon, sigma):
    # fast computation of distances between atomic coordinates
    Dij = coords.T - coords[:,:,np.newaxis]
    # apply 'minimum image convention' periodic boundary conditions
    Dij -= np.rint(Dij)

    # Rij[i,j] = distance between atomsat coordinates i and j
    Rij = np.sqrt(np.sum(Dij**2, axis=1)) * a * ncells
    # Atoms don't interact with themselves, so set Rii = 2*rc for all i (no energy contribution)
    np.fill_diagonal(Rij, rc*2)
    return (philj(epsilon, sigma, rc, Rij) / 2).sum() # sum all energies, divide by 2 to prevent double counting ij and ji

def metropolis(pe, kT):
    if (pe <= 0): return 1
    else:
        return math.exp((-pe)/kT)

def atommc(inpath='mcinput.json'):
    with open(inpath) as f:
        data = json.load(f)
    density = data["density"]
    simT = data["kT"]
    delta = data["delta"]
    totalTime = data["totalMCS"]
    epsilon = data["epsilon"]
    sigma = data["sigma"]
    rc = data["rc"]
    plotB = data["plot"]
    fileName = data["atomfile"]

    atomCoord = np.load(fileName)
    numAtoms = len(atomCoord)
    nCells = (numAtoms // 4) ** (1/3)
    latticeParam = (float(4) / density) ** (1/3)
    initPE = latticesumpbc(atomCoord, nCells, rc, latticeParam, epsilon, sigma)
    totalSteps = totalTime * numAtoms
    t = 1
    peArr = np.zeros(totalSteps)
    peArr[0] = initPE
    numSuccess = 0
 
    while (t < totalSteps):
        randAtomInd = random.randint(0, numAtoms - 1)
        (x, y, z) = atomCoord[randAtomInd]
        displaceX = random.uniform(-delta, delta)
        displaceY = random.uniform(-delta, delta)
        displaceZ = random.uniform(-delta, delta)
        oldPE = latticesumpbc(atomCoord, nCells, rc, latticeParam, epsilon, sigma)

        tmp = deepcopy(atomCoord)
        tmp[randAtomInd] = (x + displaceX, y + displaceY, z + displaceZ)
        newPE = latticesumpbc(tmp, nCells, rc, latticeParam, epsilon, sigma)
        deltaPE = newPE - oldPE

        prob = metropolis(deltaPE, simT)
        if (prob == 1):
            atomCoord[randAtomInd] = (x + displaceX, y + displaceY, z + displaceZ)
            numSuccess += 1

        else:
            prob1 = np.random.random()
            if (prob1 < prob):
                atomCoord[randAtomInd] = (x + displaceX, y + displaceY, z + displaceZ)
                numSuccess += 1
                
        peArr[t] = newPE
        t += 1

    if plotB:
        xTime = np.linspace(1, totalSteps, totalSteps)
        plt.plot(xTime, peArr)
        #plt.show()
        plt.savefig('energy.png')

        
        coords = np.asarray(atomCoord)
        plt.savefig('structure.png')
    
    successRatio = numSuccess / totalSteps
    timeAvgPE = sum(peArr) / totalSteps
    return successRatio, timeAvgPE

#atommc()

    




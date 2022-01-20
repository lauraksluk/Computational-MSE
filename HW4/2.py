import numpy as np
import json, math
import matplotlib.pyplot as plt
from copy import deepcopy

def metropolis0(dE, p0 = 1):
    if dE <= 0: return p0
    else: return 0

def delt(i, j):
    if i == j: return 1
    else: return 0

def delE(moore, newstate, Fliq = 3, Fsol = 0):
    if newstate not in moore: 
        raise AssertionError("newstate not in neighborhood")
    else:
        old = moore[1][1]
        unlike = 0
        for i in range(len(moore)):
            for j in range(len(moore)):
                if i != 1 or j != 1:
                    neigh1 = 1 - delt(moore[i][j], newstate)
                    neigh2 = 1 - delt(moore[i][j], old)
                    neigh = neigh1 - neigh2
                    unlike += neigh
        A = delt(0, old)
        B = delt(0, newstate)
        liq = (B - A) * Fliq
        solid = ((1 - B) - (1 - A)) * Fsol
        result = unlike + liq + solid
        return result

def getUniqueNeighbors(states, row, col):
    unique = set()
    startRow = row - 1
    startCol = col - 1
    endRow = row + 2
    endCol = col + 2
    for i in range(startRow, endRow):
        for j in range(startCol, endCol):
            if i != row or j != col:
                unique.add(states[i][j])
    sortedList = list(unique)
    sortedList.sort()
    return np.array(sortedList)

def mcpm_step(states, shape, rng, prob_f=metropolis0):
    row = rng.randint(1, shape[0] - 1)
    col = rng.randint(1, shape[1] - 1)
    visit = states[row][col]

    newSeed = rng.randint(0, 2**32)
    prob = np.random.random()
    
    neighbors = getUniqueNeighbors(states, row, col)
    rng.seed(newSeed)
    newState = rng.choice(neighbors)
    mooreRowS = row - 1
    mooreRowE = row + 2
    mooreColS = col - 1
    mooreColE = col + 2
    moore = states[mooreRowS : mooreRowE, mooreColS : mooreColE]
    energy = delE(moore, newState, 3, 0)
    if prob < prob_f(energy):
        states[row][col] = newState
    return states

def solidCount(data):
    count = 0
    grains = {}
    for i in range(1, len(data) - 1):
        for j in range(1, len(data[0]) - 1):
            if data[i][j] > 0:
                count += 1
                if data[i][j] in grains:
                    grains[data[i][j]] += 1
                else:
                    grains[data[i][j]] = 1
    return count, grains

def getGrains(data, id):
    count = 0
    for i in range(1, len(data) - 1):
        for j in range(1, len(data[0]) - 1):
            if data[i][j] > 0 and data[i][j] == id:
                count += 1
    return count

def grainCount(data, grains):
    data1 = deepcopy(data)
    all = data1.flatten()
    return np.bincount(all, minlength=grains+1)

def mcpm_slow(time = 80, ins = 'single3.json', plot = True, rng_seed = None):
    with open(ins) as f:
        data = json.load(f)
    rows = data["width"]
    cols = data["length"]
    init = data["initial"]
    init = np.array(init)
    rng = np.random.RandomState()
    
    solidFrac = np.zeros(time + 1)
    totalArea = (rows - 2) * (cols - 2)
    (initSolid, grains) = solidCount(init)
    solidFrac[0] = float (initSolid) / totalArea

    everyFive = (time // 5) + 1
    grainArea = np.zeros((everyFive, len(grains.keys())))
    noBuf = init[1 : rows-1, 1 : cols-1]
    grainArea[0] = grainCount(noBuf, len(grains.keys()))[1:]
    #ind = 0
    gRow = 0
    # for key in grains:
    #     num = getGrains(init, key)
    #     grainArea[gRow][key-1] = num
        #ind += 1
    #print(grainArea)
    states = deepcopy(init)
    currTime = 0
    while currTime < time:
        for i in range(totalArea):
            rng.seed(rng_seed)
            rng_seed = rng.randint(0, 2**32)
            states = mcpm_step(states, (rows, cols), rng, metropolis0)
            
        sol = solidCount(states)
        solidFrac[currTime + 1] = float (sol[0]) / totalArea

        if currTime % 5 == 0:
            gRow += 1
            noBuffer = states[1 : rows-1, 1 : cols-1]
            #grainArea[gRow] = grainCount(noBuffer, len(grains.keys()))[1:]
            # newG = sol[1]
            # index = 0
            # for k in newG:
            #     n = getGrains(states, k)
            #     grainArea[gRow][k-1] = n
       
        currTime += 1
    #print(grainArea)
    #print(solidCount(states)[0])
    if plot:
        plt.imshow(states)
        t = np.linspace(0, time, time + 1)
        #plt.plot(t, solidFrac)
        plt.show()
    
    return states, solidFrac, grainArea

#mcpm_slow()

# def test():
#     m1 = [[9,1,1],[9,1,1],[2,2,2]]
#     m = np.asarray(m1)
#     return mcpm_slow()
# test()



def mcpm_efficient(time = 5, ins = 'single3.json', plot = True, rng_seed = None):
    with open(ins) as f:
        data = json.load(f)
    rows = data["width"]
    cols = data["length"]
    init = data["initial"]
    init = np.array(init)
    rng = np.random.RandomState()
    Fliq = 3
    Fsol = 0
    
    solidFrac = np.zeros(time + 1)
    totalArea = (rows - 2) * (cols - 2)
    area = np.zeros(time + 1)   
    g = np.zeros(time + 1)
    a0 = 5
    growth = lambda t : math.pi * (t**2) + a0

    (initSolid, grains) = solidCount(init)
    solidFrac[0] = float (initSolid) / totalArea

    everyFive = (time // 5) + 1
    grainArea = np.zeros((everyFive, len(grains.keys())))
   
    gRow = 0
    noBuf = init[1 : rows-1, 1 : cols-1]
    grainArea[0] = grainCount(noBuf, len(grains.keys()))[1:]
    area[0] = initSolid

    states = deepcopy(init)
    currTime = 0
    while currTime < time:
        g[currTime] = growth(currTime)
        
        for i in range(totalArea):
            rng.seed(rng_seed)
            rng_seed = rng.randint(0, 2**32)

            row = rng.randint(1, rows - 1)
            col = rng.randint(1, cols - 1)
            newStateSeed = rng.randint(0, 2**32)
            visit = states[row][col]
            prob = np.random.random()

            neighbors = getUniqueNeighbors(states, row, col)
            rng.seed(newStateSeed)
            newState = rng.choice(neighbors)
        
            moore = states[row-1 : row+2, col-1 : col+2]
            if newState not in moore: 
                raise AssertionError("newstate not in neighborhood")
            else:
                unlike = 0
                for i in range(len(moore)):
                    for j in range(len(moore)):
                        if i != 1 or j != 1:
                            neigh1 = 1 - delt(moore[i][j], newState)
                            neigh2 = 1 - delt(moore[i][j], visit)
                            neigh = neigh1 - neigh2
                            unlike += neigh
                A = delt(0, visit)
                B = delt(0, newState)
                liq = (B - A) * Fliq
                solid = ((1 - B) - (1 - A)) * Fsol
                energy = unlike + liq + solid
            
                if prob < metropolis0(energy):
                    states[row][col] = newState

        sol = solidCount(states)
        area[currTime + 1] = sol[0]
        if currTime % 5 == 0:
            gRow += 1
            gRow += 1
            noBuffer = states[1 : rows-1, 1 : cols-1]
            #grainArea[gRow] = grainCount(noBuffer, len(grains.keys()))[1:]

        currTime += 1
    disagreement = np.zeros(time + 1)
    for i in range(len(disagreement)):
        disagreement[i] = abs(area[i] - g[i])
    maxDiff = max(disagreement)
    
    if plot:
        #plt.imshow(states)
        t = np.linspace(0, time, time + 1)
        plt.scatter(t, area, color = "none", edgecolor = "blue", label = 'CA')
        plt.plot(t, g, 'k', label = 'JMAK')
        plt.show()
    print(maxDiff)
    return states, solidFrac, grainArea

mcpm_efficient()
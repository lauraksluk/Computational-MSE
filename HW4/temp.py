import numpy as np
import json 
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

#def getGrains(data):

def mcpm_slow(time, ins, plot = True, rng_seed = None):
    # with open(ins) as f:
    #     data = json.load(f)
    # rows = data["width"]
    # cols = data["length"]
    # init = data["initial"]
    rows = 4
    cols = 5
    init = ins
    rng = np.random.RandomState()
    
    solidFrac = np.zeros(time + 1)
    totalArea = (rows - 2) * (cols - 2)

    states = deepcopy(init)
    currTime = 0
    print(states)
    while currTime < time:
        rng.seed(rng_seed)
        rng_seed = rng.randint(0, 2**32)
        
        for i in range(totalArea):
            states = mcpm_step(states, (rows, cols), rng, metropolis0)
            print(states)

        currTime += 1
    #print(states)
    if plot:
        #plt.imshow(states)
        t = np.linspace(0, time, time + 1)
        plt.plot(t, solidFrac)
        plt.show()
    
    return states, solidFrac, grainArea

def test():
    m1 = [[0,0,0,0,0], [0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,0]]
    m = np.asarray(m1)
    return mcpm_slow(15, m, True, None)
test()

 


def mcpm_efficient(time = 15, ins = 'random.json', plot = True, rng_seed = None):
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
    (initSolid, grains) = solidCount(init)
    solidFrac[0] = float (initSolid) / totalArea

    everyFive = (time // 5) + 1
    grainArea = np.zeros((everyFive, len(grains.keys())))
    ind = 0
    gRow = 0
    for key in grains:
        grainArea[gRow][ind] = grains[key]
        ind += 1

    states = deepcopy(init)
    currTime = 0
    while currTime < time:
        rng.seed(rng_seed)
        rng_seed = rng.randint(0, 2**32)
        for i in range(totalArea):
            row = rng.randint(1, rows - 1)
            col = rng.randint(1, cols - 1)
            newStateSeed = rng.randint(0, 2**32)
            visit = states[row][col]
            prob = rng.random()

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
        solidFrac[currTime + 1] = float (sol[0]) / totalArea

        if currTime % 5 == 0:
            gRow += 1
            newG = sol[1]
            index = 0
            for k in newG:
                grainArea[gRow][index] = newG[k]
                index += 1

        currTime += 1
    
    if plot:
        #plt.imshow(states)
        t = np.linspace(0, time, time + 1)
        plt.plot(t, solidFrac)
        plt.show()
    
    return states, solidFrac, grainArea
    """
    Full MC- load data, run MC for desired number of iterations.
    
    This function has the same inputs, outputs, and assumptions as mcpm_slow().
    With the same inputs, it should produce the same outputs. 
    
    Implement the MC method *without* using the functions you have written above.
    Combining all the calculations into this single function should achieve a 
    moderate speedup.
    
    This function assumes the following:
        - metropolis p0 = 1
        - Fliq = 3
        - Fsol = 0
    
    Keep track of total solid fraction and 
    area of each  grain as the simulation progresses.
    When computing the solid fraction and grain areas, only consider
    pixels that are not in the buffer regions!
    ie the system:
                [[0,0,1,0,0,0],
                 [0,1,1,2,0,0],
                 [0,0,2,2,2,0]]
    has a solid fraction of 0.75 and grain areas of [2, 1].
    
    
    The solid fraction should be tracked starting with the initial state
    (before any MC flips) and then after each subsequent MC step. 
    
    For example, if the initial solid fraction is 0.01, and increases
    by 0.01 every time step, the after t=3 the solid fraction
    array should be: [0.01, 0.02, 0.03, 0.0.04]
                  t=   0     1     2    3   
    
    The area of each grain should be tracked after every 5 MC steps.
    The area for grain with id i should at time t should be accessible by
    indexing grain_areas[t, i-1] (ie we are not tracking the liquid state)
    For example, consider a system with two initial grains [states 1 and 2].
    Grain 1 has initial area 1, and increases by 1 after every step.
    Grain 2 has initial area 2, and increases by 3 after every step.
    
    The array after 15 iterations shold be as follows:
                          [[ 1,  2],
                           [ 6, 17],
                           [11, 32],
                           [16, 47]]
    
    Note that 1 MC step consists of 1 iteration for every position,
    excluding the buffer region.Thus, for a 6 x 8 system,
    (6-2)*(8-2)= 24 random sites should be visited.
    
    To control the randomness to allow for autograding, create a np.random.RandomState instance
    (referred to as rng here). Then, inside the inner MC loop, the first two lines of code should
    be as follows:
        1) seed rng with rng_seed
        2) use rng.randint() to select an integer seed from the interval [0, 2**32) to be used as the seed
           for step 1) during the next iteration (ie the new value of rng_seed)
        3) use rng.randint() to select the index of the row to visit
        4) use rng.randint() to select the index of the column to visit
        5) use rng.randint() to select an integer seed from the interval [0, 2**32) to be used as
           the seed when picking the new state.
    
    Once again, the new states should be selected from a sorted array of unique elements in the neighborhood
    of a given state. To choose the new state:
        6) seed rng with the seed from step 5)
        7) use rng.choice() to pick the new state
        
    Parameters
    ----------
    time: int
        time (in monte carlo steps) for simulation to run
        
    ins: str or Path object
        path to read inputs from
        
    plot: Bool
        if True, display the final state as an image or heatmap
        and plot the solid fraction fmc as a function of time.
    
    rng_seed: int or None
        initial random seed, on the interval [0, 2**32). If None, a seed will
        be chosen randomly.
    
    Returns
    ----------
    state: ndarray
        array of final state of system

    fmc: ndarray
        solid fraction (format described in instructions above)
    
    grain_areas: ndarray
        area of each grain (format described in instructions above)
    """
    pass

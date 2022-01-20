# Laura (Kai Sze) Luk (kluk@andrew.cmu.edu)

import json, random, math
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt

# Returns the state of the majority of the solid neighbors
def getMajority(up, down, left, right):
    arr = [up, down, left, right]
    arr1 = [up, down, left, right]
    for dir in arr:
        if dir == 0: 
            arr1.remove(dir)
    count = [0]*(len(arr1))
    for i in range(len(arr1)):
        count[i] = arr1.count(arr1[i])
    index = count.index(max(count))
    return arr1[index]

# Calculates the initial solid area fraction 
def areaFractionJ(data):
    count = 0
    for i in range(len(data)):
        for j in range(len(data[0])):
            if data[i][j] > 0:
                count += 1
    return count

# Main method to perform a 2d Cellular Automaton solidification simulation.
def twodsol(timesteps = 10, prob = 1.0, inpath = 'single3.json', out1 = 'grains.png', out2 = 'jmak.png', plot = True):
    #load input deck
    with open(inpath) as f:
        data = json.load(f)
    rows = data["width"]
    cols = data["length"]
    init = data["initial"]

    totalArea = (rows - 2) * (cols - 2)  #exclude buffer 
    area = np.zeros(timesteps + 1)   #fCA initialize
    g = np.zeros(timesteps + 1)      #fJMAK initialize

    initialN = areaFractionJ(init)       #initial area fraction
    initialSolid = areaFractionJ(init)
    a0 = 5

    #lambda expression for JMAK area fraction equation
    jmak = lambda t : 1 - math.exp(float(-initialN) * prob * (2 * t ** 2 + 2 * t + 1) / totalArea)
    growth = lambda t : math.pi * (t**2) + a0
    #make copy of init
    states = deepcopy(init)
    #iterate thru time
    time = 0
    while (time <= timesteps):
        g[time] = growth(time)
        area[time] = initialSolid 
        final = deepcopy(states)
        for i in range(1, rows - 1):
            for j in range(1, cols - 1):
                if (states[i][j] == 0):
                    up = states[i-1][j]
                    down = states[i+1][j]
                    left = states[i][j-1]
                    right = states[i][j+1]
                    if (up > 0 or down > 0 or left > 0 or right > 0):
                        x = random.random()
                        if (x < prob):
                            final[i][j] = getMajority(up, down, left, right)
                            initialSolid += 1
        states = final                    
        time += 1
    #get fCA and fJMAK disagreement
    disagreement = np.zeros(timesteps + 1)
    for i in range(len(disagreement)):
        disagreement[i] = abs(area[i] - g[i])
    maxDiff = max(disagreement)
    
    #plot and display/save results
    finalArr = np.array(final)
    if plot:
        #plt.imshow(finalArr)
        #plt.imsave(out1, finalArr)

        t = np.linspace(0, timesteps, timesteps + 1)
        plt.scatter(t, area, color = "none", edgecolor = "blue", label = 'CA')
        plt.plot(t, g, 'k', label = 'JMAK')
        # plt.xlabel('time')
        # plt.ylabel('area fraction')
        # plt.legend()
        plt.show()
        #plt.savefig(out2)
    print(maxDiff)
    return finalArr, areaFrac, fJMAK, maxDiff


twodsol()

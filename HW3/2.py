import json, random, math
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt

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

def areaFractionJ(data):
    count = 0
    for i in range(len(data)):
        for j in range(len(data[0])):
            if data[i][j] > 0:
                count += 1
    return count

def getStop(disagree, percent):
    for i in range(len(disagree)):
        if disagree[i] > percent:
            return i

def twodsol(timesteps = 20, prob = 0.8, inpath = 'random2.json', out1 = 'grains.png', out2 = 'jmak.png', plot = False):
    
    with open(inpath) as f:
        data = json.load(f)
    rows = data["width"]
    cols = data["length"]
    init = data["initial"]

    totalArea = (rows - 2) * (cols - 2)  #exclude buffer 
    areaFrac = np.zeros(timesteps + 1)   #fCA initialize
    fJMAK = np.zeros(timesteps + 1)      #fJMAK initialize

    initialN = areaFractionJ(init)       #initial area fraction
    initialSolid = areaFractionJ(init)

    
    jmak = lambda t : 1 - math.exp(float(-initialN) * prob * (2 * t ** 2 + 2 * t + 1) / totalArea)

    states = deepcopy(init)

    time = 0
    while (time <= timesteps):
        fJMAK[time] = jmak(time)
        areaFrac[time] = float(initialSolid) / totalArea
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
        disagreement[i] = abs(areaFrac[i] - fJMAK[i])
    maxDiff = max(disagreement)
    
    #plot and display/save results
    finalArr = np.array(final)
    if plot:
        #plt.imshow(finalArr)
        #plt.imsave(out1, finalArr)

        t = np.linspace(0, timesteps, timesteps + 1)
        plt.scatter(t, areaFrac, color = "none", edgecolor = "blue", label = 'CA')
        plt.plot(t, fJMAK, 'k', label = 'JMAK')
        plt.xlabel('time')
        plt.ylabel('area fraction')
        plt.legend()
        #plt.show()
        #plt.savefig(out2)
    
    return finalArr, areaFrac, fJMAK, maxDiff
#twodsol()

def test1():
    allDiff3 = [0]*10
    for i in range(10):
        (_, _, _, diff) = twodsol(25, 0.5, 'random2.json')
        allDiff3[i] = diff
    return float(sum(allDiff3)) / 10
#print(test1())

def test():
    allDiff3 = [0]*10
    for i in range(10):
        (_, _, _, diff) = twodsol(30, 0.3, 'random2.json')
        allDiff3[i] = diff
    p3 = float(sum(allDiff3)) / 10

    allDiff5 = [0]*10
    for i in range(10):
        (_, _, _, diff5) = twodsol(25, 0.5, 'random2.json')
        allDiff5[i] = diff5
    p5 = float(sum(allDiff5)) / 10

    allDiff8 = [0]*10
    for i in range(10):
        (_, _, _, diff1) = twodsol(20, 0.8, 'random2.json')
        allDiff8[i] = diff1
    p8 = float(sum(allDiff8)) / 10
    
    dis = np.array([p3, p5, p8])
    p = np.array([0.3, 0.5, 0.8])
    plt.scatter(p, dis, color = "none", edgecolor = "blue")
    plt.xlabel('attachment probability')
    plt.ylabel('average maximum disagreement')
    plt.show()

print(test())

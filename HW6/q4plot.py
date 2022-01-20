import numpy as np
import matplotlib.pyplot as plt

def plot4():
    values = np.loadtxt("LJthermoinfo.txt")
    time = values[:,0]

    temp = values[:,1]
    pressure = values[:,3]
    density = values[:,5]
    pe = values[:,6]
    ke = values[:,7]
    E = values[:,8]

    plt.plot(time, temp, label="T")
    plt.plot(time, pressure, label="P")
    plt.plot(time, density, label="density")
    plt.plot(time, pe, label="U")
    plt.plot(time, ke, label="K")
    plt.plot(time, E, label="E")

    plt.legend()

    plt.show()

    

plot4()
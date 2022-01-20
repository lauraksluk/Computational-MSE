# Laura (Kai Sze) Luk (kluk@andrew.cmu.edu)

#!/usr/bin/python
import json
import numpy as np 
from matplotlib import pyplot as plt

def ftcs(inpath = 'sample_in.json', outpath = 'ftcs_plot.png', plot = True):
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
        plt.savefig(outpath)

    return x, tmaxC
                    
    """
    Perform the following steps.
    1) Load data from json file. Data will be a dictionary with the following key-value pairs:
        'xmax': <length of domain in x direction>
        'delx': <size of x increment>
        'tmax': <total time>
        'delt': <size of time increment>
        'alpha': <diffusion coefficient>
        'c0': <constant composition boundary condition at x=0>
        'cxmax': <constant composition boundary condition at x=xmax>
        'ic': <initial condition, as a lambda function>
    2) Set initial and boundary conditions
    3) Perform FTCS finite differencing for time = 0 to tmax, using delt as the timestep.
       Be sure to enforce the boundary conditions.
    4) If plot=True, generate and save as PNG a plot with the following format:
        plot initial concentration profile (x,c(x,t=0)) as a dashed black line
        plot final concentration profile (x, c(x,t=tmax)) as a solid blue line
        label the x axis x and the y axis c(x,t)
        create a legend to distinguish between the curves for c(x,0) and c(x,tmax)
    5) At the very end, use the Python 'return' command to return the x array and
       the final concentration profile c(x,tmax) array
    
    Parameters
    ----------
    inpath, outpath: str or Path object
        path to input deck (json) and output concentration profile plot (png)
    plot: Bool
        if True, results are plotted and saved in outpath

    Input Deck
    ----------
    xmax, delx, tmax, delt, alpha,c0,cxmax: float
        length of domain in x direction, size of x increment, total time, size of time increment,
        diffusion coefficient, composition boundary condition at x=0, composition boundary condition at x=xmax
    ic: str
        initial condition, expressed as a lambda function
        note: to convert ic from string to lambda function, can use ic=eval(inputdeck["ic"])
    
        
    Returns
    -------
    x: array
        vector of x values
    ctmax: array
        vector of concentration values at t=tmax
    
    Saves
    ------
    ftcs_plot.png: plot
        Plot of results saved as a .png file
    """

    pass
ftcs()



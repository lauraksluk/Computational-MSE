import numpy as np
import json 
import matplotlib.pyplot as plt

def pfm(alpha=5.,delta=1.,gamma=1.,delTemp=-10.,delt=0.01,inpath='random4.json',tmax=0.3,plot=True):
    """
    This function performs a 2D, non-conserved (Allen-Cahn) Phase Field simulation 
    nominally corresponding to solidification in a pure material. 
    A model double well potential is used for fv.
    The solver is FTCS finite differencing.

    Note that only solid and liquid phases are tracked; individual grains are not resolved.

    Inputs:
       tmax = number of finite differencing time steps
       delt = size of finite differencing time step
       input = initial struture file name (json file)
       delTemp = undercooling, i.e. T - Tmelt
       gamma = volumetric energy scaling constant 2
       delta = volumetric energy scaling constant 1
       alpha = interfacial energy scaling constant

    Outputs:
       an image of the evolved order parameters
       a plot comparing the results to the JMAK model
    """
    #read in initial structure
    with open(inpath) as f:
         ipd = json.load(f)
    width = ipd["width"]
    length = ipd["length"]
    initial = ipd["initial"]
    
    #precondition the initial system
    phis = np.array(initial)
    phis = phis.astype(float)
    phis[phis <= 0] = -1.
    phis[phis > 0] = 1.
    
    #initialize variables
    tsteps = np.int((tmax/delt))
    delx = 1
    dely = 1
    area = (width - 2)*(length - 2)
    L = 1.
    solarea = np.zeros(tsteps)
    solarea[0] = np.count_nonzero(initial)
    newphis = phis
    
    #evolve via phase field model
    time = 0
    for time in range(tsteps):
        for i in range(1,width-1):
            for j in range(1, length-1):
                
                #find current cell and its neighbors on the sq(1) (von Neumann) lattice
                cell = phis[i,j]
                right = phis[i,j+1]
                left = phis[i,j-1]
                above = phis[i-1,j]
                below = phis[i+1,j]
                
                #perform an FTCS FDM step                
                dfv1 = 4*delta*(cell**3 - cell)
                dfv2 = (15*gamma/8)*(1 - 2*cell**2 + cell**4)*delTemp
                dfix = (alpha/delx**2)*(left - 2*cell + right)
                dfiy = (alpha/dely**2)*(above - 2*cell + below)
                newphis[i,j] = -1*L*delt*(dfv1 + dfv2 - dfix -dfiy) + cell
        phis=newphis
        solarea[time] = (phis > 0).sum()
        
    if plot:
        plt.imshow(phis)
        plt.show()
        
        solarea=solarea/area
        t=np.arange(tsteps)
        single=(np.pi*t**2)/(area)
        nuc=np.unique(initial).size - 1
        jmak=(1-np.exp(-1*nuc*single))
        plt.plot(t,solarea,'bo',t,jmak,'k-',mfc='none')
        plt.ylabel('solid fraction')
        plt.xlabel('t')
        
        plt.show()

    return solarea

if __name__ == "__main__":
    pfm()

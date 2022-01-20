import numpy as np
import numpy.linalg as la
import json 
import matplotlib.pyplot as plt

def philj(epsilon, sigma, rc, rij):
    #if (rij == rc): return -epsilon
    #if (rij == sigma): return 0
    if (rc == 0): return 0
    if (rij >= rc): return 0
    if (rij < rc):
        first = np.power((sigma / rij), 12) 
        second = np.power((sigma / rij), 6) 
        return 4 * epsilon * (first - second)


def latticesum(a = 1.7, epsilon = 0.9, sigma = 1.1, rc = 1.01, inpath = 'atoms.npy'):
    rc = rc*a
    atomCoord = np.load(inpath)
    

    """
    Reads in a set of atomic coordinates, sums pair potential interactions,
    and returns the system potential energy.

    Perform the following steps.
    1) Load input data from npy file. Data will be triplets of normalized atomic coordinates ri = (xi,yi,zi).
        An FCC lattice is assumed.
    2) Determine the number of atoms and compute the number of unit cells in each dimension (Ncells).
    3) For each atom pair in the lattice sum, compute interatomic distance rij: 
        rij = a*Ncells*|rj - ri|  where |.| denotes the vector norm
    4) Call the function philj to compute the Lennard-Jones interatomic pair potential
    5) Return the sum of all interatomic pair potentials, i.e. the total system potential energy

    
    Parameters
    ----------
    a: float
        Lattice parameter.
    
    rc: float
        interaction cutoff radius. Note that the input rc is the proportionality constant between
        the cutoff radius and the lattice parameter. (Or rc is input in the units of lattice
        parameters.) To convert to a length, rc is multiplied by a.
    
    inpath: str
        input deck of normalized atomic coordinates.
    
    Returns
    -------
    latticesum: float
        total system potential energy
    
    """
    pass

def latticesumpbc(a=1.7,epsilon=0.9,sigma=1.1,rc=1.01,inpath='atoms.npy'):
    rc = rc*a
    """
    Identical to the above function, but includes periodic boundary conditions via
    the minimum image convention.
    """
    pass


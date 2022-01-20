# Laura (Kai Sze) Luk (kluk@andrew.cmu.edu)

import json
import numpy as np


def fem_err(inpath = 'fem_results.json'):
    # load and extract data from file
    with open(inpath) as f:
        data = json.load(f)
    x = data['x']
    c = data['comp']
    numNodes = len(x)

    # construct array of delta x_i's
    delX = np.zeros(numNodes - 1)
    for i in range(numNodes - 1):
        delX[i] = x[i + 1] - x[i]
    # check delta x_i sums to 1
    assert(sum(delX) == 1)

    # equation for exact concentration
    exact = lambda x : x - x ** 2

    # array of errors 
    error = np.zeros(numNodes - 1)
    for i in range(len(error)):
        mid = x[i] + (delX[i] / 2)
        exactVal = exact(mid)

        femVal = 0.5 * (exact(x[i]) + exact(x[i + 1]))
        error[i] = abs(exactVal - femVal)
    
    # average error
    numElem = len(error)
    average = sum(error) / numElem
 
    return numElem, average

    """
    Perform the following steps.
    1) Load data from json file. Data will be a dictionary with the following key-value pairs:
        'x': list of xi values
        'comp': list of computed concentration values at x=xi
    2) Use the formula in your homework assignment to compute the error at each element center
    3) Caculate the average error by summing all element errors and dividing by the number of elements
    4) Use the Python 'return' command to return the number of elements N and the average error Eavg
    
    Parameters
    ----------
    inpath: str or Path object
        path to input deck (json)   
        
    Returns
    -------
    N: int
        number of elements
    Eavg: float
        average error
        
    Saves
    ------
    None
    
    """

    pass
#fem_err()


# Laura (Kai Sze) Luk (kluk@andrew.cmu.edu)

import random, json
import numpy as np
from PIL import Image as im
from matplotlib import pyplot as plt
from helper import *

def q1 (inpath = 'data.npy', outpath = 'q1_results.json'):
    data = np.load(inpath)
    numRows = len(data)
    numCols = len(data[0])
    random.seed(numCols)

    # get random array
    random_ints = np.array([0]*10)
    for i in range (10):
        random_ints[i] = random.randint(0, 20)
    
    # pad matrix
    diff = numRows - numCols
    temp = data.copy()

    if diff < 0: #more columns
        diff *= -1
        while diff != 0:
            if diff > numRows:
                temp = np.append(temp, data, 0)
                diff -= numRows
            else:
                add = data[:diff]
                temp = np.append(temp, add, 0)
                diff = 0
        padded = temp
    
    elif diff > 0: #more rows
        while diff != 0:
            if diff > numCols:
                temp = np.append(temp, data, 1)
                diff -= numCols
            else:
                add = data[:,:diff]
                temp = np.append(temp, add, 1)
                diff = 0
        padded = temp
    
    # get negative elements
    count = 0
    for i in range(numRows):
        for j in range(numCols):
            if padded[i, j] < 0:
                count += 1
                padded[i, j] = 0
    n_negative = count

    # get nonneg. elements
    count = 0
    for i in range(numRows):
        for j in range(numCols):
            if padded[i, j] != 0: count += 1
    n_nonzero = count

    # get count
    count = 0
    total = numRows * numCols
    while count <= total:
        randNum = random.randint(0, numRows - 1)
        selectedRow = padded[randNum]
        for i in range(numCols):
            if padded[randNum,i] != 0: count += 1
    
    # pack results
    result = dict()
    result['random_ints'] = random_ints
    result['n_negative'] = n_negative
    result['n_nonzero'] = n_nonzero
    result['count'] = count

    with open(outpath, 'w') as fp:
        json.dump(result, fp)

    return
    """
    Perform the following steps.
    1) Load data from npy file (will always be a 2d array)
    2) Seed the random number generator with a value equal to the number of columns in the array
    3) Use the seeded random number generator to generate an array of 10 random integers between 0 and 20 
    (inclusive- should include 0 and 20.) Store this array as a variable.
    4) Pad the matrix by adding extra rows or columns until the total number of rows and columns are equal. 
       Add extra rows or colmuns AFTER (ie extra rows below, and columns to the right of) existing ones.
       The values of the padded columns should 'wrap' around the array. For example, if you're adding an 
       extra 3 rows, their values should match the first 3 rows in the original array.
    5) Determine the number of negative values in the matrix. Store this in a variable.
    6) Change any negative values in the array to 0
    7) Store the total number of nonzero elements in the padded array as a variable.
    8) Run a while loop. In each step, pick a random integer between 0 and the total number of rows in the matrix.
       Select the row corresponding to that index, and add the total number of nonzero elements in that row to a 
       counter. When the counter exceeds the total number of elements in the matrix, stop the loop.
    9) Pack the results into a dictionary with the following key-value pairs:
        'random_ints': <result from step 3>
        'n_negative': <result from step 5>
        'n_nonzero': <result from step 7>
         'count': <result from step 8>
    10) Save the dictionary as a json file.
    
    Parameters
    ----------
    inpath, outpath: str or Path object
         Path to input file (npy) to load, and output file (json) to save data to
    
    Returns
    -------
    None
    
    Saves
    ------
    q1_results: dict
        json dictionary with format described in step 10) above
    """
    pass


def q2(rows = 10, cols = 11, outpath = 'q2_results.npy', display = True):
    # initialize matrix
    matrix = np.zeros(shape = (rows, cols), dtype = np.int64)

    center = cols // 2
    matrix[0, center] = 1

    # use helper function
    for i in range(1, rows):
        for j in range(1, cols - 1):
            selectedArr = matrix[i - 1, (j - 1):(j + 2)]
            count = 0
            for k in range(3):
                if selectedArr[k] != 0: count += 1
            if (count != 0 and count != 3):
                matrix[i, j] = helper(selectedArr)
    
    # display array as image
    if display:
        displayArr = im.fromarray(matrix, 'RGB')
        img.show()

    # save array in file
    np.save(outpath, matrix)

    return

    """
    Perform the following steps:
    0) import the helper function provided in helper.py (do not just copy/paste the function
       into your code)
    1) initialize a rows x cols matrix of zeros (int or Bool)
    2) Set the center coordinate of the first row to 1
       (note, to deal with even number of columns, use the index equal to
       half the number of columns. for example, for arrays with both 10 and 11 columns
       the 5th element will be considered to be the center coordinate)
    3) Wtarting with the 2nd row, write a nested loop that does the following:
        a) Loop through all rows (excluding the first one).
        b) Within each row, loop through all the columns, excluding the first and last column.
            For example, if the matrix is [[0,0,1,0,0]
                                           [0,1,2,3,4],
                                           [5,6,7,8,9]]
            Then you should visit positions with values 1,2,3, 6,7,8.
        a) for a given element in x at index [r, c], select the 3 values above it as an array
           (ie [x(r-1, c-1), x(r-1, c), x(r-1, c+1)]). In the example above, the values corresponding
           to the element whose value is 2 should be [0,1,0].
           Note that to do this, you cannot look at elements in the first or last column. The inner
           loop should start at the second column and end at the second-to last column.
        b) count the number of nonzero elements in this set. If it is 0, or 3, leave the value as 0.
           Otherwise, set the element to be the output of the helper function. 
           (to do this, call helper function with the array containing the 3 values as its input)
    4) if the input 'display' is set to True, display the results as an image or heatmap
    5) save the final matrix as a .npy file
    
    Parameters
    ----------
    rows, cols: int
        specifies the number of rows and columns in the array
    
    outpath: str or Path object
        path to save final array file to (as a .npy file)
    
    display: Bool
        if True, the array is displayed as an image. Otherwise, no display is created       
    
    Returns
    -------
    None
    
    Saves
    -------
    q2_results: ndarray
        Matrix with pattern filled in saved as a .npy file.
    """
    pass


def q3(x, y = lambda x: x**2, outfile = 'q3_output.json', plot = True):
    numVals = len(x)

    # obtain y values by evaluating
    yvals = np.empty(numVals)
    for i in range(numVals):
        yvals[i] = y(x[i])
    
    # obtain derivative values
    deriv = np.empty(numVals - 1)
    for j in range(1, numVals):
        ydiff = (yvals[j] - yvals[j - 1])
        xdiff = (x[j] - x[j - 1])
        slope = ydiff / xdiff
        deriv[j - 1] = slope
    
    # plot results
    if plot:
        plt.plot(x, yvals, 'k', label = 'y(x)')
        x1 = x[:numVals - 1]
        plt.plot(x1, deriv, 'b', linestyle = 'dashed', label = 'dydx')
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.legend()
        plt.show()
    
    # save results to json file
    result = dict()
    result['x'] = x
    result['y'] = yvals
    result['dydx'] = deriv

    with open(outpath, 'w') as fp:
        json.dump(result, fp)
    
    return

    """
    Perform the following steps:
    1) evaluate function y with input x.
    2) compute the numerical derivative of dy(x)/dx using a forward difference method:
        dy/dx = (y2-y1)/(x2-x1).
        The derivative is calculated for every point along x except the last point 
        (ie ignore the edge case.)
    3) if plot is true, generate a plot with the following format:
        plot (x,y(x)) as a solid black line
        plot (x, dy(x)/dx) as a dashed blue line
        label the x axis x and the y axis f(x)
        create a legend to distinguish between y(x) and dydx
    4) save the results to outfile (json format) with the following key-value pairs:
       'x': list of x values
       'y': list of y(x) values
       'dydx': list of computed derivative values
    
    Parameters
    ----------
    x: ndarray
        array of x values to evaluate y at (in order)
    y: callable
        function which takes x as input and returns y values, ie yvals = y(x)
    outfile: str or Path object:
        path to save outputs to (json)
    plot: Bool
        if True, results are plotted in the manner described in step 3)
       
    Returns
    -------
    None
    
    Saves
    ------
    q3_results: dict
        json-formatted dictionary with format described above in step 4). 
    
    """
    pass

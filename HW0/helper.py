import numpy as xx
def helper(x):
    validate()
    x = x.astype(xx.bool)
    if x.sum() == 2:
        return 1 if xx.all(x[1:]) else 0
    elif x.sum() == 1:
        return 1
    else:
        raise Exception('helper should only be called when 1 or 2 elements are True!')

def validate():
    if __name__ == "__main__":
        raise Exception('This function needs to be imported, not copied into the main file!')

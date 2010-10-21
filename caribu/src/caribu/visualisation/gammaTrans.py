from math import *

def gammaTrans(values, gamma=1):
    '''    return value normalised and raised at exponent gamma
    '''
    m = min(values)
    M = max(values)
    if (m == M) :
        norm = 1
    else:
        norm = M - m
    res = map(lambda(x): ((x - m) / float(norm))**gamma,values) 
    # write the node code here.

    # return outputs
    return res,

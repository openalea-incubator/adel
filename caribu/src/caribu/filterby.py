from itertools import *

def filterby(indices, values, condition):
    '''    
    Return values whose indices match condition
    '''
    #print 'condition ', condition
    index_value = ifilter(lambda x: condition(x[0]), izip(indices, values))
    res = [value for index, value in index_value]

    # return outputs
    return res, 

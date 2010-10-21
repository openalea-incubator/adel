def PARaggregators(caribu_outdict):
    '''    returns a dict of aggregators (0/1) for summing Eabs at different levels
    '''

    aggregators = { 'Total' : [1] * len(caribu_outdict['Plt']),
		    'Green' : map(lambda(x) : x is 1, caribu_outdict['Opt']),
		    'Stem' :  map(lambda(x) : x[0] is 0 and x[1] is 1, zip(caribu_outdict['Opak'],caribu_outdict['Opt'])),
		    'Leaves' :  map(lambda(x) : x[0] is not 0 and x[1] is 1, zip(caribu_outdict['Opak'],caribu_outdict['Opt'])),
		    'Soil' :  map(lambda(x) : x == 0, caribu_outdict['Opt'])}

    # write the node code here.


    # return outputs
    return aggregators

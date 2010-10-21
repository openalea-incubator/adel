

def bindDict(d,e):
    vals=[]
    for k in d.keys():
        if (isinstance(d[k],list)) : 
            nl = d[k]
        else:
            nl = [d[k]]
        nl.append(e[k])
        vals.append(nl)


    return dict(zip(d.keys(),vals))



def reduceDict(dictlist):
    '''    apply dictionary reduction 
    '''
    reduceddict = reduce(bindDict,dictlist)
    return reduceddict

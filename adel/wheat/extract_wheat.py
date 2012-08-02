import rpy2.robjects as robj
from rpy2.robjects.numpy2ri import numpy2ri
r = robj.r

import numpy as np
from numpy import transpose

def RlistAsDict(Rlist):
    """returns a dictionary containing the elements of the Rlist"""
    if r['is.null'](Rlist.names)[0]:
        Rlist.names = r['seq'](Rlist)
    return dict(zip([n for n in r.names(Rlist)],[obj for obj in Rlist]))

def extract_leaf_info(rdata_xy, rdata_sr):
    """
    Extract geometric information for each leaf database structure:
    genotype, plant number, leaf rank, data
    database key: rank, value : x,y,s,r
    """
#    global r
    #xy = r.load(rdata_xy)[0]
    #sr = r.load(rdata_sr)[0]
    xy = r.load(rdata_xy)[0]
    sr = r.load(rdata_sr)[0]
    # rename Rlists with list index, to ensure leaf retrieval by index as key in python	
    r('names(%s) = seq(%s)'%(xy,xy))
    r('names(%s) = seq(%s)'%(sr,sr))

    xy = RlistAsDict(r(xy))
    sr = RlistAsDict(r(sr))

    rank = xy.keys()
    rank.sort()

    # Assert that the two databases are compatible.
    rank_bis = sr.keys()
    rank_bis.sort()
    assert rank == rank_bis

    leaves = {}
    for k in rank:
        xyk = RlistAsDict(xy[k])
        leaves_id = xyk.keys()
        try: 
            s, radius = transpose(np.array(sr[k]))
        except ValueError:#sr is a dataframe
            s, radius = np.array(sr[k][0]),np.array(sr[k][1])
        for leaf_id in leaves_id:
            try:
                x, y = transpose(np.array(xyk[leaf_id])) #xy is a matrix
            except ValueError:#xy is a dataframe
                x, y = np.array(xyk[leaf_id][0]),np.array(xyk[leaf_id][1])
            leaves.setdefault(k,[]).append((x,y,s,radius))

    return leaves

def test():
    rdata_xy = 'data/So99.RData'
    rdata_sr = 'data/SRSo.RData'
    leaves = extract_leaf_info(rdata_xy, rdata_sr)
    rank = leaves.keys()
    rank.sort()
    assert rank == ['1', '2', '3', '4', '5', '6']
    return leaves


if __name__ == '__main__':
    leaves = test()
    import cPickle as Pickle
    f = open('leaves_wheat.db','w')
    Pickle.dump(leaves, f)
    f.close()





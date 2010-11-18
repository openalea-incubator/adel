import rpy2.robjects as robj
r = robj.r

from numpy import transpose

def RlistAsDict(Rlist):
    """returns a dictionary containing the elements of the Rlist"""
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
    xy = r.load(rdata_xy)
    sr = r.load(rdata_sr)
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
        s, radius = transpose(sr[k])
        for leaf_id in leaves_id:
            x, y = transpose(xyk[leaf_id])
            leaves.setdefault(k,[]).append((x,y,s,radius))

    return leaves

def test():
    rdata_xy = 'So99.RData'
    rdata_sr = 'SRSo.RData'
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





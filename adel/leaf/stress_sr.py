import numpy as np
from collections import OrderedDict

class Interval(object):
    def __init__(self, min_, max_):
        self._min = min_
        self._max = max_

    @property
    def min(self):
        return self._min

    @property
    def max(self):
        return self._max

    def __contains__(self, v):
        return self.min <= v <= self.max

    def __hash__(self):
        return hash((self.min, self.max))

    def __eq__(self, obj):
        return hash(obj) == hash((self.min, self.max))

    def __cmp__(self, obj):
        # incomplete implementation
        if isinstance(obj, Interval):
            if self.max <= obj.min:
                return -1
            elif self.min >= obj.max:
                return 1
            elif self == obj:
               return 0  
                
    def __repr__(self):
        return 'Interval('+str(self.min)+','+str(self.max)+')'

def stress_sr(s, r, d):
    '''    
    '''
    intervals = dict(zip( (Interval(m,M) for m, M in d), d.itervalues())) 

    keys = list(sorted(intervals))
    # add missing intervals with factor 1
    i = keys[0]
    if i.min > 0:
        intervals[Interval(0, i.min)]=1
    for index in range(len(keys)-1):
        i, j = keys[index:index+2]
        assert i < j
        if i.max < j.min:
            intervals[Interval(i.max, j.min)]=1


    keys = list(sorted(intervals))
    print keys
    i = 0
    interval = keys[i]
    factor = intervals[interval]
    p_min = interval.min
    # p_min
    l = []
    for p in s:
        if p not in interval:
            # accumulate values by taking account of the size of the previous interval modified by the factor
            p_min += (interval.max-interval.min)*factor
            i+=1
            interval = keys[i]
            factor = intervals[interval]
        l.append(p_min + (p-interval.min) *factor)

    print s
    print l
    new_s = np.array(l)
    p_max = new_s.max()

    new_s = new_s/p_max

    return new_s, r,

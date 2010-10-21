try:
    from math import isnan
except:
    #to be back compatile with python 2.5
    def isnan(num):
        return num != num

def optical_species(label):
    return int(label[:-11])

def plant_id(label):
    return int(label[-11:-6])

def transparency(label):
    return int(bool(leaf_id(label)))

def leaf_id(label):
    return int(label[-6:-3])

def elt_id(label):
    return int(label[-3:])

def nan_to_zero(x):
    return(0 if isnan(x) else x)

def vcaribuOut_adaptor(vcdict):
    '''    adaptor from vcaribu nrj dict to genout like dict
    '''
    d = vcdict[vcdict.keys()[0]]['data']
    for k in ('Eabs','Ei_inf','Ei_sup'):
        d[k] = map(nan_to_zero,d[k])
    eabs = [e * a for e,a in zip(d['Eabs'],d['area'])]
    opt=map(optical_species,d['label'])
    opak=map(leaf_id,d['label'])
    plt=map(plant_id,d['label'])
    elt=map(elt_id,d['label'])

    godict = {'Eabs':eabs, 'Area':d['area'],'Eabsm2':d['Eabs'],'EiInf':d['Ei_inf'],'EiSup':d['Ei_sup'],'Opt':opt,'Opak':opak,'Plt':plt,'Elt':elt} 
  
    return godict

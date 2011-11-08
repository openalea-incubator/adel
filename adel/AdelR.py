#import os
#from openalea.core.pkgmanager import PackageManager
#pm = PackageManager()
#pkg = pm.get('alinea.adel') 
#path = ''
#if pkg :
#    path = pkg.path
#path_ini = os.getcwd()

#
#import os
#import openalea.adel as my_package
#my_path = os.path.dirname(my_package.__file__) 

#
#create AdelR functions into rpy2 environment
from math import sqrt
import numpy
import rpy2.robjects as robj
from rpy2.robjects.numpy2ri import numpy2ri #force conversion mode of Robject
r = robj.r

try: 
    robj.globalEnv = robj.globalenv
except:
    pass
#content of rfiles using setuptools pkg_resources utility
from pkg_resources import resource_string

# Set Numeric Locale value 
r('Sys.setlocale(category="LC_NUMERIC",locale="C")')

#r.source(os.path.join(path,'Adel.R')) : code qui suit plutot a laisser au niveau des fonctions
rcode = resource_string(__name__, 'Adel.R')
r(rcode)
rcode = resource_string(__name__, 'setAdel.R')
r(rcode)
rcode = resource_string(__name__, 'UseAdel.R')
r(rcode)
rcode = resource_string(__name__, 'ArvalisToAdel.R')
r(rcode)
rcode = resource_string(__name__, 'genString.R')
r(rcode)
#r.source(os.path.join(path,'ArvalisToAdel.R'))
RrunAdel = robj.globalEnv['runAdel']
RsetCanopy = robj.globalEnv['setCanopy']
RsetAdel = robj.globalEnv['setAdeluser']
RdevCsv = robj.globalEnv['devTcsv']
RreadCsv = robj.globalEnv['readCsv']
RgenGeoAxe = robj.globalEnv['genGeoAxe']
RgenGeoLeaf = robj.globalEnv['genGeoLeaf']
RsetAdelArv = robj.globalEnv['setAdelArv']
RgenString = robj.globalEnv['genString']
RcanL2canS = robj.globalEnv['canL2canS']
#r.load('D:\Christian\Projets\BleMaladie\ConfrontationArvalis\Calage\.RData')


def readRData(fn):
    """return a dictionary containing the Robject of the RData file"""
    
    data = r.load(fn)
    return dict(zip([nom for nom in data],[r[nom] for nom in data]))

def saveRData(Robj,name,fn):
    """ name and write the Robj to RData file """
    r.assign(name,Robj)
    r.save(list = name,file = fn)
    return(name)

def RlistAsDict(Rlist):
    """returns a dictionary containing the elements of the Rlist"""
    return dict(zip([n for n in r.names(Rlist)],[obj for obj in Rlist]))

def dataframeAsdict(df):
    """ convert an RDataframe to a python dict """
    if r['is.null'](df)[0]:
        return None
    try:
        d = dict(zip( df.colnames, numpy.array(df)))
        return d
    except:
        d = dict([(k,numpy.array(df.r[k][0])) for k in r.colnames(df)])
        return d

def dataframe(d):
    """ convert a dict of numbers to an RDataframe  """
    df = {}
    if d is None:
        return r('as.null()')
    else:
        for k, v in d.iteritems():
            df[k] = numpy2ri(numpy.array(v))
    dataf = r['data.frame'](**df)
    return dataf

def Rdflist(dictOfdict):
    """ convert a dict of dict of numpy vectors into a Rlist of Rdataframe  """
    df_tags = dictOfdict.keys()
    df_value = dictOfdict.values()
    return r.list(**dict(zip(df_tags, [dataframe(dft) for dft in df_value])))

def RdflistAsdicts(Rdflist):
    """ convert a Rlist of Rdataframe into a dict of dict of numpy vectors """
    return dict(zip([n for n in r.names(Rdflist)],[dataframeAsdict(obj) for obj in Rdflist]))

def csvAsDict(fn,type=1):
    """ returns a dictionnary with the content of csv file as numpy vectors (one colum = one key) """
    df = RreadCsv(fn,type)
    return dict([(k,numpy.array(df.r[k][0])) for k in r.colnames(df)])

def setAdelArv(Rcalage,Rfunstr,np,sdlev = 20):
    """Creates a set of parameter for simulating np plants with adel from arvalis calage Rdata (deprecated)"""
    
    r(Rfunstr)
    RFun = robj.globalEnv['Rfun']
    p = RsetAdelArv(Rcalage,np,sdlev,RFun)
    return p

def setAdel(devT,RgeoLeaf,RgeoAxe,nplants = 1,seed = None, xydb = None, srdb = None):
    """Creates a set of parameter for simulating np plants with adel from R inputs (see adeldoc.R)"""
    if seed is None:
        rseed = r('as.null()')
    else:
        rseed = seed
        
    if xydb is None:
        rxydb = r('as.null()')
    else:
        rxydb = xydb
        
    if srdb is None:
        rsrdb = r('as.null()')
    else:
        rsrdb = srdb

    RdevT = Rdflist(devT)
    p = RsetAdel(RdevT,RgeoLeaf,RgeoAxe,nplants,rseed,rxydb,rsrdb)
    return p

def canL2canS(RcanT,srdb,shrink):
    res = RcanL2canS(RcanT,srdb,shrink)
    d = dataframeAsdict(res)
    return d

def setCanopy(RcanT,nplants = 1,randomize = True, seed = None):
    """ Duplicates or subsample a canopy table to create a new canopy with specified number of plants. Radomise allows for random sample and random azimuth"""
    if seed is None:
        rseed = r('as.null()')
    else:
        rseed = seed
        
    if randomize:
        rrand = 1
    else:
        rrand = 0
        
    can = RsetCanopy(RcanT,nplants,rrand,rseed)
    return can

def RunAdel(datesTT,plant_parameters,adelpars={'senescence_leaf_shrink' : 0.5,'startLeaf' : -0.4, 'endLeaf' : 1.6, 'stemLeaf' : 1.2,'epsillon' : 1e-6}):
    """ Run Adel model for each date in datesTT according to parameter list """
    
    if (type(datesTT) is not list):
        datesTT = [datesTT]
    x = robj.FloatVector(datesTT)
    ap = robj.r['list'](**adelpars)
    #chn = RrunAdel(x,plant_parameters,ap)
    #return [c[0] for c in chn]
    res = RrunAdel(x,plant_parameters,ap)
    d = dataframeAsdict(res[0])
    return d

def devCsv(axeTfn,dimTfn,phenTfn,earTfn,ssi2senTfn):
    """ Import development parameters for adel from csv files """
    args = (axeTfn, dimTfn, phenTfn)
    if earTfn is not None:
        args += (earTfn,)
    if ssi2senTfn is not None:
        args += (ssi2senTfn,)
    Rdat = RdevCsv(*args)
    return RdflistAsdicts(Rdat)

def genString(RcanopyT):
    """ Generate an Lsystem string from a R dataframe representing the canopy """
    chn = RgenString(RcanopyT)
    return chn[0]

def genGeoAxe(azM=None,daz=None,ibmM=None,dibm=None,incT=None,dinT=None,dep=None):
    """ generate geoAxe function for Adel """
    try:
        return RgenGeoAxe(azM,daz,ibmM,dibm,incT,dincT,dep)
    except:
        return RgenGeoAxe()

def genGeoLeaf(nlim=None,dazt=None,dazb=None):
    """ generate geoLeaf function for Adel """
    try:
        return RgenGeoLeaf(nlim,dazt,dazb)
    except:
        return RgenGeoLeaf()

def freeGeoAxe(rcode):
    """ returns geoAxe object from txt definition """

    r(rcode)

    return robj.globalEnv['geoAxe']

def freeGeoLeaf(rcode):
    """ returns geoLeaf object from txt definition """

    r(rcode)
    
    return robj.globalEnv['geoLeaf']



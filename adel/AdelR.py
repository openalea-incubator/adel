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
r = robj.r

from rpy2.robjects.numpy2ri import numpy2ri
try:
    numpy2ri.activate()#force auto-conversion mode of Robject to array
except:
    pass

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
#RgenGeoAxe = robj.globalEnv['genGeoAxe']
#RgenGeoLeaf = robj.globalEnv['genGeoLeaf']
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
    return(fn)

def RlistAsDict(Rlist):
    """returns a dictionary containing the elements of the Rlist"""
    return dict(zip([n for n in r.names(Rlist)],[obj for obj in Rlist]))

def _rvect_asarray(rvect):
    """"
    convert a r_vector into an array of numeric values or into an array of string if NA are present (numpy2ri replaces NA with numeric values otherwise) . 
    Will be deprecated when numpy will offer true NA type """

    if r['length'](r['which'](r['is.na'](rvect)))[0] > 0:
        return numpy.array(r['as.character'](rvect))
    else:
        return numpy.array(rvect)

def dataframeAsdict(df):
    """ convert an RDataframe to a python dict """
    if r['is.null'](df)[0]:
        return None
    try:
        d = dict(zip( df.colnames, numpy.array(df)))
        return d
    except:
        try:
            d = dict([(k,_rvect_asarray(df.r[k][0])) for k in r.colnames(df)])
        except:
            d = dict([(k,_rvect_asarray(df.rx2(k))) for k in r.colnames(df)])# r delegator is replaced by rx in new rpy2
        return d

def dataframe(d):
    """ convert a dict of numbers to an RDataframe  """
    df = {}
    if d is None:
        return r('as.null()')
    else:
        for k, v in d.iteritems():
            df[k] = r['as.numeric'](numpy2ri(numpy.array(v)))
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
    return dataframeAsdict(df)

def setAdelArv(Rcalage,Rfunstr,np,sdlev = 20):
    """Creates a set of parameter for simulating np plants with adel from arvalis calage Rdata (deprecated)"""
    
    r(Rfunstr)
    RFun = robj.globalEnv['Rfun']
    p = RsetAdelArv(Rcalage,np,sdlev,RFun)
    return p

def setAdel(devT,RcodegeoLeaf,RcodegeoAxe,nplants = 1,seed = None, xydb = None, srdb = None):
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
    r(RcodegeoAxe)
    geoAxe = robj.globalEnv['geoAxe']
    r(RcodegeoLeaf)
    geoLeaf = robj.globalEnv['geoLeaf']
    p = RsetAdel(RdevT,geoLeaf,geoAxe,nplants,rseed,rxydb,rsrdb)
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
    if len(res) <= 0:#empty canopy
        d = None
    else:
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

def genGeoAxe(azM=75,daz=5,ibmM=2,dibm=2,incT=60,dinT=5,dep=7):
    """ generate geoAxe R code for Adel """
    rcode = """
    geoAxe <- list(
    azT = function(a) {{
          ifelse(a == 0,
          runif(1) * 360,#plant azimuth
          {azTM:.2f} + (runif(1) - .5) * {dazT:.2f})
       }},
    incT = function(a) {{
       ifelse(a == 0,
              {incBmM:.2f} + (runif(1) - .5) * {dincBm:.2f},
              {incT:.2f} + (runif(1) - .5) * {dincT:.2f})
              }},
    dredT = function(a) {{
         #1.5 is an offset to avoid tiller superposed to mainstem
          ifelse(a == 0,
                 0,
                 1.5 + runif(1) * ({depMax:.2f}-1.5))
        }}
       )
       """
    return rcode.format(azTM = azM, dazT = daz, incBmM = ibmM, dincBm = dibm, incT = incT, dincT = dinT, depMax = dep)

def genGeoLeaf(nlim=4,dazt=60,dazb=10):
    """ generate geoLeaf function for Adel """
    rcode = """
    geoLeaf <- list(
     Azim = function(a,n,ntop) {{
            ifelse(ntop <= {ntoplim:d},
            180 + {dazTop:.2f} * (runif(1) - .5),
            180 + {dazBase:.2f} * (runif(1) - .5))
            }},
     Lindex = function(a,n,ntop) {{
              ntop + 1}}
              )
        """
    return rcode.format(ntoplim = nlim, dazTop = dazt, dazBase = dazb)

def freeGeoAxe(rcode):
    """ returns geoAxe code from txt definition """

    return rcode


def freeGeoLeaf(rcode):
    """ returns geoLeaf code from txt definition """

    return rcode




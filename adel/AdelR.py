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

import pandas
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
RcheckAxeDyn = robj.globalEnv['checkAxeDyn']
RgetAxeT = robj.globalEnv['getAxeT']
RgetPhenT = robj.globalEnv['getPhenT']
#RgetLeafT = robj.globalEnv['getLeafT']

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
        return numpy.array(r['ifelse'](r['is.na'](rvect),'NA',rvect))#as.character fails if first value is less than one character  (eg (1,NA,..))
    elif isinstance(rvect,robj.vectors.FactorVector):
        return numpy.array(r['as.character'](rvect))
    else:
        return numpy.array(rvect)

def dataframeAsdict(df):
    """ convert an RDataframe to a python dict """
    if r['is.null'](df)[0]:
        return None
#    try:
#        d = dict(zip( df.colnames, numpy.array(df)))
#        return d
#    except:
    try:
        d = dict([(k,_rvect_asarray(df.r[k][0])) for k in r.colnames(df)])
    except:
        d = dict([(k,_rvect_asarray(df.rx2(k))) for k in r.colnames(df)])# r delegator is replaced by rx/rx2 in new rpy2
    return d

def _is_iterable(x):
    try:
        x = iter(x)
    except TypeError: 
        return False
    return True
    
    
def dataframe(d):
    """ convert a dict of numbers to an RDataframe  """
    df = {}
    if d is None:
        return r('as.null()')
    else:
        for k, v in d.iteritems():
            rval = numpy2ri(numpy.array(v))
            if not _is_iterable(v):
                v = [v]
            if 'NA' in v:
                df[k] = r['as.numeric'](rval)
            else :
                df[k] = rval
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

def setAdel(devT,RcodegeoLeaf,RcodegeoAxe,nplants = 1,seed = None, xydb = None, srdb = None,sample = 'random'):
    """Creates a set of parameter for simulating np plants with adel from R inputs (see adeldoc.R)"""
    if seed is None:
        rseed = r('as.null()')
    else:
        rseed = seed
        
    if xydb is None:
        rxydb = r('as.null()')
    else:
        rxydb = r.load(xydb)[0]
        rxydb = r(rxydb)
        
    if srdb is None:
        rsrdb = r('as.null()')
    else:
        rsrdb = r.load(srdb)[0]
        rsrdb = r(rsrdb)

    RdevT = Rdflist(devT)
    r(RcodegeoAxe)
    geoAxe = robj.globalEnv['geoAxe']
    r(RcodegeoLeaf)
    geoLeaf = robj.globalEnv['geoLeaf']
    p = RsetAdel(RdevT,geoLeaf,geoAxe,nplants,sample,rseed,rxydb,rsrdb)
    return p

def plantSample(setAdelPars):
    """return id of plants used by setAdel to setup the canopy """
    p = r.names(setAdelPars)
    p = r['as.numeric'](p)
    return _rvect_asarray(p)
    
def checkAxeDyn(setAdelPars, dates, plant_density = 1):
    """ Compute the number of axes present in the canopy at given dates"""
    if (type(dates) is not list):
        dates = [dates]
    d = robj.FloatVector(dates)
    df = RcheckAxeDyn(d, setAdelPars, plant_density)
    return pandas.DataFrame(dataframeAsdict(df))
   
def getAxeT(setAdelPars):
    """ return axeT table"""
    df = RgetAxeT(setAdelPars)
    return pandas.DataFrame(dataframeAsdict(df))
    
def getPhenT(setAdelPars, axe='MS'):
    """ return phenT table"""
    df = RgetPhenT(setAdelPars, axe=axe)
    return pandas.DataFrame(dataframeAsdict(df))
    
#def getLeafT(setAdelPars):
#    """ return a dict of leaves plantn_axetype_metamer : argument_dict_ for_leaf_sampler"""
#    def id_args(label):
#        plant, axe, n, nf = label.split('_')
#        id = '_'.join((plant, axe, n))
#        args = {'axe':axe, 'n': int(n), 'nf': int(nf)}
#        return (id,args)
#    return dict([id_args(s) for s in RgetLeafT(setAdelPars)])
   
def canL2canS(RcanT,srdb,shrink=1):
    """outputs are : 
       - Ll : length of the blade
       - Lv : visible (emerged) length of blade (green + senesced, rolled + unrolled)
       - Lr : rolled part of the blade
       - Lsen : length of the senescent part of the blade (hidden + visible)       
       - L_shape : Mature length of the blade used to compute blade shape
       - Lw_shape : Maximal width of the blade used to compute blade shape
       - Linc : relative inclination of the base of leaf blade at the top of leaf sheath (deg) 
       - Laz : relative azimuth of the leaf
       - Lsect : the number of sectors per leaf blade
       - Gl : length of the sheath (hidden + visible)
       - Gv : emerged length of the sheath
       - Gsen : senescent length of the sheath (hidden + visible)
       - Gd : apparent diameter of the sheath
       - Ginc : relative inclination of the sheath
       - El: length of the internode (hidden + visible)
       - Ev: emerged length of the internode 
       - Esen: senescent length of the internode (hidden + visible)
       - Ed: diameter of the internode
       - Einc : relative inclination of the internode
 
    """

    sr = r.load(srdb)[0]
    sr=r(sr)
    res = RcanL2canS(RcanT,sr,shrink)
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

def RunAdel(datesTT,plant_parameters,adelpars={'senescence_leaf_shrink' : 0.5,'startLeaf' : -0.4, 'endLeaf' : 1.6, 'endLeaf1': 1.6, 'stemLeaf' : 1.2,'epsillon' : 1e-6, 'HSstart_inclination_tiller': 1, 'rate_inclination_tiller': 30}):
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

def devCsv(axeTfn,dimTfn,phenTfn,earTfn=None,ssi2senTfn=None):
    """ Import development parameters for adel from csv files and/or pandas dataframes """
    args = [axeTfn, dimTfn, phenTfn]
    if earTfn is not None:
        args += [earTfn,]
    if ssi2senTfn is not None:
        args += [ssi2senTfn,]
    # use tempdir + csv to pass correctly NA from pandas to R, and to let RdevCsv make the columns renaming, if any
    if any([isinstance(f,pandas.DataFrame) for f in args]) :
        import tempfile, shutil
        tmp_dir = tempfile.mkdtemp()
        for i in range(len(args)):
            df = args[i]
            if isinstance(df, pandas.DataFrame):
                f = tmp_dir + '/args%d.csv'%(i)
                df.to_csv(f, na_rep='NA', index=False)
                args[i] = f
        Rdat = RdevCsv(*args)
        shutil.rmtree(tmp_dir)
    else :        
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
          ifelse(a == 'MS',
          runif(1) * 360,#plant azimuth
          {azTM:.2f} + (runif(1) - .5) * {dazT:.2f})
       }},
    incT = function(a) {{
       ifelse(a == 'MS',
              {incBmM:.2f} + (runif(1) - .5) * {dincBm:.2f},
              {incT:.2f} + (runif(1) - .5) * {dincT:.2f})
              }},
    dredT = function(a) {{
         #1.5 is an offset to avoid tiller superposed to mainstem
          ifelse(a == 'MS',
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
     Azim = function(a,n,nf) {{
            ntop = nf - n
            ifelse(ntop <= {ntoplim:d},
            180 + {dazTop:.2f} * (runif(1) - .5),
            180 + {dazBase:.2f} * (runif(1) - .5))
            }},
     Lindex = function(a,n,nf) {{
              nf - n + 1}}
              )
        """
    return rcode.format(ntoplim = nlim, dazTop = dazt, dazBase = dazb)

def freeGeoAxe(rcode):
    """ returns geoAxe code from txt definition """

    return rcode


def freeGeoLeaf(rcode):
    """ returns geoLeaf code from txt definition """

    return rcode




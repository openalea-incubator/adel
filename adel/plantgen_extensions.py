#
#       Adel.plantgen_extensions
#
#       Copyright 2006-2014 INRIA - CIRAD - INRA
#
#       File author(s): Christian Fournier <Christian.Fournier@supagro.inra.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################

"""
Extensions of plantgen module developed for reconstructing echap canopies
"""

import pandas
import numpy
import random
import operator
from itertools import chain

from alinea.adel.exception import *
from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data, plantgen2adel
import alinea.adel.plantgen.tools as tools
from alinea.adel.AdelR import devCsv




def pgen_tables(pgen):
    axeT_, dimT_, phenT_, phenT_abs, dimT_abs, dynT_, phenT_first, HS_GL_SSI_T, tilleringT, cardinalityT, config = gen_adel_input_data(**pgen)
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
   
    return {'adelT': (axeT, dimT, phenT), 'phenT_abs':phenT_abs, 'dimT_abs':dimT_abs, 'phenT_first':phenT_first, 'HS_GL_SSI_T':HS_GL_SSI_T, 'tilleringT':tilleringT, 'cardinalityT':cardinalityT, 'config':config}
    

    
def adelT_to_devT(pgen):
    """ Creates devT tables from plantgen dict
    """
    devT = devCsv(*adelT)
    return devT, 


def flat_list(nested_list):
    return list(chain.from_iterable(nested_list))
   
   
def _parent(axe):
    return '.'.join(axe.rsplit('.')[:-1])  
    
def modalities(nff):
    m1,m2 = int(nff), int(nff) + 1
    p = m1 + 1 - nff
    return {m1: p, m2: 1 - p}
    
def cardinalities(proba, n, parents = None):
    if parents is not None:#filter botanical impossibilities
        proba = {k:v for k,v in proba.iteritems() if _parent(k) in parents}
        proba = {k: v / sum(proba.values()) for k,v in proba.iteritems()} 
        card  = {k:min(int(v*n), parents.get(_parent(k),n)) for k,v in proba.iteritems()}
    else:
        card  = {k:int(v*n) for k,v in proba.iteritems()}
    if parents is not None:   #filter saturated parents after int rounding
        proba = {k:v for k,v in proba.iteritems() if (card[k] < parents.get(_parent(k),n) or _parent(k) == '')}
    
    missing = int(n - sum(card.values()))
    new=[1]
    while (missing > 0 and len(new) > 0):
        sorted_p = sorted(proba.iteritems(), key=operator.itemgetter(1), reverse=True)
        new = [sorted_p[i][0] for i in range(min(len(sorted_p),missing))]
        for k in new:
            card[k] += 1
            if parents is not None:
                if (card[k] >= parents.get(_parent(k),n) and _parent(k) != ''):
                    proba.pop(k)
        missing = int(n - sum(card.values()))   
            
    res = {k:v for k,v in card.iteritems() if v > 0}
    if parents is not None:
        parents.update(res)
    return res
 

 
def axis_list(emited_cohorts, theoretical_probabilities,  nplants = 2):
    """ compute cardinalities of axis in a stand of n plants
    The strategy used here is based on deterministic rounding, and differs from the one used in plantgen (basd on random sampling). Difference are expected for small plants numbers
    """
     
    df = emited_cohorts
    df = df.set_index('cohort') 
    cohort_cardinalities = {c:round(df.ix[c,'total_axis'] * nplants) for c in df.index}
    
    modal_proba = {c:modalities(df.ix[c,'nff']) for c in df.index}
    
    cohort_modalities = {k:cardinalities(modal_proba[k],v) for k,v in cohort_cardinalities.iteritems()}
    cohort_mods = {k:flat_list(map(lambda x: [x[0]] * int(x[1]),v.items())) for k, v in cohort_modalities.iteritems()}

    
    p = theoretical_probabilities
    axis_p = [(k[0],(k[1],v)) for k,v in p.iteritems()] 
    axis_proba = {k:dict([a[1] for a in axis_p if a[0] == k]) for k in dict(axis_p)}
    axis_proba = {k:{kk:vv/sum(v.values()) for kk,vv in v.iteritems()} for k,v in axis_proba.iteritems()}
    
    parents={'':0}
    cohort_axis = {k:cardinalities(axis_proba[k], sum(v.values()),parents) for k, v in cohort_modalities.iteritems()}
    cohort_ax = {k:flat_list(map(lambda x: [x[0]] * int(x[1]),v.items())) for k, v in cohort_axis.iteritems()}
    
    axis = {k:zip(cohort_mods[k],cohort_ax[k]) for k in cohort_mods}
    axis_list = [map(lambda x: (k,x[0],x[1]), v) for k,v in axis.iteritems()]
    
    return flat_list(axis_list)
    
    
def plant_list(axis, nplants = 2):

# TO DO for more robust name selection : test if any parent of same cohort is present on the plant to filter candidates,and, in update find the matching prent and then choose a compatible name : debug possible avec get_reconstruction Tremie 12

    def _choose_plant(axe_name, plantlist):   
        candidates = filter(lambda x: (axe_name not in x) and (_parent(axe_name) in x or _parent(axe_name) == ''), plantlist)
        if len(candidates) == 0:
            raise AdelImplementationError(' Unable to build plants from cardinalities of axis...')
        return random.sample(candidates,1)[0]
        
    def _update(plantlist, axe):
        plant = _choose_plant(axe[-1],plantlist)
        plant.update({axe[-1] : axe[:-1]})
        return plantlist
    
    plants = []
    for i in range(nplants):
        plants.append({})   
    plants= reduce(_update, axis, plants)
    
    return plants
    
def axeT_user(plants):

    iplant = flat_list([[i+1]*len(x) for i,x in enumerate(plants)])
    df = pandas.DataFrame(index=range(len(iplant)),
                                           columns=['id_plt', 'id_cohort', 'id_axis', 'N_phytomer_potential', 'N_phytomer', 'HS_final', 'TT_stop_axis', 'TT_del_axis', 'id_dim', 'id_phen', 'id_ear', 'TT_app_phytomer1', 'TT_col_phytomer1', 'TT_sen_phytomer1', 'TT_del_phytomer1'],
                                           dtype=float)

    df['id_plt'] = iplant
    df['id_axis'] = flat_list(map(lambda x: x.keys(),plants))
    df['id_cohort'] = flat_list(map(lambda x: zip(*x.values())[0],plants))
    df['N_phytomer_potential'] = flat_list(map(lambda x: zip(*x.values())[1],plants))
    df['id_phen'] = df['id_cohort'] * 100 + df['N_phytomer_potential']
    
    df= df.sort(['id_plt','id_cohort','id_axis'])
    return df
    
def dynT_user(MS_parameters = {'a_cohort':1. / 110.,'TT_hs_0':160.,'TT_hs_N_phytomer_potential':1100,'n0':4.5,'n1':2.5,'n2':5}, primary_tillers = ['T%d'%(i) for i in range (1,4)]):
    p = primary_tillers
    p.sort()
    idaxis = ['MS'] + p
    df = pandas.DataFrame(index=idaxis,
                          columns=['id_axis','a_cohort','TT_hs_0','TT_hs_N_phytomer_potential','n0','n1','n2'],
                          dtype=float)
    df.ix['MS'] = pandas.Series(MS_parameters)
    df['id_axis'] = idaxis
    df['TT_hs_N_phytomer_potential'] = df['TT_hs_N_phytomer_potential'][0]
    df = df.reset_index(drop=True)
    return df
 
class GL_model(object):
    """
    An object interface to a plantgen Green Leaf model variant
    This variant handles cases where GL=f(HS) is fixed whatever nff
    """

    def __init__(self, HS, GL, nff=12, n0=4.4, n1=1.5, n2=5, hs_t1=8):
        self.n0 = n0
        self.n1 = n1 
        self.hs_t1 = hs_t1
        self.HS = HS
        self.GL = GL
        self.nff = nff
        self.hs_t2 = nff
        self.n2 = n2
               
    def GL_number(self):
        GLpol = pandas.DataFrame({'HS':self.HS, 'GL':self.GL})
        GLpol = GLpol.ix[GLpol['HS'] > self.hs_t2,:]
        return GLpol
     
    def polyfit(self, a_start = -4e-9):
        GLpol = self.GL_number()
        c = (self.n2 - self.n1) / (self.hs_t2 - self.hs_t1) - 1
        fixed_coefs = [0.0, c, self.n2]
        a, rmse = tools.fit_poly((GLpol['HS'] - self.hs_t2), GLpol['GL'], fixed_coefs, a_start)
        return numpy.poly1d([a] + fixed_coefs), rmse
        
    def curve(self, step = 0.1, a_start = -4e-9):
        pol, rmse = self.polyfit(a_start) 
        lin = pandas.DataFrame({'HS':[0, self.n0, self.hs_t1, self.hs_t2],
                               'GL':[0, self.n0, self.n1, self.n2]})
        xpol = numpy.arange(self.hs_t2, 2 * self.nff, step)
        ypol = pol(xpol - self.hs_t2)
        dpol = pandas.DataFrame({'HS':xpol,'GL':ypol})
        dpol = dpol.ix[dpol['GL'] >= 0,:]
        return pandas.concat([lin, dpol])
       
def time_of_death(nplants, density_table, relative_density  = False):
    """
    return n times of death for an effective of nplants that should suit density time course given in density_data
    density data is a TT, density pandas DataFrame
    """
    df = density_table.sort('TT')
    if relative_density:
        card = df['density'] * nplants
    else:
        card = df['density'] * 1. / df['density'].iloc[0] * nplants
    ndead = card.iloc[0] - round(min(card))
    tdeath = numpy.interp(range(int(card.iloc[0] - ndead), int(card.iloc[0])),card[::-1],df['TT'][::-1])
    return tdeath

def kill_axis(devT, who, when, TT_stop_del = 2.8 * 110):
    """
    update devT tables by killing axis at pre-defined time
    who is a list of boolean mask for rows of devT:axeT table identifying the axis that should die
    when is a prralell list of time of death for axis dentified by who
    """
    df = pandas.DataFrame(devT['axeT'])
    dp = pandas.DataFrame(devT['phenT'])
    tdel = pandas.Series([float(v) if v != 'NA' else 1e6 for v in df['TT_del_axis'] ])
    tstop = pandas.Series([float(v) if v != 'NA' else 1e6 for v in df['TT_stop_axis'] ])
    for i in range(len(who)):
        to_kill = who[i]
        td = tdel[to_kill]
        tdd = pandas.Series([when[i]] * len(td))
        td.index = tdd.index
        ts = tstop[to_kill]
        tss = tdd - TT_stop_del
        ts.index = tss.index
        newt = [str(round(t)) if t != 1e6 else 'NA' for t in pandas.DataFrame([td, tdd]).min()]
        newstop = [str(round(t)) if t != 1e6 else 'NA' for t in pandas.DataFrame([ts, tss]).min()]
        df['TT_del_axis'][to_kill] = newt
        df['TT_stop_axis'][to_kill] = newstop
        #maj HS_final
        d = df[to_kill]
        newhs = []
        for a in range(len(d)):
            da = d.iloc[a]
            phen = dp[dp['id_phen'] == da['id_phen']]
            x = da['TT_col_phytomer1'] + phen['dTT_col_phytomer']
            y = da['N_phytomer'] * phen['index_rel_phytomer']
            newhs.append(numpy.interp(float(da['TT_stop_axis']), x, y))
        df['HS_final'][to_kill] = newhs
    devT['axeT'] = df.to_dict('list')
    
    return devT
    
def adjust_density(devT, density, TT_stop_del = 2.8 * 110, adjust_tiller = None):
    nplants = len(set(devT['axeT']['id_plt']))
    tdeath = time_of_death(nplants, density)
    dead = random.sample(set(devT['axeT']['id_plt']),len(tdeath))
    df = pandas.DataFrame(devT['axeT'])
    who = [df['id_plt'] == dead[i] for i in range(len(dead))]
    devT = kill_axis(devT, who, tdeath, TT_stop_del = TT_stop_del)
        
    return devT
    
def adjust_tiller_survival(devT, cohort_survival, TT_stop_del = 2.8 * 110):
    """
    Make (all) tillers of a living plant die along survival table
    """
    df = pandas.DataFrame(devT['axeT'])
    df['iaxe'] = range(len(df))
    for cohort in cohort_survival:
        naxes = len(df[df['id_cohort'] == cohort])
        tdeath = time_of_death(naxes, cohort_survival[cohort], relative_density=True)
        dead = random.sample(df['iaxe'][df['id_cohort'] == cohort], len(tdeath))
        who = [df['iaxe'] == dead[i] for i in range(len(dead))]
        devT = kill_axis(devT, who, tdeath, TT_stop_del = TT_stop_del)
        
    return devT
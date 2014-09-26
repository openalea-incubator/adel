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
from alinea.adel.AdelR import devCsv



def plantgen_to_devT(pgen):
    """ Creates devT tables from plantgen dict
    """
       
    axeT_, dimT_, phenT_, phenT_abs, dimT_abs, dynT_, phenT_first, HS_GL_SSI_T, tilleringT, cardinalityT, config = gen_adel_input_data(**pgen)
    
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    devT = devCsv(axeT, dimT, phenT)
    return devT, {'phenT_abs':phenT_abs, 'dimT_abs':dimT_abs, 'phenT_first':phenT_first, 'HS_GL_SSI_T':HS_GL_SSI_T, 'tilleringT':tilleringT, 'cardinalityT':cardinalityT, 'config':config}


def flat_list(nested_list):
    return list(chain.from_iterable(nested_list))
   
def modalities(nff):
    m1,m2 = int(nff), int(nff) + 1
    p = m1 + 1 - nff
    return {m1: p, m2: 1 - p}
    
def cardinalities(proba, n):
    card  = {k:int(v*n) for k,v in proba.iteritems()}
    sorted_p = sorted(proba.iteritems(), key=operator.itemgetter(1), reverse=True)
    missing = int(n - sum(card.values()))
    new = [sorted_p[i][0] for i in range(missing)]
    for k in new:
        card[k] += 1
    return {k:v for k,v in card.iteritems() if v > 0}
 
def _parent(axe):
    return '.'.join(axe.rsplit('.')[:-1])  
 
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
    
    
    cohort_axis = {k:cardinalities(axis_proba[k], sum(v.values())) for k, v in cohort_modalities.iteritems()}
    cohort_ax = {k:flat_list(map(lambda x: [x[0]] * int(x[1]),v.items())) for k, v in cohort_axis.iteritems()}
    
    axis = {k:zip(cohort_mods[k],cohort_ax[k]) for k in cohort_mods}
    axis_list = [map(lambda x: (k,x[0],x[1]), v) for k,v in axis.iteritems()]
    
    return flat_list(axis_list)
    
    
def plant_list(axis, nplants = 2):

# TO DO for more robust name selection : test if any parent of same cohort is present on the plant to filter candidates,and, in update find the matching prent and then choose a compatible name

    def _choose_plant(axe_name, plantlist):   
        candidates = filter(lambda x: (axe_name not in x) and (_parent(axe_name) in x or _parent(axe_name) == ''), plantlist)
        if len(candidates) == 0 and _parent(axe_name) == '':
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
    
def dynT_user(MS_parameters = {'a_cohort':1. / 110.,'TT_col_0':160.,'TT_col_N_phytomer_potential':1100,'n0':4.5,'n1':2.5,'n2':5}, primary_tillers = ['T%d'%(i) for i in range (1,4)]):
    p = primary_tillers
    p.sort()
    idaxis = ['MS'] + p
    df = pandas.DataFrame(index=idaxis,
                          columns=['id_axis','a_cohort','TT_col_0','TT_col_N_phytomer_potential','n0','n1','n2'],
                          dtype=float)
    df.ix['MS'] = pandas.Series(MS_parameters)
    df['id_axis'] = idaxis
    df['TT_col_N_phytomer_potential'] = df['TT_col_N_phytomer_potential'][0]
    df = df.reset_index(drop=True)
    return df
  
    
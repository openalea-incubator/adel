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
from scipy.interpolate import interp1d
import random

from alinea.adel.exception import *
from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data, plantgen2adel
from alinea.adel.plantgen import tools, params
from alinea.adel.AdelR import devCsv

#define classes for handling models

class TillerEmission(object):
    """ An tiller emission model based on plantgen random sampling strategy
    """
    
    def __init__(self, primary_tiller_probabilities = {'T1':0.95, 'T2':0.85, 'T3':0.75}, MS_leaves_number_probabilities={'11': 0.3, '12': 0.7}, inner_parameters ={}):
        self.child_cohort_delay = inner_parameters.get('FIRST_CHILD_DELAY', params.FIRST_CHILD_DELAY)
        self.primary_tiller_probabilities = {k:v for k,v in primary_tiller_probabilities.iteritems() if v > 0 }
        self.a1_a2 = inner_parameters.get('SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS',params.SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS)
        self.MS_leaves_number_probabilities = MS_leaves_number_probabilities
    

    def nff(self, cohort, nff_MS):
        if cohort == 1:
            nff = nff_MS
        else:
            nff = tools.calculate_tiller_final_leaves_number(nff_MS, cohort, self.a1_a2)
        nff_int = int(nff)
        if random.random() <= (nff - nff_int):
            nff += 1
        return int(nff)
        
        
    def new_plant(self,id_plant=1):
        """ Return an axis sample for one  plant, following probabilities given as input
        """
        #id_plant = range(start, start + nplants + 1)
        cohort_probabilities = tools.calculate_decide_child_cohort_probabilities(self.primary_tiller_probabilities)
        cohorts =  tools.decide_child_cohorts(cohort_probabilities, self.child_cohort_delay)
        id_cohort, id_axis = zip(*cohorts)
        nff_MS = tools.calculate_MS_final_leaves_number(self.MS_leaves_number_probabilities)
        nff = [self.nff(i,nff_MS) for i in id_cohort]
        return pandas.DataFrame({'id_plt': id_plant, 'id_cohort': id_cohort, 'id_axis' : id_axis, 'N_phytomer_potential': nff})

    def emission_parameters(self):
        return {'decide_child_axis_probabilities':self.primary_tiller_probabilities, 'MS_leaves_number_probabilities':self.MS_leaves_number_probabilities}

class TillerRegression(object):
    """ Tiller regression model
    """
        
    def __init__(self, ears_per_plant = 2.5, n_elongated_internode = 4, delta_stop_del = 2):
        self.ears_per_plant = ears_per_plant
        self.n_elongated_internode = n_elongated_internode
        self.delta_stop_del = delta_stop_del
        self.delta_reg = float(params.DELAIS_REG_MONT) / 110
        

    def hs_debreg(self, nff):
        return nff - self.n_elongated_internode + self.delta_reg
        
class HaunStage(object):
    """ Handle HaunStage = f (ThermalTime) fits
    """
    
    def __init__(self, a_cohort = 1. / 110., TT_hs_0 = 0):
        self.a_cohort = a_cohort
        self.TT_hs_0 = TT_hs_0
        
    def __call__(self, TT):# HS
        return (numpy.array(TT) - self.TT_hs_0) * self.a_cohort
        
    def TT(self, HS):
        return self.TT_hs_0 + numpy.array(HS) / self.a_cohort
        
    def TTem(self, TT):
        return numpy.array(TT) - self.TT_hs_0
        
    def phyllochron(self):
        return 1. / self.a_cohort

class GreenLeaves(object):
    """
    An object interface to a plantgen derived Green Leaf model
    This variant is for GL=f(HS) fits (instead of GL=f(TT)) with varying nff
    """

    def __init__(self, GL_start_senescence=4.8, GL_bolting=3.2, GL_HS_flag=5.8, n_elongated_internode= 4, curvature = -0.01):
        # rename parameters using plantgen terminology
        self.n0 = GL_start_senescence
        self.n1 = GL_bolting
        self.n_elongated_internode = n_elongated_internode
        self.n2 = GL_HS_flag
        self.a = curvature
      
    def hs_t1(self, nff=12):
        return nff - self.n_elongated_internode
    
    def linear_fit(self, nff=12):
        hs_t2 = nff
        return interp1d([0, self.n0, self.hs_t1(nff), hs_t2],[0, self.n0, self.n1, self.n2], bounds_error=False, fill_value=0)
    
    def polynomial_fit(self, nff=12):
        hs_t2 = nff
        c = (self.n2 - self.n1) / (hs_t2 - self.hs_t1(nff)) - 1
        pol = numpy.poly1d([self.a, 0.0, c, self.n2])
        def _fit(hs):
            gl = numpy.where(hs <= hs_t2, 0, pol(hs - hs_t2))
            return numpy.where(gl >=0, gl, 0)
        return _fit
        
    def curve(self, nff=12):
        lin = self.linear_fit(nff)
        pol = self.polynomial_fit(nff)
        def _curve(hs):
            return lin(hs) + pol(hs)            
        return _curve
        
    def HS_GL_sample(self, nff=12):
        curve = self.curve(nff)
        hs = numpy.linspace(nff, 2*nff,20)
        df = pandas.DataFrame({'HS':hs,'GL':curve(hs)})
        return df.loc[df['GL'] > 0,:]
            
      
class WheatDimensions(object):
    """ A dimension generator model based on simple scaling of a reference dataset
    """
    
    def __init__(self):
        self.ref = pandas.DataFrame({'index_phytomer': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
                                   'L_blade': [9.45, 9.27,8.04, 9.6, 11.26, 12.33, 14.07, 17.25, 17.75, 16.27, 10.74],
                                   'W_blade': [0.38, 0.37, 0.45, 0.58, 0.75, 0.99, 1.02, 1.06, 1.08, 1.1, 1.2],
                                   'L_sheath': [2.84, 2.93, 3.24, 3.89, 4.48, 6.81, 8.89, 9.58, 10.12, 11.07, 14.72],
                                   'W_sheath': [0.16, 0.18, 0.21, 0.23, 0.26, 0.28, 0.31, 0.34, 0.36, 0.39, 0.39],
                                   'L_internode': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.55, 4.45, 6.43, 11.3, 15.35], 
                                   'W_internode': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.06, 0.13, 0.19, 0.3]})
                        
        irel = [0]  + (self.ref['index_phytomer'] / 11.).tolist()
        self.predict = {k:interp1d(irel,[0] + self.ref[k].tolist()) for k in self.ref.columns}
    
   
    def dimT_user_table(self, nff=12, scale=1.):
        irel = numpy.arange(1, nff + 1) / float(nff)
        df = pandas.DataFrame({k:self.predict[k](irel) * scale for k in self.predict})
        df['index_phytomer'] = range(1, nff + 1)
        return df
        
        
      
class PlantGen(object):
    """ A class interface for generating plants one by one with plantgen
    """
    
    def __init__(self, HSfit=None, GLfit=None, Dimfit=None, Regression = None, Emission=None, inner_parameters={}):
    
        #defaults for fitted models
        if HSfit is None:
            HSfit = HaunStage(a_cohort = 1. / 110., TT_hs_0 = 0)           
        if GLfit is None:
            GLfit = GreenLeaves(GL_start_senescence=4.4, GL_bolting=1.5, GL_HS_flag=5, n_elongated_internode= 4, curvature = -0.01)        
        if Dimfit is None:
            Dimfit = WheatDimensions()            
        if Regression is None:
            Regression = TillerRegression(ears_per_plant=2.5, n_elongated_internode=4, delta_stop_del= 2)            
        if Emission is None:
            Emission = TillerEmission({'T1':0.95, 'T2':0.85, 'T3':0.75}, {'11':0.3, '12':0.7}, inner_parameters)
            
        self.inner_parameters = inner_parameters
        self.HSfit = HSfit
        self.GLfit = GLfit
        self.Dimfit = Dimfit
        self.Emission = Emission
        self.Regression = Regression
        
        #setup base configuration for plantgen interface for one plant
        self.base_config = {}
        self.base_config['inner_params'] = inner_parameters
        self.base_config['inner_params']['EMF_1_MS_STANDARD_DEVIATION'] = 0 #strictly respect TT_hs_0. 
        self.base_config.update({'plants_number': 1,
                            'plants_density': 1.})
                  
        self.count=0
                       
    def reset(self):
        self.count = 0
     
    def config(self, axeT = None):
        """ Generate plantgengen configuration dict for a plant decribed in axeT
            if axeT is None, a random plant with id is generated using the plant_generator
        """
           
        config = {}
        config.update(self.base_config) #avoid side effects on base_config
        
        if axeT is None:
            self.count += 1
            axeT = self.Emission.new_plant(self.count)
          
        # consign emision parameters for information/plantgen checkings only (will not be used as axeT will be send as input)
        emission_parameters = self.Emission.emission_parameters()
        config.update({'decide_child_axis_probabilities': emission_parameters['decide_child_axis_probabilities'],
                        'MS_leaves_number_probabilities': emission_parameters['MS_leaves_number_probabilities']})
                                
        #nff dependent fits
        nff = axeT['N_phytomer_potential'].max()
        
        # Tillering
        hs_deb_reg = self.Regression.hs_debreg(nff)
        axeT_user = self.axeT_user_table(axeT)
        config.update({'ears_density' : self.Regression.ears_per_plant, 
                        'delais_TT_stop_del_axis': self.Regression.delta_stop_del * self.HSfit.phyllochron(),
                        'TT_regression_start_user': self.HSfit.TT(hs_deb_reg),
                        'axeT_user':axeT_user
                       })
                       
        # Dimensions
        config['dimT_user'] = self.Dimfit.dimT_user_table(nff)
        
        # Dynamic              
        config['dynT_user'] = self.dynT_user_table(nff)
        
        GL = self.GLfit.HS_GL_sample(nff)
        GL['TT'] = self.HSfit.TT(GL['HS'])
        GL = dict(zip(GL['TT'],GL['GL']))
        config['GL_number'] = GL
        
        hs_t1 = self.GLfit.hs_t1(nff)
        config['TT_t1_user'] = self.HSfit.TT(hs_t1)
        
        return config
        
            
    def axeT_user_table(self, axeT):
    
        df = pandas.DataFrame(index=range(len(axeT['id_plt'])),
                                               columns=['id_plt', 'id_cohort', 'id_axis', 'N_phytomer_potential', 'N_phytomer', 'HS_final', 'TT_stop_axis', 'TT_del_axis', 'id_dim', 'id_phen', 'id_ear', 'TT_app_phytomer1', 'TT_col_phytomer1', 'TT_sen_phytomer1', 'TT_del_phytomer1'],
                                               dtype=float)

        df['id_plt'] = axeT['id_plt']
        df['id_axis'] = axeT['id_axis']
        df['id_cohort'] = axeT['id_cohort']
        df['N_phytomer_potential'] = axeT['N_phytomer_potential']
        df['id_phen'] = df['id_cohort'] * 100 + df['N_phytomer_potential']
        
        df= df.sort(['id_plt','id_cohort','id_axis'])
        return df
      
    def dynT_user_table(self, nff):
    
        MS_parameters = {'a_cohort': self.HSfit.a_cohort,
                         'TT_hs_0': self.HSfit.TT_hs_0,
                         'TT_hs_N_phytomer_potential': self.HSfit.TT(nff),
                         'n0': self.GLfit.n0,
                         'n1': self.GLfit.n1,
                         'n2': self.GLfit.n2}
        tillers_probabilities = self.Emission.primary_tiller_probabilities
        cohort_probabilities = tools.calculate_decide_child_cohort_probabilities(tillers_probabilities)
        primary_tillers = tillers_probabilities.keys()
        primary_tillers.sort()
        idaxis = ['MS'] + primary_tillers
        df = pandas.DataFrame(index=idaxis,
                              columns=['id_axis','a_cohort','TT_hs_0','TT_hs_N_phytomer_potential','n0','n1','n2'],
                              dtype=float)
        df.ix['MS'] = pandas.Series(MS_parameters)
        df['id_axis'] = idaxis
        cohorts = cohort_probabilities.keys()
        cohorts.sort()
        hs_flag = df['TT_hs_N_phytomer_potential'][0] + (numpy.array(cohorts) - 1) * 1. / (4 * self.HSfit.a_cohort) #cf pgen core line 975
        df['TT_hs_N_phytomer_potential'] = [df['TT_hs_N_phytomer_potential'][0]] + hs_flag.tolist()
        df = df.reset_index(drop=True)
        return df
    
    def dimT_user_table(self, nff=11, scale=1.):
        ref = pandas.DataFrame({'index_phytomer': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
                                   'L_blade': [9.45, 9.27,8.04, 9.6, 11.26, 12.33, 14.07, 17.25, 17.75, 16.27, 10.74],
                                   'W_blade': [0.38, 0.37, 0.45, 0.58, 0.75, 0.99, 1.02, 1.06, 1.08, 1.1, 1.2],
                                   'L_sheath': [2.84, 2.93, 3.24, 3.89, 4.48, 6.81, 8.89, 9.58, 10.12, 11.07, 14.72],
                                   'W_sheath': [0.16, 0.18, 0.21, 0.23, 0.26, 0.28, 0.31, 0.34, 0.36, 0.39, 0.39],
                                   'L_internode': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.55, 4.45, 6.43, 11.3, 15.35], 
                                   'W_internode': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.06, 0.13, 0.19, 0.3]})
        
        #to do scale dimT
        return ref
    
    def tables(self, axe_T=None):
        config = self.config(axe_T)
        axeT_, dimT_, phenT_, phenT_abs, dimT_abs, dynT_, phenT_first, HS_GL_SSI_T, tilleringT, cardinalityT, config = gen_adel_input_data(**config)
        axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
       
        return {'adelT': (axeT, dimT, phenT), 'phenT_abs':phenT_abs, 'dimT_abs':dimT_abs, 'phenT_first':phenT_first, 'HS_GL_SSI_T':HS_GL_SSI_T, 'tilleringT':tilleringT, 'cardinalityT':cardinalityT, 'config':config}
 

    

# extensions
 
def adelT_to_devT(pgen):
    """ Creates devT tables from plantgen dict
    """
    devT = devCsv(*pgen['adelT'])
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
 

 
def axis_list(TilleringModel,  nplants = 2):
    """ compute cardinalities of axis in a stand of n plants
    The strategy used here is based on deterministic rounding, and differs from the one used in plantgen (basd on random sampling). Difference are expected for small plants numbers
    """
    df = TilleringModel.emited_cohorts()
    df = emited_cohorts
    df = df.set_index('cohort') 
    cohort_cardinalities = {c:round(df.ix[c,'total_axis'] * nplants) for c in df.index}
    
    modal_proba = {c:modalities(df.ix[c,'nff']) for c in df.index}
    
    cohort_modalities = {k:cardinalities(modal_proba[k],v) for k,v in cohort_cardinalities.iteritems()}
    cohort_mods = {k:flat_list(map(lambda x: [x[0]] * int(x[1]),v.items())) for k, v in cohort_modalities.iteritems()}

    
    p = TilleringModel.theoretical_probabilities()
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
 
 
def axeT_user(nplants, TilleringModel):
    
    axis = axis_list(TilleringModel, nplants)
    plants = plant_list(axis, nplants) 
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
    update devT tables by killing (deleting) axis at pre-defined time, hence stoping them growing at when - TT_stop_del
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
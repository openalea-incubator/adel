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

# some new tools

def flat_list(nested_list):
    return list(chain.from_iterable(nested_list))

def random_round(decimal):  
    """ ramdom rounding of a decimal to int so that mean(rounded) = decimal for a large series
    """
    rounded = int(decimal)
    if random.random() <= (decimal - rounded):
        rounded += 1
    return rounded    
    
def _parent(axe):
    return '.'.join(axe.rsplit('.')[:-1])  
    
def modalities(nff):
    m1,m2 = int(nff), int(nff) + 1
    p = m1 + 1 - nff
    return {m1: p, m2: 1 - p}
    
def cardinalities(proba, n, parents = None):
    """ Return optimised cardinalities per modality from a dict {modality:probability} and a global carinality n
        parents allows for control of botanical constraint in the specific case of tiller sampling
    """
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
   
def card2list(card,shuffle=True):
    liste = flat_list([[int(key)] * val for key,val in card.iteritems()])
    if shuffle:
        random.shuffle(liste)
    return liste

def plant_list(axis, nplants = 2):
    """ contruct plants by picckling/grouping axis
    """
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
    
#define classes for structuring/handling the different botanical models found in pgen
                
class TillerEmission(object):
    """ A class interface to plantgen tiller emmission model
    """
    
    def __init__(self, primary_tiller_probabilities = {'T1':0.95, 'T2':0.85, 'T3':0.75, 'T4':0.4}, inner_parameters ={}):
        
        self.primary_tiller_probabilities = {k:v for k,v in primary_tiller_probabilities.iteritems() if v > 0 }
        self.cohort_probabilities = tools.calculate_decide_child_cohort_probabilities(self.primary_tiller_probabilities)
              
        self.child_cohort_delay = inner_parameters.get('FIRST_CHILD_DELAY', params.FIRST_CHILD_DELAY)
        self.a1_a2 = inner_parameters.get('SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS',params.SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS)
        self.tiller_delays = inner_parameters.get('LEAF_NUMBER_DELAY_MS_COHORT',params.LEAF_NUMBER_DELAY_MS_COHORT)    

        
    def cohort_delays(self):
        """ compute cohort emergence delays
        emergence delay = delays HS_0 MS -> HS_0 tillers HS=0 <=> emergence leaf
        """
        cohort_delays = {(int(k.lstrip('T')) + 3):v for k,v in self.tiller_delays.iteritems()}
        cohort_delays[1] = 0
        return cohort_delays
            
    def theoretical_probabilities(self):
        """ theoretical probability of emergence of tillers/cohorts"""

        _,res = tools.calculate_theoretical_cardinalities(1, 
                                                      self.cohort_probabilities,
                                                      self.primary_tiller_probabilities,
                                                      self.child_cohort_delay)
        return res
        
    def final_leaf_numbers(self, MS_nff=12.):
        """ returns the decmal final leaf numbers of cohorts from the decimal number of leaves on main stem
        """
        ms_nff = [(1,MS_nff)]
        cohort_nff = [(cohort,tools.calculate_tiller_final_leaves_number(MS_nff, cohort, self.a1_a2)) for cohort in self.cohort_probabilities]
        return dict(ms_nff + cohort_nff)  

    def emission_table(self, MS_nff = 12., hs_debreg=None, max_order=None):
        """ Computes cohort number, botanical position and theoretical probabilities of emergence all tillers using botanical position and probabilities of emergence of primary ones
        """
        proba = self.theoretical_probabilities()
        emergence_delays = self.cohort_delays()
        nff = self.final_leaf_numbers(MS_nff)
        data = zip(*[(k[0],k[1],v) for k,v in proba.iteritems()])
        df = pandas.DataFrame({'axis': data[1],
                               'cohort': data[0],
                               'nff': [nff[c] for c in data[0]],
                               'delay':[emergence_delays[c] for c in data[0]],
                               'probability':data[2]})
        df = df.sort(['cohort','axis'])
        
        if max_order != None:
            def _order(axe):
                return len(axe.rsplit('.'))
            df = df[map(lambda x: _order(x) <= max_order,df['axis'])]
           
        if hs_debreg != None:
            df = df[df['delay'] <= hs_debreg]
                   
        return df
    
    def cohort_emission_table(self, emission_table):
        how = {'cohort':numpy.mean,'delay':numpy.mean, 'nff':numpy.mean, 'probability':numpy.sum}
        df = emission_table.groupby('cohort', group_keys=False).agg(how)
        return df
        
    #deprecated
   
    def emission_parameters(self):
        return {'decide_child_axis_probabilities':self.primary_tiller_probabilities, 'MS_leaves_number_probabilities':self.MS_leaves_number_probabilities}

class TillerRegression(object):
    """ Tiller regression model
    """
        
    def __init__(self, ears_per_plant = 2.5, n_elongated_internode = 4, delta_stop_del = 2, inner_parameters ={}):
        self.ears_per_plant = ears_per_plant
        self.n_elongated_internode = n_elongated_internode
        self.delta_stop_del = delta_stop_del
        # Express time-delta between tiller regression and bolting in phyllochronic units
        self.delta_reg = float(inner_parameters.get('DELAIS_REG_MONT', params.DELAIS_REG_MONT)) / 110

    def decimal_elongated_internode_number(self, internode_ranks, internode_lengths):
        """ estimate the (decimal) number of elongated internode as determined in Plantgen from internode dimension of most frequent axis
        """
        df = pandas.DataFrame({'rank':internode_ranks, 'length':internode_lengths})
        # keep only the non-zero lengths
        df = df[df['length'] > 0]
        # Fit a polynomial of degree 2 to, and get the coefficient of degree 0. 
        n = df['rank'].max() - numpy.polyfit(df['length'].values, df['rank'].values, 2)[2]
        return n   

    def hs_debreg(self, nff):
        return nff - self.n_elongated_internode + self.delta_reg
        
    def regression_table(self, cohort_emission_table):
        """ compute regression table for the differents cohorts
        cohort_emission_table is a emission table as generated by en emission model when grouped by cohort
        """
        
        # regression period for the disparition of axes
        nff = cohort_emission_table['nff'].max()
        hs_deb = self.hs_debreg(nff) + self.delta_stop_del
        hs_end = nff + self.delta_stop_del
        
        #regression rate
        cohorts = cohort_emission_table
        d_deb = cohorts['probability'].sum()
        d_end = self.ears_per_plant
        regression_rate = (d_end - d_deb) / (hs_end - hs_deb)
        
        # regression parameters per cohort
        # time of disparition of cohorts in case of complete regression
        # hypothesis : disparition of cohorts is sequential from last appeared to first
        reverse_cohorts = cohorts.sort_index(by=['delay'], ascending = False)
        loss_acc = reverse_cohorts['probability'].cumsum().values
        t_disp = hs_deb + loss_acc / abs(regression_rate)
        #actual fraction lost at t_disp
        f_disp = t_disp.copy()
        f_disp[t_disp < hs_end] = 1
        f_disp[t_disp >= hs_end] = 0
        #compute last fraction
        card = reverse_cohorts['probability'].values
        last = min(numpy.where(t_disp >= hs_end)[0])
        last_loss = (d_deb - d_end) - sum(card * f_disp) #total_lost - loss by completly disapeared cohorts
        f_disp[last] = last_loss / card[last]
        #
        t_disp = numpy.minimum(t_disp, hs_end)
        reg_table = pandas.DataFrame({'cohort':reverse_cohorts['cohort'].values,
                                      't_disp': t_disp,
                                      'f_disp': f_disp,
                                      't_start':[hs_deb] + t_disp.tolist()[:-1]})
        
        return reg_table[reg_table['f_disp'] > 0]    
        
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
        
class HS_flag(object):
    """ models for predicting HS flag of axes from most frequent Main stem
    """
    def __init__(self, most_frequent_MS_nff = 12, most_frequent_MS_a_cohort = 1. / 110, inner_params={}):
        self.most_frequent_nff = most_frequent_MS_nff
        self.a_cohort = most_frequent_MS_a_cohort
        #
        self.a1_a2 = inner_parameters.get('SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS',params.SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS)
    #  
    def dTT_MS_cohort(self, cohort):
        """ 
        'model' for thermal time delay between Most frequent MS flag leaf appearance and and its relative tillers
        (plantgen doest provide a model, but consider this as inpus)     
        id axis is tiller axis name (T0, T1...)
        """
        
        if cohort == 1:
            return 0
        else:
            return 60 + (cohort - 3) #valeurs maxwell mariem
        
    def dTT_nff(self, cohort, nff):
        """ delay between flag leaf appearance of the most frequent axis of a cohort and an axis of the same cohort with nff leaves 
        """
        if cohort == 1:#MS
            most_frequent_nff = self.most_frequent_nff
        else:
            most_frequent_nff = tools.calculate_tiller_final_leaves_number(self.most_frequent_nff, cohort, self.a1_a2)
            
        return 1. * (nff - most_frequent_nff) / (4 * self.a_cohort)
        
        
    def dTT_MS(self, axis_id, nff):
        """ delay between most frequent MS flag leaf ligulation and axis axis_id with nff leaves
        """
        if axis_id == "MS":
            cohort = 1
        else:
            cohort = int(id_axis.lstrip('T')) + 3
        return self.dTT_MS_cohort(cohort) + self.dTT_nff(cohort, nff)

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
    

class AxePop(object):
    """ An extended axe population generator based on plantgen axe generation methods
    """
        
    def __init__(self, MS_leaves_number_probabilities={'11': 0.3, '12': 0.7}, Emission=None, Regression=None, Damages=None):

        if Regression is None:
            Regression = TillerRegression(ears_per_plant=2.5, n_elongated_internode=4, delta_stop_del= 2)            
        if Emission is None:
            Emission = TillerEmission({'T1':0.95, 'T2':0.85, 'T3':0.75, 'T4':0.4})
            
        self.Emission = Emission
        self.Regression = Regression
        self.Damages = Damages
        self.MS_probabilities = MS_leaves_number_probabilities
      
    
    def mean_nff(self):
        return sum([int(k) * v for k,v in self.MS_probabilities.iteritems()])
        
    def mode_nff(self):
        v = self.MS_probabilities.values()
        nff = self.MS_probabilities.keys()
        return int(nff[v.index(max(v))])        
      
    def random_population(self, nplants=1):
        """ Ramdomly generate a population of axes (reproduce plantgen original behavior)
        """
        plants = []
        cohort_decimal_nff={int(k):self.Emission.final_leaf_numbers(int(k)) for k in self.MS_probabilities}
        for i in range(nplants): 
            id_plant = i + 1
            # Random sampling of axes
            axes =  tools.decide_child_cohorts(self.Emission.cohort_probabilities, self.Emission.child_cohort_delay)
            id_cohort, id_axis = zip(*axes)
            # random sampling of nff for main stem
            nff_MS = tools.calculate_MS_final_leaves_number(self.MS_probabilities)
            # nfff on tillers : random rounding of decimal nff of cohort for that nff
            nff = [random_round(cohort_decimal_nff[nff_MS][c]) for c in id_cohort]
            plants.append(pandas.DataFrame({'id_plt': id_plant, 'id_cohort': id_cohort, 'id_axis' : id_axis, 'N_phytomer_potential': nff}))
        return pandas.concat(plants) 

    def control_population(self, nplants=1, max_order=None):
        """ Generate an axe population using deteministic rounding
        """
        nff = self.mean_nff()
        #regression start is set/defined at population level
        hs_debreg = self.Regression.hs_debreg(nff)
        table = self.Emission.emission_table(nff, hs_debreg=hs_debreg, max_order=max_oder)               
        cohort_table = self.Emission.cohort_emission_table(table)
        
        cohort_cardinalities = numpy.round(cohort_table['probability'] * nplants)
        
        groups = table.groupby('cohort')
        axis_proba = groups.apply(lambda x: dict(zip(x['axis'].values, x['probability'].values / sum(x['probability'])))).to_dict()
        parents={'':0}
        axis_cardinalities = {k:cardinalities(axis_proba[k], v,parents) for k, v in cohort_cardinalities.iteritems()}
    
        axis_list = flat_list([flat_list(map(lambda x: [(k,x[0])] * int(x[1]),v.items())) for k, v in axis_cardinalities.iteritems()])
        plist = plant_list(axis_list, nplants)
        plant_axes = map(lambda x: [(v[0],k) for k,v in x.iteritems()],plist)
        
        nff_MS_cardinalities = cardinalities(self.MS_probabilities, nplants)
        nff_MS = card2list(nff_MS_cardinalities)
        cohort_decimal_nff = {int(k):self.Emission.final_leaf_numbers(int(k)) for k in self.MS_probabilities}
        cohort_nff_modalities = {k:{kk:modalities(vv) for kk,vv in v.iteritems()} for k,v in cohort_decimal_nff.iteritems()}
        cohort_nff_cardinalities = {int(k):{kk:cardinalities(vv, int(float(nff_MS_cardinalities[k]) / nplants * cohort_cardinalities[kk]) + 1) for kk,vv in cohort_nff_modalities[int(k)].iteritems()} for k in nff_MS_cardinalities}
        cohort_nff = {k:{kk:card2list(vv) for kk,vv in v.iteritems()} for k,v in cohort_nff_cardinalities.iteritems()}
        
        plants = []
        for i in range(nplants): 
            id_plant = i + 1
            axes =  plant_axes[i]
            id_cohort, id_axis = zip(*axes)
            # deterministic sampling of nff for main stem 
            nff_p = nff_MS.pop()
            # nfff on tillers
            nff = [cohort_nff[nff_p][c].pop() for c in id_cohort]
            plants.append(pandas.DataFrame({'id_plt': id_plant, 'id_cohort': id_cohort, 'id_axis' : id_axis, 'N_phytomer_potential': nff}))
        df=pandas.concat(plants)
        df= df.sort(['id_plt','id_cohort','id_axis'])
        return df

        
class PlantGen(object):
    """ A class interface for generating plants one by one with plantgen
    """
    
    def __init__(self, HSfit=None, GLfit=None, Dimfit=None, Regression = None, Emission=None, inner_parameters={}):
    
        #defaults for fitted models
        if HSfit is None:#HS fit is for the most frequent main stem
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
        self.HS_flag = HS_flag(most_frequent_MS_nff = self.Emission.most_frequent_nff_MS(), most_frequent_MS_a_cohort = self.HSfit.a_cohort, inner_params=self.inner_parameters)
        
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
          
        # emision parameters for plantgen (used for hs_flag model among others)
        emission_parameters = self.Emission.emission_parameters()
        config.update({'decide_child_axis_probabilities': emission_parameters['decide_child_axis_probabilities'],
                        'MS_leaves_number_probabilities': emission_parameters['MS_leaves_number_probabilities']})
                                
        #nff dependent fits: TO DO : compute HSfit for current plant nff (a_cohort should be re-computed to match hs_flag predicted by hs_flag model)
        # alternatively, HSfit may contain parameter forfitted plant main stem, but then ms_probalities_nff  for pgen and should be one and for HSflag model too
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
        hs_flag = df['TT_hs_N_phytomer_potential'][0] #ICI : ajouter prediction HSflag
        df['TT_hs_N_phytomer_potential'] = [df['TT_hs_N_phytomer_potential'][0]] + hs_flag.tolist()
        df = df.reset_index(drop=True)
        return df
    
    
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



  
   


 
def axis_list(TilleringModel,  nplants = 2):
    """ compute cardinalities of axis in a stand of n plants
    The strategy used here is based on deterministic rounding, and differs from the one used in plantgen (basd on random sampling). Difference are expected for small plants numbers
    """
    
    df = TilleringModel.emited_cohorts()
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
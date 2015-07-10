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
import alinea.adel.data_samples as adel_data

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
    
def get_normal_dist(nb_plants=10, sigma=30.):
    """ Calculate the "best possible" sample for a given number of plants 
    as a function of the standard deviation sigma of a normal centered ditribution. """
    N = 1000
    norm_list = numpy.random.normal(scale=sigma, size=N)
    h = numpy.histogram(norm_list, bins=nb_plants)
    classes = h[1]
    distri_plants = h[0]*nb_plants/float(N)
    round_distri = [numpy.round(d) for d in distri_plants]
    missing_pl = nb_plants - sum(round_distri)
    while missing_pl > 0:
        gap = round_distri - distri_plants
        round_distri[numpy.argmin(gap)] += 1
        missing_pl = nb_plants - sum(round_distri)
    while missing_pl < 0:
        gap = distri_plants - round_distri
        round_distri[numpy.argmin(gap)] -= 1
        missing_pl = nb_plants - sum(round_distri)
    return numpy.hstack([numpy.linspace(h[1][i], h[1][i+1], d+2)[1:-1] for i, d in enumerate(round_distri) if d>0])

def _order(axe):
    if axe is 'MS':
        return 0
    else:
        return len(axe.rsplit('.'))    
    
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
        proba = {k:v for k,v in proba.iteritems() if (card[k] <= parents.get(_parent(k),n) or _parent(k) == '')}
    
    missing = int(n - sum(card.values()))
    while (missing > 0):
        # current frequencies
        freq = {k:float(v) / n for k,v in card.iteritems() if k in proba}
        # diff with probabilities
        dp = {k:abs(freq[k] - proba[k]) for k in freq}
        sorted_p = sorted(dp.iteritems(), key=operator.itemgetter(1), reverse=True)
        k = sorted_p[0][0]
        card[k] += 1
        if parents is not None:
            if (card[k] >= parents.get(_parent(k),n) and _parent(k) != ''):
                proba.pop(k)
        missing -= 1   
            
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
    
def t_death(n, t_start, t_end):
    """ regular sampling of time of death between t_start and t_end
    """
    step = 1. / n
    at_n = numpy.linspace(step / 2., n - step / 2., n)
    return numpy.interp(at_n, [0, n], [t_start, t_end]).round(2).tolist()
    
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
            df = df[map(lambda x: _order(x) <= max_order,df['axis'])]
           
        if hs_debreg != None:
            df = df[df['delay'] <= hs_debreg]
                   
        return df
    
    def cohort_emission_table(self, emission_table):
        how = {'cohort':numpy.mean,'delay':numpy.mean, 'nff':numpy.mean, 'probability':numpy.sum}
        df = emission_table.groupby('cohort', group_keys=False).agg(how)
        return df

    def curve_emission_table(self, emission_table):
        emission_table['primary'] = map(lambda x: not '.' in x, emission_table['axis'])
        def _fun(x):
            d = {'cohort': x['cohort'].mean(), 
                 'delay': x['delay'].mean(),
                 'primary_axis':x[x['primary']]['probability'].sum(),
                 'other_axis' : x[~x['primary']]['probability'].sum(),
                 'nff' : x['nff'].mean()
                 }
            out = pandas.DataFrame(d,index = [x.index[0]])
            out['total_axis'] = out['primary_axis'] + out['other_axis']
            return out
        return emission_table.groupby(['cohort'],group_keys=False).apply(_fun)
        
    def emission_curves(self, emission_table, include_MS=True, delta=0.1):
        """ return interpolation functions for axe emission (primary, other, total) as a function of haun stage
        """
        cohorts = self.curve_emission_table(emission_table)        
        if not include_MS:
            cohorts = cohorts.loc[cohorts['cohort'] > 1,:] 
            
        hs = reduce(lambda x,y:x+y,[[hs - delta / 2,hs + delta / 2] for hs in cohorts['delay']])        
        curves = {}
        for w in ('primary','other', 'total'):
            card = cohorts[w + '_axis'].values
            cum = card.cumsum()
            em = reduce(lambda x,y:x+y,[[cum[i] - card[i],cum[i]] for i in range(len(card))])
            curves[w] = interp1d([-delta] + hs + [20],[0] + em + [em[-1]])
        return curves

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

    def regression_curves(self, cohort_emission_table, curve_emission_table, damage_table=None):
        """ interpolation function for natural axe regression (number of axes lost) as a function of haun stage
        """
        reg_table = self.regression_table(cohort_emission_table)
        cohorts = curve_emission_table
        regressing_cohorts = reg_table.merge(cohorts)
        #
        #  Reduce regression to compensate for damage loss (damages shoud not influence fitted final axe density nor final ear number)
        damages = damage_table
        if damages is not None:
            regressing_cohorts = regressing_cohorts.set_index('cohort')
            damages = damages.set_index('cohort')
            cohorts = cohorts.set_index('cohort')
            # for damaged regressing cohort, reduce f_disp to compensate for part/all of f_comp
            for c in damages.index:
                if c in regressing_cohorts.index:
                    comp = min(damages['f_damaged'][c], regressing_cohorts['f_disp'][c])
                    regressing_cohorts['f_disp'][c] = regressing_cohorts['f_disp'][c] - comp
                    damages['f_damaged'][c] = damages['f_damaged'][c] - comp
            # remaining compensation obtained by reducing f_disp of regressing cohorts, older first
            if sum(damages['f_damaged']) > 1e-6:
                damages = damages[damages['f_damaged'] > 0]
                regressing_cohorts = regressing_cohorts.sort_index(by=['delay'], ascending = True) 
                for c in damages.index:
                    f_d = damages['f_damaged'][c]
                    for c_r in regressing_cohorts.index:
                        f_r = f_d * 1. * cohorts['total_axis'][c] / cohorts['total_axis'][c_r]
                        comp = min(f_r, regressing_cohorts['f_disp'][c_r])
                        regressing_cohorts['f_disp'][c_r] = regressing_cohorts['f_disp'][c_r] - comp
                        comp_d = comp * 1. * cohorts['total_axis'][c_r] / cohorts['total_axis'][c]
                        f_d = f_d - comp_d
                        if f_d < 1e-6:
                            break
                    assert f_d < 1e-6, 'Damages are too important to be compensated by reggressing tillers !'
                        
        #curve
        regressing_cohorts = regressing_cohorts.sort_index(by=['delay'], ascending = False)           
        hs = [0, regressing_cohorts['t_start'].tolist()[0]] +  regressing_cohorts['t_disp'].tolist() + [20]
        curves = {}            
        for w in ('primary','other', 'total'):
            card = regressing_cohorts[w + '_axis'] * regressing_cohorts['f_disp']
            n_d = card.cumsum()
            loss = [0] * 2 + card.cumsum().tolist() + [card.sum()]
            curves[w] = interp1d(hs, loss)
        return curves

        
               
class HaunStage(object):
    """ Handle HaunStage = f (ThermalTime) fits for mean plant and its nff modalities
    """
    # TO DO : add method for thermal time date of leaf emergence/ligulation
    
    def __init__(self, a_cohort = 1. / 110., TT_hs_0 = 0, std_TT_hs_0 = 0, nff = 12, dTT_nff = 1. / 110. / 4., dTT_cohort={'first': 30, 'increment': 10}):
        self.a_cohort = a_cohort
        self.TT_hs_0 = TT_hs_0
        self.nff = nff# nff of the mean plant
        self.dTT_nff = dTT_nff
        self.dTT_cohort = dTT_cohort
        self.std_TT_hs_0 = std_TT_hs_0
        
    def __call__(self, TT, nff=None):
        return self.HS(TT, nff)
        
    def HS(self, TT, nff=None):
        return (numpy.array(TT) - self.TT_hs_0) * self.a_nff(nff)
        
    def TT(self, HS, nff=None):
        return self.TT_hs_0 + numpy.array(HS) / self.a_nff(nff)
        
    def TTem(self, TT):
        return numpy.array(TT) - self.TT_hs_0
        
    def TTflag(self, nff = None, cohort = 1):
        if nff is None: 
            nff = self.nff
        return self.TT_hs_0 + self.nff / self.a_cohort + self.dTT_nff * (nff - self.nff) + self.dTT_MS_cohort(cohort)
    
    def HSflag(self, nff = None):
        if nff is None:
            return self.nff
        else:
            return nff
    
    def a_nff(self, nff=None):
        if nff is None:
            return self.a_cohort
        else:
            return nff / (self.TTflag(nff) - self.TT_hs_0)
    
    def phyllochron(self, nff=None):
        return 1. / self.a_nff(nff)
        
    def dHS_nff(self):
        return self.dTT_nff / self.phyllochron()
        
    def dTT_MS_cohort(self, cohort=1):
        """ delay between main stem mean flag leaf emergenece and mean flag leaf emergence of a cohort
        """
        if cohort == 1:
            return 0
        else:
            return self.dTT_cohort['first'] + self.dTT_cohort['increment'] * (cohort - 3)
    
    def curve(self, nff=None):
        if nff is None: 
            ymax = self.nff
        else:
            ymax = nff
        return interp1d([-1000., self.TT_hs_0, self.TTflag(nff), 3000.],[0., 0., ymax, ymax]) 

class GreenLeaves(object):
    """
    An object interface to a plantgen derived Green Leaf model
    This variant is for GL=f(HS_since_flag_leaf) fits (instead of GL=f(TT)) that allow to predict a curve for different nff
    """

    def __init__(self, GL_start_senescence=4.8, GL_bolting=3.2, GL_flag=5.8, n_elongated_internode= 4, curvature = -0.01):
        """ n_elongated_internode is the HS interval between bolting and flag leaf ligulation
        """
        # rename parameters using plantgen terminology
        self.n0 = GL_start_senescence
        self.n1 = GL_bolting
        self.n_elongated_internode = n_elongated_internode
        self.n2 = GL_flag
        self.a = curvature
    
    def fit_a(self, HS_since_flag, GL):
        """ Fit curvature coefficient from a HS, GL dataset
        """
        GLpol = pandas.DataFrame({'HS':HS_since_flag, 'GL':GL})
        GLpol = GLpol.ix[GLpol['HS'] > 0,:]
        c = (self.n2 - self.n1) / (self.n_elongated_internode) - 1
        fixed_coefs = [0.0, c, self.n2]
        a, rmse = tools.fit_poly(GLpol['HS'], GLpol['GL'], fixed_coefs, a_starting_estimate=self.a)
        self.a = a
        return a,rmse
      
      
    def hs_t1(self, hs_flag=12):
        return hs_flag - self.n_elongated_internode
    
    def linear_fit(self, hs_flag=12):
        hs_t2 = hs_flag
        return interp1d([0, self.n0, self.hs_t1(hs_flag), hs_t2],[0, self.n0, self.n1, self.n2], bounds_error=False, fill_value=0)
    
    def polynomial_fit(self, hs_flag=12):
        hs_t2 = hs_flag
        c = (self.n2 - self.n1) / (hs_t2 - self.hs_t1(hs_flag)) - 1
        pol = numpy.poly1d([self.a, 0.0, c, self.n2])
        def _fit(hs):
            gl = numpy.where(hs <= hs_t2, 0, pol(hs - hs_t2))
            return numpy.where(gl >=0, gl, 0)
        return _fit
        
    def curve(self, hs_flag=12):
        lin = self.linear_fit(hs_flag)
        pol = self.polynomial_fit(hs_flag)
        def _curve(hs):
            return lin(hs) + pol(hs)            
        return _curve
        
    def HS_GL_sample(self, hs_flag=12):
        curve = self.curve(hs_flag)
        hs = numpy.linspace(hs_flag, 2*hs_flag,20)
        df = pandas.DataFrame({'HS':hs,'GL':curve(hs)})
        return df.loc[df['GL'] > 0,:]
            
      
class WheatDimensions(object):
    """ A dimension generator model based on scaling of a reference dataset
    """
    
    def __init__(self, hsfit=None, dHS_max=2.4, nff=12, dHS_flag=0.25, phyllochron=110, scale=1):
        """
        dHS_flag, phyllochron and nff are read from hsfit when given, or from  explicit arguments if hsfit is None
        dHS_max : delta Haunstage between ligulation of longest leaf and flag leaf ligulation.
        scale allows an initial glogal re-scale of the dimensions. If scaleis a dict, then only specified dimensions are scaled
        """

        #defaults
        if hsfit is None:
            hsfit = HaunStage(a_cohort=1. / phyllochron, nff=nff, dTT_nff=dHS_flag * 1. / phyllochron)
        
        if not isinstance(scale,dict):
            scale={}
        
        self.hsfit = hsfit
        self.TTmax = hsfit.TT(hsfit.nff - dHS_max)
        self.ref = adel_data.wheat_dimension_profiles()
        self.xnref = self.ref['xn'].tolist()
        self.ref = self.ref.drop('xn',axis=1)
        self.predict = {k:interp1d(self.xnref, numpy.array(self.ref[k].tolist()) * scale.get(k,1)) for k in self.ref.columns}

    def fit_dimensions(self, data, fit_TTmax = True, fit_nff = True):
        """ 
        if fit_max is True, use blade length profile to re-adjust dHS_max
        if fit_nff is True, use dimension per nff to fit nff scale factors scale = 1 + scale_nff * delta_nff
        """
        
        if (fit_TTmax and 'L_blade' in data.columns):
            dat = data.dropna(subset=['L_blade'])
            xn = self.xn(dat['rank'])
            fit = numpy.polyfit(xn,dat['L_blade'],6)
            x = numpy.linspace(min(xn), max(xn), 500)
            y = numpy.polyval(fit,x)
            xnmax = x[numpy.argmax(y)]
            self.TTmax = self.TTxn(xnmax)
            
        xn = self.xn(data['rank'])   
        scales = {k: numpy.mean(data[k] / self.predict[k](xn)) for k in self.predict if k in data.columns}
        for k in scales:
            self.predict[k] = interp1d(self.xnref, numpy.array(self.ref[k].tolist()) * scales[k])
        return scales
    
    def xn(self, ranks, nff=None):
        """
        xn = (TTem_leaf - TTem_leafmax) / mean_phyllochron
        """
        return (self.hsfit.TT(ranks, nff) - self.TTmax) * self.hsfit.a_cohort
        
    def TTxn(self, xn, nff=None):
        return self.TTmax + xn / self.hsfit.a_cohort
    
    def dimT(self, nff=None, scale=1.0):
    
        if not isinstance(scale,dict):
            scale={}
        
        if nff is None:
            nff = self.hsfit.nff
            ranks = numpy.linspace(0, nff, 20)
        else:
            ranks = numpy.arange(1, nff + 1)
            
          
        xn = self.xn(ranks,nff)        
        df = pandas.DataFrame({k:self.predict[k](xn) * scale.get(k,1) for k in self.predict})
        df['rank'] = ranks
        return df

    def dimT_user_table(self, nff=12, scale=1.0):
        df = self.dimT(nff,scale)
        df = df.rename(columns={'rank':'index_phytomer'})
        df = df.drop('H_col',axis=1)
        return df  

class AxePop(object):
    """ An axe population generator based on plantgen axe generation methods and new extension (smart pop/damages)
    """
        
    def __init__(self, MS_leaves_number_probabilities={'11': 0.3, '12': 0.7}, Emission=None, Regression=None, tiller_damages=None, max_order=None, std_em=0):

        if Regression is None:
            Regression = TillerRegression(ears_per_plant=2.5, n_elongated_internode=4, delta_stop_del= 2)            
        if Emission is None:
            Emission = TillerEmission({'T1':0.95, 'T2':0.85, 'T3':0.75, 'T4':0.4})
            
        self.Emission = Emission
        self.Regression = Regression
        self.tiller_damages = tiller_damages
        self.MS_probabilities = MS_leaves_number_probabilities
        self.max_order = max_order
        self.std_em = std_em
          
    def mean_nff(self):
        return sum([int(k) * v for k,v in self.MS_probabilities.iteritems()])
        
    def mode_nff(self):
        v = self.MS_probabilities.values()
        nff = self.MS_probabilities.keys()
        return int(nff[v.index(max(v))])

    def hs_debreg(self):
        nff = self.mean_nff()
        return self.Regression.hs_debreg(nff)
        
    def axis_emission_table(self):
        nff = self.mean_nff()
        #regression start is set/defined at population level
        hs_debreg = self.Regression.hs_debreg(nff)
        table = self.Emission.emission_table(nff, hs_debreg=hs_debreg, max_order=self.max_order)
        return table

    def cohort_emission_table(self):
        axis_table = self.axis_emission_table()
        cohort_table = self.Emission.cohort_emission_table(axis_table)
        return cohort_table
        
    def curve_emission_table(self):
        axis_table = self.axis_emission_table()
        return self.Emission.curve_emission_table(axis_table)
        
    def regression_table(self):
        cohort_table = self.cohort_emission_table()
        return self.Regression.regression_table(cohort_table)
     
    def damage_table(self):
        """ compute damage table for the differents cohorts
        """
        damages = self.tiller_damages
        if damages is not None:
            when_start,when_end = damages['when']
            f_damaged = tools.calculate_decide_child_cohort_probabilities(damages['damage'])# tools function convert 'tiller name' kays into cohort index keys
            damages = pandas.DataFrame({'cohort':f_damaged.keys(), 'f_damaged':f_damaged.values()})
            # compute start/end for the different cohorts
            delays = pandas.DataFrame({'cohort':self.Emission.cohort_delays().keys(), 'delay':self.Emission.cohort_delays().values()})
            damages = damages.merge(delays)
            damages['start_damages'] = numpy.maximum(damages['delay'],when_start)
            damages['end_damages'] = numpy.maximum(damages['delay'], when_end)
            #filter tillers emited after end
            damages = damages.loc[damages['end_damages'] > damages['start_damages'],:]
            damages = damages.drop('delay',axis=1)
        return damages
     
    def emission_curves(self, include_MS=True):
        axis_table = self.axis_emission_table()
        return self.Emission.emission_curves(axis_table, include_MS=include_MS)
        
    def regression_curves(self):
        return self.Regression.regression_curves(self.cohort_emission_table(), self.curve_emission_table(), self.damage_table())
     
    def damage_curves(self):
        """ interpolation function for damages to tilers (number of axes lost) as a function of haun stage
        """
        curves = {}
        cohorts = self.curve_emission_table()
        damages = self.damage_table()
        if damages is None:
            hs = [0,20]
            damages = [0,0]
            no_damage = interp1d(hs, damages)       
            for w in ('primary','other', 'total'):
                curves[w] = no_damage
        else:   
            damages = damages.merge(cohorts)
            damages = damages.set_index('cohort')
            for w in ('primary','other', 'total'):
                # compute curves for each cohort individualy
                cfits = {}
                for c in damages.index:
                    d = damages.loc[c,:]
                    card = d[w + '_axis'] * d['f_damaged']
                    hs = [0,d['start_damages'], d['end_damages'],20]
                    loss = [0] * 2 + [card] * 2
                    cfits[c] = {'hs':hs,'loss':loss}
                #merge cohorts
                hs = [0] + numpy.unique(damages.ix[:,('start_damages','end_damages')].values).tolist() + [20]
                loss = numpy.array([0] * len(hs))
                for c in cfits:
                    loss = loss + numpy.interp(hs,cfits[c]['hs'],cfits[c]['loss'])
                curves[w] = interp1d(hs, loss)              
        return curves
     
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
        
    def axis_dynamics(self, plant_density = 1, include_MS = True):
        """ Compute axis density table (growing + stopped but not disapeared) = f (HS_mean_MS)
        
            Parameters:
                - plant_density : plant per square meter or callable returning plant density as a function of haun stage
        """
        hs = numpy.arange(0,18,0.1)
        
        if not callable(plant_density):
            plant_density = interp1d([0,20], [plant_density] * 2)
            
        emission = self.emission_curves(include_MS=include_MS)
        regression = self.regression_curves()
        damages = self.damage_curves()
        dynamics = {'HS':hs}
        for w in ('primary','other', 'total'):
            dynamics[w] = plant_density(hs) * (emission[w](hs) - regression[w](hs) - damages[w](hs))

        #approx nt3F : max = nb axe > 2phyllo at start of regression (stop of dvpt)
        hs_debreg = self.hs_debreg()
        max3F = numpy.interp(hs_debreg - 2, hs, dynamics['total'])
        tt3F = dynamics['total'].copy()
        tt3F[tt3F > max3F] = max3F
        tt3F_em = numpy.interp(hs, hs + 2, tt3F)
        tt3F[hs < hs_debreg + 2] = tt3F_em[hs < hs_debreg + 2]
        dynamics['3F'] = tt3F
        
        return pandas.DataFrame(dynamics)


    def smart_population(self, nplants=1):
        """ Generate an axe population using deteministic rounding
        """
        table = self.axis_emission_table()              
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
        #count cohort cardinalities per nff
        def _dfc(x,y):
            df = pandas.DataFrame(x)
            df = df.rename(columns={0:'cohort',1:'axe'})
            df['nff'] = y
            return df
        df = pandas.concat(map(lambda x: _dfc(x[0],x[1]), zip(plant_axes, list(reversed(nff_MS)))))# plants are constructed by poping nff_MS, ie from last to first
        cardnff = df.groupby(('nff','cohort')).agg('count')
        #find nff per cohort
        cohort_decimal_nff = {int(k):self.Emission.final_leaf_numbers(int(k)) for k in self.MS_probabilities}
        cohort_nff_modalities = {k:{kk:modalities(vv) for kk,vv in v.iteritems()} for k,v in cohort_decimal_nff.iteritems()}
        cohort_nff_cardinalities = {}
        for nff in nff_MS_cardinalities:
            d = {}
            for c in cardnff.loc[int(nff)].index:
                d[int(c)] = cardinalities(cohort_nff_modalities[int(nff)][int(c)], int(cardnff.loc[int(nff),c]))
            cohort_nff_cardinalities[int(nff)] = d
        cohort_nff = {k:{kk:card2list(vv) for kk,vv in v.iteritems()} for k,v in cohort_nff_cardinalities.iteritems()}
        
        plants = []
        TTem = get_normal_dist(nb_plants=nplants, sigma=self.std_em)
        for i in range(nplants): 
            id_plant = i + 1
            axes =  plant_axes[i]
            id_cohort, id_axis = zip(*axes)
            # deterministic sampling of nff for main stem 
            nff_p = nff_MS.pop()
            # nfff on tillers
            nff = [cohort_nff[nff_p][c].pop() for c in id_cohort]
            plants.append(pandas.DataFrame({'id_plt': id_plant, 'id_cohort': id_cohort, 'id_axis' : id_axis, 'N_phytomer_potential': nff, 'TTem':TTem[i]}))
        df=pandas.concat(plants)
        df= df.sort(['id_plt','id_cohort','id_axis'])
        return df

    def disparition_times(self, population):
        """ Compute disparition times  of a population of axis, expressed in haun stage of mean plant
        """
              
        regression_table = self.regression_table()
        regression_table.set_index('cohort',inplace=True)
        fdisp = regression_table['f_disp'].to_dict()
        
        damages = self.damage_table()
        
        cards = population.groupby('id_cohort').count()['id_axis'].to_dict()
        
        nreg = {k: round(fdisp[k] * cards[k]) for k in fdisp if k in cards}
        nreg = {k:v for k,v in nreg.iteritems() if v > 0}
        treg = {k: t_death(v, regression_table['t_start'][k], regression_table['t_disp'][k]) for k,v in nreg.iteritems()}
        
        if damages is None:
            tdisp = treg
        else:
            damages.set_index('cohort', inplace=True)
            fdamaged = damages['f_damaged'].to_dict()
            ndamaged = {k: round(v * cards[k]) for k,v in fdamaged.iteritems()}
            ndamaged = {k:v for k,v in ndamaged.iteritems() if v > 0}
            assert sum(ndamaged.values()) <= sum(nreg.values()), 'Damages are too important to be compensated by reggressing tillers !'
            tdamaged = {k: t_death(v, damages['start_damages'][k], damages['end_damages'][k]) for k,v in ndamaged.iteritems()}
            tdisp = {}
            # for damaged regressing tillers, choose min(tdamage, treg)
            for k in tdamaged:
                if k in treg:
                    tdisp[k] = []
                    for i in range(min(len(treg[k]), len(tdamaged[k]))):
                        tdisp[k].append(min(treg[k].pop(), tdamaged[k].pop()))
            # remaining damages are compensated by unregressing oldest tillers (compensation is for getting the right nears/plant)
            for k in tdamaged:
                if len(tdamaged[k]) > 0:# more than the regressing tillers are damaged
                    if not k in tdisp:
                        tdisp[k]= []
                    tdisp[k].extend(tdamaged[k])
                    # compensation : unregressing oldest tillers
                    for i in range(len(tdamaged[k])):
                        for j in sorted(treg.keys()):
                            if len(treg[j]) > 0:
                                treg[j].pop()
                                break
            # add remaining regression
            for k in treg:
                if len(treg[k]) > 0:
                    if not k in tdisp:
                        tdisp[k] = []
                    tdisp[k].extend(treg[k])
        
        disp_times = {k: tdisp[k] + [None] * (cards[k] - len(tdisp[k])) if k in tdisp else [None] * cards[k] for k in cards}
        for k in disp_times:
            random.shuffle(disp_times[k])
        
        return [disp_times[c].pop() for c in population['id_cohort']]
        
    def plant_list(self, nplants=1):
        pop = self.smart_population(nplants)
        pop['hs_disparition'] = self.disparition_times(pop)
        pop['hs_stop'] = pop['hs_disparition'] - self.Regression.delta_stop_del
        return dict(list(pop.groupby('id_plt'))).values()
        

        
class PlantGen(object):
    """ A class interface to plantgen
    """
    
    def __init__(self, HSfit=None, GLfit=None, Dimfit=None, inner_parameters={}):
    
        #defaults for fitted models
        if HSfit is None:#HS fit is for the main stem of the plant (or mean plant but then difference for flag leaf emergence should be added in equations)
            HSfit = HaunStage()           
        if GLfit is None:
            GLfit = GreenLeaves()        
        if Dimfit is None:
            Dimfit = WheatDimensions()            
        self.inner_parameters = inner_parameters
        self.HSfit = HSfit
        self.GLfit = GLfit
        self.Dimfit = Dimfit
        
        #setup base configuration for plantgen interface for one plant
        self.base_config = {}
        self.base_config['inner_params'] = inner_parameters
        self.base_config['inner_params']['EMF_1_MS_STANDARD_DEVIATION'] = 0 #strictly respect TT_hs_0. 
        self.base_config.update({'plants_number': 1,
                            'plants_density': 1.})
                  
    
     
    def config(self, plant):
        """ Generate plantgen configuration dict for one plant without any axis regression
        """
           
        config = {}
        config.update(self.base_config) #avoid side effects on base_config
        
        nff = plant['N_phytomer_potential'].max()
               
        # Tillering
          
        axeT_user = self.axeT_user_table(plant)
        config.update({'ears_density' : None, #no regression 
                        'axeT_user':axeT_user,
                        'MS_leaves_number_probabilities':{str(nff):1.0},
                        'decide_child_axis_probabilities':{k:1.0 for k in plant['id_axis'] if  0 < _order(k) <= 1}
                        })
                       
        # Dimensions
        config['dimT_user'] = self.Dimfit.dimT_user_table(nff)
        
        # Leaf appearance / Green leaves          
        config['dynT_user'] = self.dynT_user_table(plant)
        
        hs_flag = nff
        GL = self.GLfit.HS_GL_sample(hs_flag)
        GL['TT'] = self.HSfit.TT(GL['HS'], nff)
        GL = dict(zip(GL['TT'],GL['GL']))
        config['GL_number'] = GL
        
        hs_t1 = self.GLfit.hs_t1(hs_flag)
        config['TT_t1_user'] = self.HSfit.TT(hs_t1, nff)
        
        return config

    def pgen_tables(self, plant):
        config = self.config(plant)
        axeT_, dimT_, phenT_, phenT_abs, dynT_, phenT_first, HS_GL_SSI_T, tilleringT, cardinalityT, config = gen_adel_input_data(**config)
        axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
        # include regression
        axeT['TT_stop_axis'] = self.HSfit.TT(plant['hs_stop'])# use hs of mean plant
        axeT['TT_del_axis'] = self.HSfit.TT(plant['hs_disparition'])
        for w in ('TT_em_phytomer1','TT_col_phytomer1','TT_sen_phytomer1','TT_del_phytomer1'):
            axeT[w] = axeT[w] + plant['TTem']
        for i in axeT.index:
            TT_stop = axeT['TT_stop_axis'][i]
            TTem = plant['TTem'].values[0]
            if not numpy.isnan(TT_stop):
                phen = phenT_abs[phenT_abs['id_phen'] == axeT['id_phen'][i]]
                axeT['HS_final'][i] = numpy.interp(TT_stop + TTem, phen['TT_col_phytomer'], phen['index_phytomer'])
                axeT['id_ear'][i] = numpy.nan
                
        # include plant number in ids to allow future concatenation with other plants
        id_plt = plant['id_plt'][0]
        axeT['id_dim'] = id_plt * 1e5 + axeT['id_dim']
        dimT['id_dim'] = id_plt * 1e5 + dimT['id_dim']
        axeT['id_phen'] = id_plt * 1e5 + axeT['id_phen']
        phenT['id_phen'] = id_plt * 1e5 + phenT['id_phen']
        
        return {'adelT': (axeT, dimT, phenT), 'phenT_abs':phenT_abs, 'phenT_first':phenT_first, 'HS_GL_SSI_T':HS_GL_SSI_T, 'tilleringT':tilleringT, 'cardinalityT':cardinalityT, 'config':config}
        
    def adelT(self, plants):
        """ Compute Adel input tables (axeT, phenT, dimT) for a collection of plants
        """        
        tables = map(lambda x: self.pgen_tables(x)['adelT'], plants)
        if len(plants) > 1:
            axeT, dimT, phenT =  map(pandas.concat, zip(*tables))
        else:
            axeT, dimT, phenT = tables[0]
        return axeT, dimT, phenT
     
    def axeT_user_table(self, axeT):
    
        df = pandas.DataFrame(index=range(len(axeT['id_plt'])),
                                               columns=['id_plt', 'id_cohort', 'id_axis', 'N_phytomer_potential', 'N_phytomer', 'HS_final', 'TT_stop_axis', 'TT_del_axis', 'id_dim', 'id_phen', 'id_ear', 'TT_app_phytomer1', 'TT_col_phytomer1', 'TT_sen_phytomer1', 'TT_del_phytomer1'],
                                               dtype=float)

        df['id_plt'] = axeT['id_plt']
        df['id_axis'] = axeT['id_axis']
        df['id_cohort'] = axeT['id_cohort']
        df['N_phytomer_potential'] = axeT['N_phytomer_potential']
        df['id_phen'] = (df['id_cohort'] * 100 + df['N_phytomer_potential']) * 10 + 1#id for axe with ear
        
        df= df.sort(['id_plt','id_cohort','id_axis'])
        return df
      
    def dynT_user_table(self, plant):
    
        nff = plant['N_phytomer_potential'].max()
        MS_parameters = {'a_cohort': self.HSfit.a_nff(nff),
                         'TT_hs_0': self.HSfit.TT_hs_0,
                         'TT_hs_N_phytomer_potential': self.HSfit.TTflag(nff),
                         'n0': self.GLfit.n0,
                         'n1': self.GLfit.n1,
                         'n2': self.GLfit.n2}
        idaxis = plant['id_axis'].values
        df = pandas.DataFrame(index=idaxis,
                              columns=['id_axis','a_cohort','TT_hs_0','TT_hs_N_phytomer_potential','n0','n1','n2'],
                              dtype=float)
        df.ix['MS'] = pandas.Series(MS_parameters)
        df['id_axis'] = idaxis
        cohort = plant['id_cohort'].values
        dTT_MS_cohort = numpy.frompyfunc(self.HSfit.dTT_MS_cohort, 1, 1) # make function operates on 'vectorised' arrays
        df['TT_hs_N_phytomer_potential'] = self.HSfit.TTflag(nff) + dTT_MS_cohort(cohort)
        df = df.reset_index(drop=True)
        return df
    
    
 

    

# extensions
 
def adelT_to_devT(pgen):
    """ Creates devT tables from plantgen dict
    """
    devT = devCsv(*pgen['adelT'])
    return devT, 



  
   
#deprecated

 
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
#
#       Adel.WheatTillering
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

""" Wraping /extending  Wheat Tillering model used in Alinea.plantgen for calibration purposes
    Author : Christian Fournier
"""
import pandas
import numpy

from alinea.adel.plantgen import params, tools
import alinea.adel.plantgen_extensions as pgen_ext

#
# New propositions (C. Fournier, January 2014)
#
#index cohort delays with cohorts number (instead of primary_axis id of the same cohort) and add MS
cohort_delays = {(int(k.lstrip('T')) + 3):v for k,v in params.LEAF_NUMBER_DELAY_MS_COHORT.iteritems()}
cohort_delays[1] = 0


#
# Add  defaults to handle 'less than MIN' datasets / set some defaults differents from those of plantgen
#
# Defaults probabilities of appearance of tillers
primary_tiller_probabilities = {'T0':0,'T1':0.95, 'T2':0.85, 'T3':0.75, 'T4':0.4, 'T5':0.2, 'T6':0.1}
# relative probability  of emergence (ie actual to theoretical) of secondary tillers for the different cohorts (may be usfull for fitting?)
#secondary_tiller_probabilities = {i:1.0 for i in range(11)}
# Express time-delta between tiller regression and bolting in phyllochronic units
delta_reg = float(params.DELAIS_REG_MONT) / 110
# number of elongated internodes (used to compute hs_regression_start with the same formula as pgen: hs_start_reg =  hs_bolting + delta_reg  & hsbolting = hs_end - n_elongated_internode
n_elongated_internode = 4
# mean number of ears per plant. This is better than ear density as it allows separate fitting of tillering and of and global scaling of axis per plant  with ear density
ears_per_plant = 2.5
# mean number of leaves on main stem
nff = 12.0
# delay between end of growth and disparition of an axe
delta_stop_del = 2.


def decimal_elongated_internode_number(internode_ranks, internode_lengths):
    """ estimate the (decimal) number of elongated internode as determined in Plantgen from internode dimension of most frequent axis
    """
    df = pandas.DataFrame({'rank':internode_ranks, 'length':internode_lengths})
    # keep only the non-zero lengths
    df = df[df['length'] > 0]
    # Fit a polynomial of degree 2 to, and get the coefficient of degree 0. 
    n = df['rank'].max() - numpy.polyfit(df['length'].values, df['rank'].values, 2)[2]
    return n

  

class WheatTillering(object):
    """Model of tiller population dynamics on wheat"""
    
    def __init__(self, primary_tiller_probabilities = primary_tiller_probabilities,
                       ears_per_plant = ears_per_plant,
                       nff = nff,
                       n_elongated_internode = n_elongated_internode,
                       child_cohort_delay = params.FIRST_CHILD_DELAY,
                       #secondary_tiller_probabilities =secondary_tiller_probabilities,
                       cohort_delays = cohort_delays, #delays HS_0 MS -> HS_0 tillers HS=0 <=> emergence leaf 1
                       delta_reg = delta_reg,
                       a1_a2 = params.SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS,
                       delta_stop_del = delta_stop_del,
                       max_order = None,
                       tiller_survival=None):
        """ Instantiate model with default parameters """
        self.primary_tiller_probabilities = primary_tiller_probabilities
        self.ears_per_plant = ears_per_plant
        self.nff = nff
        self.n_elongated_internode = n_elongated_internode
        self.child_cohort_delay = child_cohort_delay
        #self.secondary_tiller_probabilities = secondary_tiller_probabilities
        self.cohort_delays = cohort_delays
        self.delta_reg = delta_reg
        self.a1_a2 = a1_a2
        self.delta_stop_del=delta_stop_del
        self.max_order = max_order
        self.tiller_survival = tiller_survival #parameter used to simulate plant damages (fly, freeze...) that make tillers die earlier than expected with the tiller regression model

        
       
    def theoretical_probabilities(self):
        """ Computes cohort number, botanical position and theoretical probabilities of emergence all tillers using botanical position and probabilities of emergence of primary ones
        """

        axis_probabilities = {k:v for k,v in self.primary_tiller_probabilities.iteritems() if v > 0 }
        cohort_probabilities = tools.calculate_decide_child_cohort_probabilities(axis_probabilities)
        _,res = tools.calculate_theoretical_cardinalities(1, 
                                                      cohort_probabilities,
                                                      axis_probabilities,
                                                      self.child_cohort_delay)
        return res
         
         
    def final_leaf_numbers(self):
        """ returns the mean final leaf numbers of cohorts
        """
        axis_probabilities = {k:v for k,v in self.primary_tiller_probabilities.iteritems() if v > 0 }
        cohort_probabilities = tools.calculate_decide_child_cohort_probabilities(axis_probabilities)
        ms_nff = [(1,self.nff)]
        cohort_nff = [(cohort,tools.calculate_tiller_final_leaves_number(self.nff, cohort, self.a1_a2)) for cohort in cohort_probabilities]
        return dict(ms_nff + cohort_nff)   
        
    def emited_cohort_density(self, plant_density = 1):
        proba = self.theoretical_probabilities()
        nffs = self.final_leaf_numbers()
        data = zip(*[(k[0],k[1],v) for k,v in proba.iteritems()])
        df = pandas.DataFrame({'axis': data[1],
                               'cohort' : data[0] ,
                               'probability':data[2], 
                               'primary':map(lambda x: not '.' in x, data[1])})
        if self.max_order != None:
            def _order(axe):
                return len(axe.rsplit('.'))
            df = df[map(lambda x: _order(x) <= self.max_order,df['axis'])]
 
        grouped = df.groupby(['cohort'],group_keys=False)
        def _fun(x):
            cohort = x['cohort'].values[0]
            d = {'cohort': cohort, 
                 'delay': self.cohort_delays[cohort],
                 'primary_axis':x[x['primary']]['probability'].sum() * plant_density,
                 'other_axis' : x[~x['primary']]['probability'].sum() * plant_density   
                 #'other_axis' : x[~x['primary']]['probability'].sum() * plant_density * self.secondary_tiller_probabilities[cohort],
                 }
            out = pandas.DataFrame(d,index = [x.index[0]])
            out['total_axis'] = out['primary_axis'] + out['other_axis']
            out['nff'] = nffs[cohort]
            return out
        return grouped.apply(_fun)
        
    def axis_list(self, nplants = 2):
        """ compute cardinalities of axis in a stand of n plants
        The strategy used here is based on deterministic rounding, and differs from the one used in plantgen (basd on random sampling). Difference are expected for small plants numbers
        """
        cohorts = self.emited_cohort_density()
        # filter, if any, cohorts emited after regression start (not realy handled in pgen, as no model of regression is available for these tillers)
        cohorts = cohorts[cohorts['delay'] <= self.hs_debreg()] 
        
        axis=pgen_ext.axis_list(cohorts, self.theoretical_probabilities(), nplants)
                    
        return axis
        
    def plant_list(self, nplants = 2):

        axis = self.axis_list(nplants)
        plants = pgen_ext.plant_list(axis, nplants)
        return plants
        
    def to_pgen(self, nplants=2, density = 250, phyllochron = 110, TTem=0, pgen_base ={}):
        plants = self.plant_list(nplants)
        axeT = pgen_ext.axeT_user(plants)
        mods = pgen_ext.modalities(self.nff)
        nff_plants = {k: axeT['id_plt'][(axeT['id_axis'] == 'MS') & (axeT['N_phytomer_potential'] == k)].values.astype(int).tolist() for k in mods}
        
        base={}
        base.update(pgen_base)# avoid altering pgen_base
        base.update({'decide_child_axis_probabilities' : self.primary_tiller_probabilities,
              'plants_density': density,
              'ears_density' : density * self.ears_per_plant, 
              'delais_TT_stop_del_axis': self.delta_stop_del * phyllochron,
              'TT_regression_start_user': TTem + self.hs_debreg() * phyllochron
              })
        
        pgen = {k:{'MS_leaves_number_probabilities': {str(k):1.0},
                   'axeT_user':axeT.ix[axeT['id_plt'].isin(v),:],
                   'plants_number':len(v)} for k,v in nff_plants.iteritems() if len(v) > 0}
        
        for k in pgen:
            pgen[k].update(base)
            
        return pgen
    
    def hs_debreg(self):
        hs_bolting = self.nff - self.n_elongated_internode
        return hs_bolting + self.delta_reg
        
    def cohort_survival(self):
        cohort_survival = None
        if self.tiller_survival is not None:
            cohort_survival = tools.calculate_decide_child_cohort_probabilities(self.tiller_survival)# tools function convert 'tiller name' kays into cohort index keys
        return cohort_survival
    
    def axis_dynamics(self, plant_density = 1, hs_bolting = None, include_MS = True):
        """ Compute axis density (growing + stopped but not deleted) = f (HS_mean_MS)
            Parameters:
                - plant_density : plant per square meter
                - hs_bolting : Force Haun Stage of mean_MS at bolting (facultative).If None, hs_bolting is estimated using the number of elongated internodes
        """
        
        hs_max = self.nff
        hs_debreg = self.hs_debreg() 
        
        cohorts = self.emited_cohort_density(plant_density = plant_density)
        ear_density = self.ears_per_plant * plant_density
                
        # filter, if any, cohorts emited after regression start
        cohorts = cohorts[cohorts['delay'] <= hs_debreg]
        
        # early loss (if any) occuring before debreg
        early_loss = 0
        early_frac = 0
        when_lost = -1 #means never
        if self.tiller_survival is not None:
            cohort_survival = self.cohort_survival()
            start_loss = {k:max(float(cohorts['delay'][cohorts['cohort'] == k]),v['HS'][0]) for k,v in cohort_survival.iteritems()}
            start_loss = {k:v for k,v in start_loss.iteritems() if v < hs_debreg}
            lost = cohorts['cohort'].isin(start_loss)
            # apply early loss at weighted mean of losses (to be replaced by differentila rate/date of loss for every cohorts)
            start = sum([start_loss[k] * float(cohorts['total_axis'][cohorts['cohort']==k]) for k in start_loss]) / cohorts['total_axis'][lost].sum()
            when_lost = numpy.mean((start,hs_debreg))
            # ammount lost before startreg
            fraction_lost = {k:(1 - float(cohort_survival[k]['density'][1])) for k in start_loss}
            end_loss = {k:float(cohort_survival[k]['HS'][1]) for k in start_loss}
            total_lost = sum([fraction_lost[k] * float(cohorts['total_axis'][cohorts['cohort']==k]) for k in start_loss])
            end = sum([end_loss[k] * float(cohorts['total_axis'][cohorts['cohort']==k]) for k in start_loss]) / cohorts['total_axis'][lost].sum()
            early_loss = total_lost * min(1,(hs_debreg - start) / (end - start))
            early_frac = early_loss / cohorts['total_axis'].sum()
            
        #regression rate (all cohorts)
        dmax = cohorts['total_axis'].sum() - early_loss
        regression_rate = (ear_density - dmax) / (hs_max - hs_debreg)
        # compute loss values that correspond to 100 percent loss for a cohort
        reverse_cohorts = cohorts.sort_index(by=['delay'], ascending = False)
        cohorts['cumulative_loss'] = reverse_cohorts['total_axis'].cumsum()       

        
        def _density(x, delays, cardinalities, total, cumulative_loss):
            if x < when_lost:#when_lost is time of deletion
                d = cardinalities[delays <= x].sum()
            elif x < (hs_debreg + self.delta_stop_del):#hs_debreg is time of stop
                d = cardinalities[delays <= x].sum() - early_frac * cardinalities.sum()
            else :
                hs = min(x,hs_max + self.delta_stop_del)
                loss = - regression_rate * (hs - hs_debreg - self.delta_stop_del)#dmax - (dmax + reg)
                fraction_lost = (loss - cumulative_loss + total) / total
                fraction_lost = numpy.maximum(0,numpy.minimum(1,fraction_lost))
                d = (cardinalities * (1 - fraction_lost)).sum() - early_frac * cardinalities.sum()
            return d
               
        hs = numpy.arange(0,1.2 * (hs_max + self.delta_stop_del),0.1)
        primary = numpy.array(map(lambda x: _density(x, cohorts['delay'], cohorts['primary_axis'], cohorts['total_axis'],cohorts['cumulative_loss']),hs))        
        others = numpy.array( map(lambda x: _density(x, cohorts['delay'], cohorts['other_axis'], cohorts['total_axis'],cohorts['cumulative_loss']),hs))
        total = primary + others
        
        max3F = numpy.interp(hs_debreg - 2, hs, total)
        tt3F = total.copy()
        tt3F[tt3F > max3F] = max3F
        tt3F_em = numpy.interp(hs, hs + 2, tt3F)
        tt3F[hs < hs_debreg + 2] = tt3F_em[hs < hs_debreg + 2]
        
        df = pandas.DataFrame({'HS':hs, 'primary':primary, 'others' : others, 'total': total, '3F':tt3F})
        if not include_MS:
            for w in ['primary', 'total', '3F']:
                df[w] -= 1
        return df
        

        
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
# Express time-delta between tiller regression and bolting in phyllochronic units
delta_reg = float(params.DELAIS_REG_MONT) / 110

#
# Add  defaults to handle 'less than MIN' datasets
#
# Defaults probabilities of appearance of tillers
primary_tiller_probabilities = {'T0':0,'T1':0.95, 'T2':0.85, 'T3':0.75, 'T4':0.4, 'T5':0.2, 'T6':0.1}
# relative probability  of emergence (ie actual to theoretical) of secondary tillers for the different cohorts (may be usfull for fitting?)
#secondary_tiller_probabilities = {i:1.0 for i in range(11)}
# number of elongated internodes (used in case hs_bolting is unknown, hs_bolting = hs_end - n_elongated_internode)
n_elongated_internode = 4
# mean number of ears per plant. This is better than ear density as it allows separate fitting of tillering and of and global scaling of axis per plant  with ear density
ears_per_plant = 2.5
# mean number of leaves on main stem
nff = 12.0

# converter to pgen
#


class WheatTillering(object):
    """Model of tiller population dynamics on wheat"""
    
    def __init__(self, primary_tiller_probabilities = primary_tiller_probabilities,
                       ears_per_plant = ears_per_plant,
                       nff = nff,
                       n_elongated_internode = n_elongated_internode,                   child_cohort_delay = params.FIRST_CHILD_DELAY,
                       #secondary_tiller_probabilities =secondary_tiller_probabilities,
                       cohort_delays = cohort_delays, #delays HS_0 MS -> HS_0 tillers HS=0 <=> emergence leaf 1
                       delta_reg = delta_reg,
                       a1_a2 = params.SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS, delta_stop_del = 600. / 110,
                       max_order = None):
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
        #plantgen copies of primary proba emission per axis (0 filtered) and per cohort
        self.axis_probabilities = {k:v for k,v in self.primary_tiller_probabilities.iteritems() if v > 0 }
        self.cohort_probabilities = tools.calculate_decide_child_cohort_probabilities(self.axis_probabilities)
        
       
    def theoretical_probabilities(self):
        """ Computes cohort number, botanical position and theoretical probabilities of emergence all tillers using botanical position and probabilities of emergence of primary ones
        """
        _,res = tools.calculate_theoretical_cardinalities(1, 
                                                      self.cohort_probabilities,
                                                      self.axis_probabilities,
                                                      self.child_cohort_delay)
        return res
         
         
    def final_leaf_numbers(self):
        """ returns the mean final leaf numbers of cohorts
        """
        ms_nff = [(1,self.nff)]
        cohort_nff = [(cohort,tools.calculate_tiller_final_leaves_number(self.nff, cohort, self.a1_a2)) for cohort in self.cohort_probabilities]
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
        
    def to_pgen(self, nplants=2, density = 250, force_start_reg = False):
        plants = self.plant_list(nplants)
        axeT = pgen_ext.axeT_user(plants)
        mods = pgen_ext.modalities(self.nff)
        nff_plants = {k: axeT['id_plt'][(axeT['id_axis'] == 'MS') & (axeT['N_phytomer_potential'] == k)].values.astype(int).tolist() for k in mods}
        
        pgen_base = {'decide_child_axis_probabilities' : self.primary_tiller_probabilities,
              'plants_density': density,
              'ears_density' : density * self.ears_per_plant,
              'TT_col_break' : 0.0
              }
        if force_start_reg:
            pgen_base.update({'hs_deb_reg': self.hs_debreg()})
        pgen = {k:{'MS_leaves_number_probabilities': {str(k):1.0},
                   'axeT_user':axeT.ix[axeT['id_plt'].isin(v),:],
                   'plants_number':len(v)} for k,v in nff_plants.iteritems() if len(v) > 0}
        
        for k in pgen:
            pgen[k].update(pgen_base)
            
        return pgen
    
    def hs_debreg(self):
        hs_bolting = self.nff - self.n_elongated_internode
        return hs_bolting + self.delta_reg
    
    def axis_dynamics(self, plant_density = 1, hs_bolting = None, include_MS = True):
        """ Compute axis density = f (HS_mean_MS)
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
        # compute loss values that correspond to 100 percent loss for a cohort
        reverse_cohorts = cohorts.sort_index(by=['delay'], ascending = False)
        cohorts['cumulative_loss'] = reverse_cohorts['total_axis'].cumsum()       
        dmax = cohorts['total_axis'].sum()
        regression_rate = (ear_density - dmax) / (hs_max - hs_debreg)
        
        def _density(x, delays, cardinalities, total, cumulative_loss):
            if x < (hs_debreg + self.delta_stop_del):
                d = cardinalities[delays <= x].sum()
            else :
                hs = min(x,hs_max + self.delta_stop_del)
                loss = - regression_rate * (hs - hs_debreg - self.delta_stop_del) #dmax - (dmax + reg)
                fraction_lost = (loss - cumulative_loss + total) / total
                fraction_lost = numpy.maximum(0,numpy.minimum(1,fraction_lost))
                d = (cardinalities * (1 - fraction_lost)).sum()           
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
        

        
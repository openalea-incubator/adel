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
from scipy.interpolate import interp1d

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
                       tiller_damages=None,
                       hs_bolting=None):
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
        self.tiller_damages = tiller_damages #parameter used to simulate plant damages (fly, freeze...) that make tillers die earlier than expected with the tiller regression model
        if hs_bolting is None:
            self.hs_bolting = self.nff - self.n_elongated_internode
        else:
            self.hs_bolting = hs_bolting
        
       
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
        
    def emited_cohort_density(self, plant_density=1, include_MS=True):
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
        cohorts = grouped.apply(_fun)
        
        # filter, if any, cohorts emited after regression start (not realy handled in pgen, as no model of regression is available for these tillers)
        cohorts = cohorts[cohorts['delay'] <= self.hs_debreg()]
        
        #filter MS if needed
        if not include_MS:
            cohorts = cohorts.loc[cohorts['cohort'] > 1,:] 
            
        return cohorts
        
    def axis_list(self, nplants = 2):
        """ compute cardinalities of axis in a stand of n plants
        The strategy used here is based on deterministic rounding, and differs from the one used in plantgen (basd on random sampling). Difference are expected for small plants numbers
        """
        cohorts = self.emited_cohort_density()
  
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
        return self.hs_bolting + self.delta_reg
        
    def cohort_survival(self):# futur deprecated
        cohort_survival = None
        if self.tiller_survival is not None:
            cohort_survival = tools.calculate_decide_child_cohort_probabilities(self.tiller_survival)# tools function convert 'tiller name' kays into cohort index keys
        return cohort_survival
        
    def cohort_damages(self, plant_density=1):#futur deprecated
        damages = self.tiller_damages
        if damages is not None:
            when = damages['when']
            f_damaged = tools.calculate_decide_child_cohort_probabilities(damages['damage'])# tools function convert 'tiller name' kays into cohort index keys
            hs_deb, hs_end, reg_rate, reg_table = self.regression_parameters()
            cohorts = self.emited_cohort_density(plant_density = plant_density)
            cohorts = cohorts.merge(reg_table)
            damages = cohorts[cohorts['cohort'].isin(f_damaged)]
            damages['start_damages'] = numpy.maximum(damages['delay'],when[0])
            damages['end_damages'] = numpy.maximum(damages['delay'],when[1])
            damages = damages.merge(pandas.DataFrame({'cohort':f_damaged.keys(), 'f_damaged':f_damaged.values()}))
            #filter tillers emited after end
            damages = damages.loc[damages['end_damages'] > damages['start_damages'],:]
        return damages
        
    def damage_table(self):
        """ compute damage table for the differents cohorts
        """
        damages = self.tiller_damages
        if damages is not None:
            when_start,when_end = damages['when']
            f_damaged = tools.calculate_decide_child_cohort_probabilities(damages['damage'])# tools function convert 'tiller name' kays into cohort index keys
            damages = pandas.DataFrame({'cohort':f_damaged.keys(), 'f_damaged':f_damaged.values()})
            # compute start/end for the different cohorts
            delays = pandas.DataFrame({'cohort':self.cohort_delays.keys(), 'delay':self.cohort_delays.values()})
            damages = damages.merge(delays)
            damages['start_damages'] = numpy.maximum(damages['delay'],when_start)
            damages['end_damages'] = numpy.maximum(damages['delay'], when_end)
            #filter tillers emited after end
            damages = damages.loc[damages['end_damages'] > damages['start_damages'],:]
            damages = damages.drop('delay',axis=1)
        return damages
        
    def emission_curves(self, plant_density=1, include_MS=True, delta = 30. / 110.):
        """ return interpolation functions for axe emission as a function of haun stage
        """
        cohorts = self.emited_cohort_density(plant_density = plant_density, include_MS=include_MS)
            
        hs = reduce(lambda x,y:x+y,[[hs - delta / 2,hs + delta / 2] for hs in cohorts['delay']])
        
        curves = {}
        for w in ('primary','other', 'total'):
            card = cohorts[w + '_axis'].values
            cum = card.cumsum()
            em = reduce(lambda x,y:x+y,[[cum[i] - card[i],cum[i]] for i in range(len(card))])
            curves[w] = interp1d([-delta] + hs + [20],[0] + em + [em[-1]])
        return curves
     
    def regression_table(self):
        """ compute regression table for the differents cohorts
        """
        
        # regression period for the disparition of axes
        hs_deb = self.hs_debreg() + self.delta_stop_del
        hs_end = self.nff + self.delta_stop_del
        
        #regression rate
        cohorts = self.emited_cohort_density(plant_density=1, include_MS = True)
        d_deb = cohorts['total_axis'].sum()
        d_end = self.ears_per_plant
        regression_rate = (d_end - d_deb) / (hs_end - hs_deb)
        
        # regression parameters per cohort
        # time of disparition of cohorts in case of complete regression
        # hypothesis : disparition of cohorts is sequential from last appeared to first
        reverse_cohorts = cohorts.sort_index(by=['delay'], ascending = False)
        loss_acc = reverse_cohorts['total_axis'].cumsum().values
        t_disp = hs_deb + loss_acc / abs(regression_rate)
        #actual fraction lost at t_disp
        f_disp = t_disp.copy()
        f_disp[t_disp < hs_end] = 1
        f_disp[t_disp >= hs_end] = 0
        #compute last fraction
        card = reverse_cohorts['total_axis'].values
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
               
     # ici to do: calculer damage table et en tenir compte pour mettre a jour la table de regresssion (en transferant des quantite), puis enchainer sur curve   
    def regression_curves(self, plant_density=1):
        """ interpolation function for natural axe regression (number of axes lost) as a function of haun stage
        """
        reg_table = self.regression_table()
        cohorts = self.emited_cohort_density(plant_density = plant_density)
        regressing_cohorts = reg_table.merge(cohorts)
        #
        # compensate for damages occuring before resgression starts (damages shoud not influence fitted final axe density nor final ear number)
        damages = self.damage_table()
        if damages is not None:
            # compute fraction to be compensated (ie occuring before regression starts)
            deb = min(regressing_cohorts['t_start'])
            early_frac =  1. * (deb - damages['start_damages']) / (damages['end_damages'] - damages['start_damages'])
            early_frac = 1 #all should be compensated for preserving axe density after regression
            damages['f_comp'] = damages['f_damaged'] * numpy.maximum(0,numpy.minimum(1, early_frac))
            # for damaged regressing cohort, reduce f_disp to compensate for part/all of f_comp
            regressing_cohorts = regressing_cohorts.set_index('cohort')
            damages = damages.set_index('cohort')
            cohorts = cohorts.set_index('cohort')
            for c in damages.index:
                if c in regressing_cohorts.index:
                    comp = min(damages['f_comp'][c], regressing_cohorts['f_disp'][c])
                    regressing_cohorts['f_disp'][c] = regressing_cohorts['f_disp'][c] - comp
                    damages['f_comp'][c] = damages['f_comp'][c] - comp
            # remaining compensation obtained by reducing f_disp of regressing cohorts, older first
            if sum(damages['f_comp']) > 1e-6:
                damages = damages[damages['f_comp'] > 0]
                regressing_cohorts = regressing_cohorts.sort_index(by=['delay'], ascending = True) 
                for c in damages.index:
                    f_d = damages['f_comp'][c]
                    for c_r in regressing_cohorts.index:
                        f_r = f_d * 1. * cohorts['total_axis'][c] / cohorts['total_axis'][c_r]
                        comp = min(f_r, regressing_cohorts['f_disp'][c_r])
                        regressing_cohorts['f_disp'][c_r] = regressing_cohorts['f_disp'][c_r] - comp
                        comp_d = comp * 1. * cohorts['total_axis'][c_r] / cohorts['total_axis'][c]
                        f_d = f_d - comp_d
                        if f_d < 1e-6:
                            break
                        
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
        
    # ici, to do : calculer damage table et en deduire courbe, a priori juste when a ajuster    
    def damage_curves(self, plant_density=1):
        """ interpolation function for damages to tilers (number of axes lost) as a function of haun stage
            As damages may influence regression curve, an anti-regression curve is also return, so that is total loss = regression + damage - anti_regression
        """
        #default curves for zero damages
        hs = [0,20]
        damages = [0,0]
        no_damage = interp1d(hs, damages)       
        curves = {}
        #anti_curves = {}
        for w in ('primary','other', 'total'):
            curves[w] = no_damage
            #anti_curves[w] = no_damage
        
        damages = self.damage_table()
        if damages is not None:
            cohorts = self.emited_cohort_density(plant_density = plant_density)
            damages = damages.merge(cohorts)
            damages = damages.set_index('cohort')
            for w in ('primary','other', 'total'):
                # compute curves for each cohort individualy
                cfits = {}
                for c in damages.index:
                    d = damages.loc[c,:]
                    card = d[w + '_axis']
                    hs = [0,d['start_damages'], d['end_damages'],20]
                    loss = [0] * 2 + [d['f_damaged'] * card] * 2
                    cfits[c] = {'hs':hs,'loss':loss}
                #merge cohorts
                hs = [0] + numpy.unique(damages.ix[:,('start_damages','end_damages')].values).tolist() + [20]
                loss = numpy.array([0] * len(hs))
                for c in cfits:
                    loss = loss + numpy.interp(hs,cfits[c]['hs'],cfits[c]['loss'])
                curves[w] = interp1d(hs, loss)
                #compute compensation
                # cfits = {}
                # for c in damages.index:
                    # d = damages.loc[c,:]
                    # card = d[w + '_axis']
                    # if d['f_disp'] <= 0 or d['start_damages'] >= d['t_disp']: # no compensation needed
                        # hs = [0,20]
                        # loss = [0,0]
                    # else: #compensation up to ammount lost by damages
                        # f_rec = min((d['f_damaged'],d['f_disp']))
                        # t_disp = d['t_start'] + d['f_disp'] * (d['t_disp'] - d['t_start'])
                        # t_end = d['t_start'] + f_rec * 1. /  d['f_disp'] * (t_disp - d['t_start'])
                        # hs = [0, d['t_start'], t_end, 20]
                        # loss = [0] * 2 + [f_rec * card] * 2
                    # cfits[c] = {'hs':hs,'loss':loss}
                # hs = [0] + numpy.unique(damages.ix[:,('start_damages','end_damages', 't_start','t_disp')].values).tolist() + [20]
                # loss = numpy.array([0] * len(hs))
                # for c in cfits:
                    # loss = loss + numpy.interp(hs,cfits[c]['hs'],cfits[c]['loss'])
                # anti_curves[w] = interp1d(hs, loss)                    
               
        return curves
        
    def new_axdyn(self, plant_density = 1, hs_bolting = None, include_MS = True):
        
        hs = numpy.arange(0,18,0.1)
        emission = self.emission_curves(plant_density=plant_density, include_MS=include_MS)
        regression = self.regression_curves(plant_density=plant_density)
        damages= self.damage_curves(plant_density=plant_density)
        primary = emission['primary'](hs) - regression['primary'](hs) - damages['primary'](hs)
        others = emission['other'](hs) - regression['other'](hs) - damages['other'](hs)
        total = emission['total'](hs) - regression['total'](hs) - damages['total'](hs)

        #approx nt3F : max = nb axe > 2phyllo at start of regression (stop of dvpt)
        hs_debreg = self.hs_debreg()
        max3F = numpy.interp(hs_debreg - 2, hs, total)
        tt3F = total.copy()
        tt3F[tt3F > max3F] = max3F
        tt3F_em = numpy.interp(hs, hs + 2, tt3F)
        tt3F[hs < hs_debreg + 2] = tt3F_em[hs < hs_debreg + 2]
        
        df = pandas.DataFrame({'HS':hs, 'primary':primary, 'others' : others, 'total': total, '3F':tt3F})

        return df
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
            
        #regression rate (all cohorcohortsts)
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
        

        
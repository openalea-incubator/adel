# -*- python -*-
#
#       Adel.PlantGen
#
#       Copyright 2006-2012 INRIA - CIRAD - INRA
#
#       File author(s): Camille Chambon <camille.chambon@grignon.inra.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################
'''
Generic tools used in the :mod:`alinea.adel.plantgen` package. These routines can 
also be used by other packages. 

Authors: Mariem Abichou, Camille Chambon, Bruno Andrieu
'''


import random
import math
import types

import numpy as np
import pandas
from scipy.optimize import leastsq

def decide_child_cohorts(decide_child_cohort_probabilities, parent_cohort_index=None, parent_cohort_position=None, first_child_delay=2):
    '''
    Decide (recursively) of the child cohorts actually produced by a parent cohort, 
    according to the *decide_child_cohort_probabilities* and the *parent_cohort_index*. 
    The main stem always exists.
    
    :Parameters:
    
        - `decide_child_cohort_probabilities` (:class:`dict`) - the probabilities of the 
          child cohorts.
        - `parent_cohort_index` (:class:`int`) - the index of the parent cohort. 
          ``None`` (the default) means that there isn't any parent cohort. 
        - `parent_cohort_position` (:class:`str`) - the position of the parent cohort. 
          ``None`` (the default) means that there isn't any parent cohort.
        - `first_child_delay` (:class:`int`) - the delay between the parent cohort and 
          the first child cohort. This delay is expressed in number of cohorts.

    :Returns:
        The indices of the child cohorts and their positions in the tree.
    
    :Returns Type:
        list of tuples
    
    '''
    child_cohorts = []
    if parent_cohort_index is None:
        first_possible_cohort_number = 1
    else:
        first_possible_cohort_number = parent_cohort_index + first_child_delay
    if first_possible_cohort_number == 1:
        # The main stem always exists: add it.
        cohort_position = 'MS'
        child_cohorts.append((first_possible_cohort_number, 'MS'))
        child_cohorts.extend(decide_child_cohorts(decide_child_cohort_probabilities, 
                                                  first_possible_cohort_number, 
                                                  cohort_position))
    else:
        # Find the children of the secondary stem.
        for cohort_number_str, cohort_probability in decide_child_cohort_probabilities.iteritems():
            cohort_number = int(cohort_number_str)
            if cohort_number >= first_possible_cohort_number:
                if cohort_probability >= random.random():
                    child_cohort_position = cohort_number - parent_cohort_index - first_child_delay
                    if parent_cohort_position == 'MS':
                        cohort_position = 'T%s' % child_cohort_position
                    else:
                        cohort_position = '%s.%s' % (parent_cohort_position, child_cohort_position)
                    child_cohorts.append((cohort_number, cohort_position)) 
                    child_cohorts.extend(decide_child_cohorts(decide_child_cohort_probabilities, 
                                                              cohort_number, 
                                                              cohort_position))
    return child_cohorts


def calculate_MS_final_leaves_number(MS_leaves_number_probabilities):
    '''
    Calculate the final number of leaves of a main stem. This is done by randomly 
    drawing a number in a probability distribution. Uses the probabilities 
    of the main stem leaves number. 
    
    :Parameters:
    
        - `MS_leaves_number_probabilities` (:class:`dict`) - the probabilities 
          of the main stem leaves number.
          
    :Returns:
        The final number of leaves of the main stem.
    
    :Returns Type:
        :class:`float`
        
    '''
    random_value = random.random()
    probabilities_sum = 0.0
    MS_final_leaves_number = None
    for leaves_number_str, leaves_probability in MS_leaves_number_probabilities.iteritems():
        probabilities_sum += leaves_probability
        if random_value <= probabilities_sum:
            MS_final_leaves_number = float(leaves_number_str)
            break
    return MS_final_leaves_number


def calculate_tiller_final_leaves_number(MS_final_leaves_number, cohort_number, secondary_stem_leaves_number_coefficients):
    '''
    Calculate the final number of leaves of a tiller.  
    Uses the final number of leaves of  the main stem, the index of the cohort to 
    which belongs the tiller, and specific coefficients.  
    
    :Parameters:
    
        - `MS_final_leaves_number` (:class:`float`) - the final number of leaves of the 
          main stem.
        - `cohort_number` (:class:`int`) - the index of cohort.
        - `secondary_stem_leaves_number_coefficients` (:class:`dict`) - The coefficients 
          a_1 and a_2 to calculate the final number of leaves on tillers from 
          the final number of leaves on main stem. Calculation is done as follow::
        
              tiller_final_leaves_number 
                  = a_1 * MS_final_leaves_number - a_2 * cohort_number

    :Returns:
        The final number of leaves of a tiller.
    
    :Returns Type:
        :class:`float`
        
    '''
    a_1 = secondary_stem_leaves_number_coefficients['a_1']
    a_2 = secondary_stem_leaves_number_coefficients['a_2']
    return a_1* MS_final_leaves_number - a_2 * cohort_number
    

def decide_time_of_death(max_axes_number, min_axes_number, TT_em_phytomer1, TT_bolting, TT_flag_leaf_ligulation):
    '''
    Decide the thermal times (relative to canopy emergence) when the axes stop 
    growing. Uses a linear function which describes the decay of the global population. 

    :Parameters:
    
        - `max_axes_number` (:class:`int`) - the maximum number of existing axes.
        - `min_axes_number` (:class:`int`) - the minimum number of existing axes. 
        - `TT_em_phytomer1` (:class:`list`) - Thermal times (relative to canopy emergence) 
          of tip emergence of the first true leaf (not coleoptile or prophyll)
        - `TT_bolting` (:class:`float`) - date in thermal time at which the bolting starts.
        - `TT_flag_leaf_ligulation` (:class:`float`) - the thermal time of the flag leaf ligulation.

    :Returns: 
        the thermal times (relative to canopy emergence) when the axes stops growing.
    
    :Returns Type: 
        :class:`list`
        
    .. warning:: 
    
        * *min_axes_number*, *max_axes_number*, *TT_bolting* and *TT_flag_leaf_ligulation* 
          must be positive or null.
        * *TT_bolting* must be smaller (or equal) than *TT_flag_leaf_ligulation*.
        * *min_axes_number* must be smaller (or equal) than *max_axes_number*.

    '''
    
    checkValidity(max_axes_number >= 0 and min_axes_number >=0 and TT_bolting >= 0 and TT_flag_leaf_ligulation >= 0)
    checkValidity(TT_bolting <= TT_flag_leaf_ligulation)
    checkValidity(min_axes_number <= max_axes_number)
    
    polynomial_coefficient_array = np.polyfit([TT_flag_leaf_ligulation, TT_bolting], [min_axes_number, max_axes_number], 1)
                
    remaining_axes_number = max_axes_number
    T_em_leaf1_tuples = zip(TT_em_phytomer1[:], range(len(TT_em_phytomer1)))
    T_em_leaf1_tuples.sort()
    T_stop_axis_tuples = []
    for tt in range(int(TT_bolting), int(TT_flag_leaf_ligulation) + 1):
        simulated_axes_number = int(np.polyval(polynomial_coefficient_array, tt))
        axes_to_delete_number = remaining_axes_number - simulated_axes_number
        while axes_to_delete_number > 0:
            max_emf_1, axis_row_number = T_em_leaf1_tuples.pop()
            T_stop_axis_tuples.append((axis_row_number, tt))
            axes_to_delete_number -= 1
            remaining_axes_number -= 1
        if remaining_axes_number == 0:
            break 
    T_stop_axis_tuples.sort()
    T_stop_axis_row_number_list = [T_stop_axis_tuple[0] for T_stop_axis_tuple in T_stop_axis_tuples]
    TT_stop_axis_list = [T_stop_axis_tuple[1] for T_stop_axis_tuple in T_stop_axis_tuples]
    for i in range(len(TT_em_phytomer1)):
        if i not in T_stop_axis_row_number_list:
            TT_stop_axis_list.insert(i, np.nan)
    return TT_stop_axis_list 


def fit_poly(x_meas_array, y_meas_array, fixed_coefs, a_starting_estimate):
    '''
    Calculate the best-fit parameter *a*, where *a* is the coefficient of highest 
    degree of a polynomial. The other polynomial coefficients are supposed to be 
    known and are given in *fixed_coefs*.
    We first define a function to compute the residuals. Then we use the least-squares 
    fit routine of scipy to find the best-fit parameter *a*, selecting *a_starting_estimate* 
    as starting position and using (*x_meas_array*, *y_meas_array*) as the measured 
    data to fit. Finally, we calculate the *RMSE* to check the validity of the fit. 

    .. seealso:: :func:`scipy.optimize.leastsq` 

    :Parameters:
    
        - `x_meas_array` (:class:`np.ndarray`) - the x-coordinates. These data are 
          measured.
        - `y_meas_array` (:class:`np.ndarray`) - the y-coordinates. These data are 
          measured.
        - `fixed_coefs` (:class:`list`) - the other coefficients of the polynomial to fit 
          (*x_meas_array*, *y_meas_array*) to. 
          These coefficients are not fitted. They are given from highest degree 
          to lowest degree ("descending powers").
        - `a_starting_estimate` (float) - the starting estimate for the minimization.
          
    :Returns: 
        the best-fit coefficient *a* and the *RMSE* of the fit.
    
    :Returns Type: 
        :class:`tuple` of :class:`float`
        
    '''
    def residuals(p, y, x):
        a, = p
        err = y - peval(x, a)
        return err
    def peval(x, a):
        return np.poly1d([a] + fixed_coefs)(x)
    p, cov, infodict, mesg, ier = leastsq(residuals, [a_starting_estimate], args=(y_meas_array, x_meas_array), full_output=1)
    # RMSE_gl
    chisq = (infodict['fvec']**2).sum()
    dof = len(x_meas_array) - 1 # dof is degrees of freedom
    rmse = np.sqrt(chisq / dof)
    return p[0], rmse


def calculate_theoretical_cohorts_cardinalities(plant_number, 
                                                decide_child_cohort_probabilities, 
                                                first_child_delay):
    '''
    Calculate the theoretical cardinality of each cohort. 
    
    :Parameters:
    
        - `plant_number` (:class:`int`) - the number of plants. 
        - `decide_child_cohort_probabilities` (:class:`dict`) - the probabilities of the child 
          cohorts. 
        - `first_child_delay` (:class:`int`) - The delay between 
          a parent cohort and its first possible child cohort. This delay is 
          expressed in number of cohorts.
          
    :Returns:
        a dictionary which contains, for each cohort, the cardinality of the cohort. 
        Keys are the index of the cohorts (int), values are the cardinalities of 
        the cohorts (float).
    
    :Returns Type:
        :class:`dict`
        
    '''
    decide_cohort_probabilities = decide_child_cohort_probabilities.copy()
    # the cohort '1' always exists, so its probability is 1.0.
    decide_cohort_probabilities['1'] = 1.0
    possible_cohorts = np.array(decide_cohort_probabilities.keys()).astype(int)
    decide_cohort_probability_values = np.array(decide_cohort_probabilities.values())
    possible_child_cohorts = np.array(decide_child_cohort_probabilities.keys()).astype(int)
    possible_parent_cohorts = possible_child_cohorts - first_child_delay
    commons = np.where(np.intersect1d(possible_cohorts, possible_parent_cohorts))
    possible_parent_cohorts = possible_cohorts[commons]
    decide_parent_cohort_probability_values = decide_cohort_probability_values[commons]
    
    cohort_cardinalities_dataframe = pandas.DataFrame(index=range(len(decide_cohort_probabilities)), 
                                         columns=['cohort', 
                                                  'theoretical_cardinality'])
    theoretical_probabilities_series = pandas.Series(index=cohort_cardinalities_dataframe.index)
    idx = 0
    for (cohort_str, decide_probability) in decide_cohort_probabilities.iteritems():
        cohort_int = int(cohort_str)
        cohort_cardinalities_dataframe['cohort'][idx] = cohort_int
        if cohort_str == '1':
            theoretical_probabilities_series[idx] = decide_cohort_probabilities['1']
        else:
            first_possible_parent = cohort_int - first_child_delay
            curr_possible_parent_indexes = np.where(possible_parent_cohorts <= first_possible_parent)
            curr_decide_parent_cohort_probabilities = decide_parent_cohort_probability_values[curr_possible_parent_indexes]
            theoretical_probabilities_series[idx] = (curr_decide_parent_cohort_probabilities * decide_probability).sum()
        idx += 1
    cohort_cardinalities_dataframe['theoretical_cardinality'] = theoretical_probabilities_series * plant_number
    
    return dict(zip(cohort_cardinalities_dataframe.values[:, 0], cohort_cardinalities_dataframe.values[:, 1],))
    

def checkValidity(is_valid):
    '''
    Raise an InputError exception when an invalid input is detected. 
    
    :Parameters:
    
        - `is_valid` (bool) - the result of the expression which has been used 
          to test the validity of an input.
    
    '''
    if not is_valid:
        raise InputError()


class Error(Exception):
    '''Base class for exceptions in :mod:`alinea.adel.plantgen.axeT`.'''
    pass


class InputError(Error):
    '''Exception raised when an invalid input is detected.'''
    def __init__(self):
        self.message = '''Invalid input detected ! Look at the traceback to find 
the invalid input.'''
        
    def __str__(self):
        return self.message


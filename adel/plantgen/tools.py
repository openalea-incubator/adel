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
import warnings

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
        for cohort_id, cohort_probability in decide_child_cohort_probabilities.iteritems():
            if cohort_id >= first_possible_cohort_number:
                if cohort_probability >= random.random():
                    child_cohort_position = cohort_id - parent_cohort_index - first_child_delay
                    if parent_cohort_position == 'MS':
                        cohort_position = 'T%s' % child_cohort_position
                    else:
                        cohort_position = '%s.%s' % (parent_cohort_position, child_cohort_position)
                    child_cohorts.append((cohort_id, cohort_position)) 
                    child_cohorts.extend(decide_child_cohorts(decide_child_cohort_probabilities, 
                                                              cohort_id, 
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
    
    if max_axes_number < 0:
        raise InputError("max_axes_number negative")
    if min_axes_number < 0:
        raise InputError("min_axes_number negative")
    if TT_bolting < 0:
        raise InputError("TT_bolting negative")
    if TT_flag_leaf_ligulation < 0:
        raise InputError("TT_flag_leaf_ligulation negative")
    
    if TT_bolting > TT_flag_leaf_ligulation:
        raise InputError("TT_bolting greater than TT_flag_leaf_ligulation")
    
    if min_axes_number > max_axes_number:
        raise InputError("min_axes_number greater than max_axes_number")
    
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


def calculate_theoretical_cardinalities(plant_number, 
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
        a dictionary which contains the theoretical cardinality of each cohort. 
        Keys are the index of the cohorts (str), values are the cardinalities of 
        the cohorts (float).
    
    :Returns Type:
        :class:`dict`
        
    '''
    decide_cohort_probabilities = decide_child_cohort_probabilities.copy()
    # the cohort 1 always exists, so its probability is 1.0.
    decide_cohort_probabilities[1] = 1.0
    possible_cohorts = np.array(decide_cohort_probabilities.keys())
    decide_cohort_probability_values = np.array(decide_cohort_probabilities.values())
    possible_child_cohorts = np.array(decide_child_cohort_probabilities.keys())
    possible_parent_cohorts = possible_child_cohorts - first_child_delay
    commons = np.where(np.intersect1d(possible_cohorts, possible_parent_cohorts))
    possible_parent_cohorts = possible_cohorts[commons]
    decide_parent_cohort_probability_values = decide_cohort_probability_values[commons]
    
    cohort_cardinalities = pandas.DataFrame(index=range(len(decide_cohort_probabilities)), 
                                            columns=['id_cohort', 'theoretical_cardinality'])
    theoretical_probabilities_series = pandas.Series(index=cohort_cardinalities.index)
    idx = 0
    for (cohort, decide_probability) in decide_cohort_probabilities.iteritems():
        cohort_cardinalities['id_cohort'][idx] = cohort
        if cohort == 1:
            theoretical_probabilities_series[idx] = decide_cohort_probabilities[1]
        else:
            first_possible_parent = cohort - first_child_delay
            curr_possible_parent_indexes = np.where(possible_parent_cohorts <= first_possible_parent)
            curr_decide_parent_cohort_probabilities = decide_parent_cohort_probability_values[curr_possible_parent_indexes]
            theoretical_probabilities_series[idx] = (curr_decide_parent_cohort_probabilities * decide_probability).sum()
        idx += 1
    cohort_cardinalities['theoretical_cardinality'] = theoretical_probabilities_series * plant_number
    
    return dict(zip(cohort_cardinalities.values[:, 0], cohort_cardinalities.values[:, 1],))
    
    
def calculate_decide_child_cohort_probabilities(decide_child_axis_probabilities):
    '''
    For each primary tiller in *decide_child_axis_probabilities*, calculate the corresponding 
    cohort number, and return a dictionary which keys are cohort number and values remain 
    the same.
    
    :Parameters:
    
        - `decide_child_axis_probabilities` (:class:`dict`) - the probability for each 
          primary tiller to have a child. Keys are the botanical positions (e.g. "T1", "T2",...), 
          values are the probabilities (float).
          
    :Returns:
        the probability for each cohort to have a child. Keys are the indexes of the cohorts 
        (e.g. 3, 4,...), values are the probabilities (float).
          
    :Returns Type:
        :class:`dict`
    
    :Examples:
    
        >>> decide_child_axis_probabilities = {'T0': 0.0, 'T1': 0.900, 'T2': 0.983, 'T3': 0.817, 'T4': 0.117}
        >>> calculate_decide_child_cohort_probabilities(decide_child_axis_probabilities)
        {3: 0.0, 4: 0.900, 5: 0.983, 6: 0.817, 7: 0.117}
        
    '''
    id_axis_array = np.array(decide_child_axis_probabilities.keys())
    id_cohort_array = np.char.lstrip(id_axis_array, 'T').astype(int) + 3
    decide_child_cohort_probabilities = dict(zip(id_cohort_array, 
                                                 decide_child_axis_probabilities.values()))
    return decide_child_cohort_probabilities


def get_primary_axis(id_axis, first_child_delay):
    '''
    Calculate the primary axis of *id_axis*.
    
    :Parameters:
    
        - `id_axis` (:class:`str`) - the botanical position of the axis.
        - `first_child_delay` (:class:`int`) - the delay between the axis and its 
          first child.
          
    :Returns:
        the primary axis of *id_axis*. 
          
    :Returns Type:
        :class:`str`
    
    :Examples:
    
        >>> get_primary_axis('T1.0', 2)
        'T3'
        >>> get_primary_axis('T1.0.0', 2)
        'T5'
        
    '''
    id_axis = id_axis[1:]
    while '.' in id_axis:
        id_axis_split = id_axis.rsplit('.', 2)
        last_pos = int(id_axis_split.pop())
        last_but_one_pos = int(id_axis_split.pop())
        new_last_pos = last_but_one_pos + last_pos + first_child_delay
        id_axis = '.'.join(id_axis_split + [str(new_last_pos)])
    primary_id_axis = 'T' + id_axis
    return primary_id_axis


def get_real_roots(poly):
    '''
    Get the real roots of polynomial *poly*.
    
    :Parameters:
    
        - `poly` (:class:`numpy.lib.polynomial.poly1d`) - a one-dimensional polynomial.
          
    :Returns:
        the real roots of *poly* 
          
    :Returns Type:
        :class:`numpy.array`
        
    :Examples:
    
        >>> p = numpy.poly1d([1, 2, 1])
        array([-1., -1.])
        >>> p = numpy.poly1d([1, 2, 3])
        array([], dtype=float64)
    
    '''
    roots_array = poly.r
    real_roots = filter(lambda x: x.imag == 0.0, roots_array)
    real_roots = map(lambda x: x.real, real_roots)
    return np.array(real_roots)
            

class InputError(Exception):
    '''Exception raised when an invalid input is detected.'''
    pass
    

class InputWarning(UserWarning):
    '''Warning issued when an input is dubious and may lead to an error.'''
    pass
    

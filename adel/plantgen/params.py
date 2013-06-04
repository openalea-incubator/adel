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
This module stores the constant parameters used in the :mod:`alinea.adel.plantgen` package. The user 
can set their value from here.

Authors: Mariem Abichou, Camille Chambon, Bruno Andrieu
'''

secondary_stem_leaves_number_coefficients = {'a_1': 0.9423, 'a_2': 0.555}
'''The coefficients *a_1* and *a_2* to calculate the final number of leaves on 
tillers from the final number of leaves on main stem. 

Calculation is done as follow::
    
    tiller_final_leaves_number 
        = a_1 * MS_final_leaves_number - a_2 * cohort_number
    
'''


emf_1_MS_standard_deviation = 30.0
'''the standard deviation in the thermal of emergence of plants in the plot.

This parameter is used to calculate main stem emf_1 value.'''


leaf_number_delay_MS_cohort = {3: 1.6173, 4: 2.5181, 5: 3.4189, 6: 4.5576, 7: 5.8097}
'''Delays between the emergence of the main stem and the emergence of each cohort.

The delays are expressed in main stem phyllochron unit.
One value is given for each cohort. It specifies the delay between a main 
stem having the most frequent number of leaves (for a main stem) and a tiller having the most 
frequent number of leaves (for that cohort).

This parameter is a Python dictionary. 
The keys represent the cohort indexes and the values represent the delays. 
The keys are integers ans the values are floats.
'''


n2_MS_div_n2_cohort = 0.85
'''Ratio between the maximum number of green leaves on the tillers and the 
maximum green leaves on the main stem.

Value is given for the axes with the most frequent leaves number.
'''


delais_phyll_col_tip = 1.6
'''Delay between tip emergence and collar emergence.

The delay is given in phyllochron unit and is taken the same for all leaf 
positions.
'''


delais_phyll_sen_disp = 3.0
'''The time during which a fully senesced leaf on a non-elongated internode 
remains on the plant. 

The delay is given in phyllochron unit. 
'''


TT_del_Fhaut = 3000
'''The thermal time at which leaves on elongated internode disappear.

The thermal time is given in degree.day. 
'''


first_child_delay = 2
'''The delay between a parent cohort and its first possible child cohort. 
This delay is expressed in number of cohorts.'''


regression_of_dimensions = {'L_blade': (0.06, 0.1, 0.15), 
                            'L_internode': (0.06, 0.1, 0.15), 
                            'L_sheath': (0.06, 0.1, 0.15), 
                            'W_blade': (0.06, 0.1, 0.15), 
                            'W_internode': (0.06, 0.1, 0.15), 
                            'W_sheath': (0.06, 0.1, 0.15)}
'''The regression of the dimensions for the last 3 phytomers of each organ.

Each key "dimension_organ" is associated to a 3-tuple. In each 3-tuple, the first 
value is the reduction of the phytomer N-2, the second value is the reduction of 
the phytomer N-1, and the third value is the reduction of the phytomer N, where 
N is the final number of leaves of an axis. The reduction are given in percent, 
from 0.0 to 1.0. 
'''

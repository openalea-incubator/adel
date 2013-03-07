
   
++++++++++++++++++++

.. sidebar:: Summary

    :Topic: *Alinea.Adel Documentation*
    :Release: |release|
    :Date: |today|
    :Authors: Chrisitian FOURNIER, Mariem ABICHOU, Christophe PRADAL, Bruno ANDRIEU and Camille CHAMBON
    :Target: users, developers and administrators
 
.. contents:: **Alinea.Adel Documentation**
   

General introduction (write authors names)
===========================================

ADEL-Wheat (Architectural model of DEvelopment based on L-systems) is designed for
simulating the 3D architectural development of the aerial part of wheat plants. The model has been
coupled with a light model adapted to field conditions (Caribu model), in order to be able to calculate accurately the
light interception by phytoelements during the crop cycle. Applications of such a tool range from the
interpretation of remote-sensing signals, to the estimation of crop light use efficiency, or to the
assistance in the ecophysiological analysis of plant response to light conditions. The model is based on
an analysis of developmental and geometrical similarities that exists among phytomers, to allow a
concise parameterisation. Finally, parameterising the model requires only: (i) the calibration of simple
functions that describe the mature size and geometry of successive organs along the main stem; (ii) the
timing of collar emergence on all axes; and (iii) a (estimated) developmental delay between main stem
and every reproductive tiller. Simulations were realistic (Fig 1). The model
correctly reproduced important agronomic features such as kinetics of LAI and of plant height.



.. _adel_input:

Description of Adel's inputs (write authors names)
==========================================================

Adel has three kinds of inputs available for the user:

 * inputs charactering the development of plants in the canopy (topology, developmental rate and size of plants)
 * inputs characterising the geometry of axes and of leaves
 * inputs charactering the simulation (time step..) and the plot configuration (number of plant, position)

Beside these user inputs, Adel also has 'fixed' parameters and methodologies (set and 'readable' in Adel.R, but not documented) describing coordination of leaves, dynamics of geometry, computation of visibility and progression of senescence as a function of ssi. See Robert et al (Functional Plant Biology 2008) for description.

.. image:: image/exemple1.png

.. _development_input:

Inputs describing development
********************************

Development is parameterised in a Rlist containing 3 tables (axeT,dimT and phenT) . 

These tables have dependencies (cross references). However some may be compatible with others if cross references are. This allows for recombination of parameters

.. _axeT:

axeT : *description of  plant variability within the canopy*
------------------------------------------------------------    

This table characterises the plant to plant variability within the canopy. If only one plant is given, then adel will clone that plant. For wheat, to have a correct simulation of tiller dynamics, a minimum of 30 plants is recommended.

The table contents explicit parameters (eg number of leaves on axes) and indexes that reference a particular developmental pattern (parameter set) described in the other tables (dimT and phenT) to allows factorisation.

There is one line per axe (see axeTdemo.csv). Columns are :

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **plant**
      - the plant number the axe belongs to
    * - **axe**
      - axe number (starting at 0 for main stem) within the plant
    * - **nf**
      - final number of leaves on the axe
    * - **end**
      - *NA* or date of end of growth of the axe for regressing ones
    * - **disp**
      - *NA* or date of disparition the axe
    * - **dimIndex**
      - reference (int) to a parameter set in dimT
    * - **phenIndex**
      - reference to a parameter set in phenT
    * - **earIndex**
      - reference to a parameter set in earT 
    * - **emf1**
      - date of tip emergence
    * - **ligf1**       
      - collar emergence                              
    * - **senf1**
      - end of senescence of the first leaf of the axe
    * - **dispf1**
      - disparition of the first leaf of the axe     

.. Note :: This allows for parameterising heterogeneous emergence date.


.. _dimT:

dimT : *description of dimensions of axes (or a subset of axes)*
----------------------------------------------------------------


Dimensions could be given for all axes (using one dimIndex per axe) or a subset of axes. Dimension are given for relative position within the axe (=position/final number of leaves). Actual dimension of axe are hence dependent on dimIndex and nf.

.. Note :: add scale factor (or max dim) in axeT and describe here double normalised functions)

There is one line per phytomer (see dimTdemo.csv) of all indexed axes. Columns are :

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **index**
      - the index refered to in axeT
    * - **nrel** 
      - normalised phytomer position, starting from the base)
    * - **Ll**
      - blade length
    * - **Lw**
      - blade max width
    * - **Gl** 
      - sheath length
    * - **Gd** 
      - sheath diameter
    * - **El** 
      - internode length
    * - **Ed** 
      - internode diameter


.. _phenT:

phenT : *description of phenology of axes (or a subset of axes)*
-----------------------------------------------------------------


Phenology controls the rate of plant development (hence extension rates of organs), the dynamics of leaf appearance and the dynamics of senescence. 

Positions are normalised to final leaf numbers to allows sharing of data between axes of axeT table.

Dates of developmental events are given relative taking as origin the date of the event on leaf 1 of the axe. Actual development is computed from this table and the date concerning leaf 1 in axeT. 

There is one line per phytomer (see phenTdemo.csv) of all indexed axes. Columns are :

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **index** 
      - the index referred to in axeT
    * - **nrel** 
      - normalised phytomer position, starting from 0 (to allow extrapolation)
    * - **tip** 
      - date (origin tip leaf 1) of tip emergence of the phytomer
    * - **col** 
      - date (origin col leaf 1) of collar emergence of the phytomer
    * - **ssi** 
      - date (origin sen leaf 1) of full senescence of the phytomer (ssi)
    * - **disp** 
      - date (origin disp leaf 1) of leaf disappearance.Blade disappear at disp. Sheath disappear when leaf above it disappear

earT : *description of ear dimension and phenology*
----------------------------------------------------


There is one line per ear type (refered by ear Index in axeT)

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **index** 
      - the index refered to in axeT
    * - **em_ear** 
      - delay between flag leaf ligulation and ear (tip of highest spike without awn) appearance
    * - **em_ped** 
      - delay between flag leaf ligulation and peduncle (tip = base of th ear) appearance
    * - **end_gf** 
      - delay between flag leaf ligulation and end of grain filling (full senescence of the ear+stem)
    * - **l_ped** 
      - length of the peduncle
    * - **d_ped** 
      - diameter of the peduncle
    * - **l_ear** 
      - length of the ear (without awns)
    * - **Sp_ear** 
      - projected area of ear without awn
    * - **l_ear_awn** 
      - length of the ear+awns



ssi2sen : *desciption of progression of upper leaf senescence as a function of ssi*
------------------------------------------------------------------------------------

By default, leaves start senescence 1 ssi unit before ssi = leaf number and complete senecence when ssi = leaf number.
For the *ndel* upper leaves, senecence start at :math:`sssi = final\ leaf\ number - ndel(t0)` at a slower rate specified in this file (*rate*), 
and accelerate *dssit1* ssi unit after :math:`t0`, and leaves are fully senesced *dssit2* ssi unit after :math:`t0`.

The table allows for definition of *rate*, *dssit1* and *dssit2* for the ndel upper leaves.

*ndel* is given by the number of lines of the file.


Inputs describing geometry
*****************************

Geometry of leaves is defined by two lists of lists of matrices describing midrib curvature and leaf width variation with distance to the base of the leaf:

    * the first level in the list is for collection index
    * the second level is for matrix index. see alea


Beside leaf shapes two lists of R function should be provided as inputs.

The first list should provide 3 R function of axe number (0 = main stem) that returns:
    * **azT** : the azimuth(deg) of the first leaf of the axe with reference to the azimuth of the parent leaf
    * **incT** : the inclination (deg) of the base of the tiller compared with main stem
    * **dredT** : the distance (at maturity) between tiller and main stem

These functions could be created with the *genGeoAxe* node (with constraints) or freely defined with *freeGeoAxe*. 
A sample definition may be :

.. code-block:: r

	geoAxe <- list(
	  azT = function(a) {
	    ifelse(a == 0, 0, 75 + (runif(1)-0.5)*5) 
	  },
	  incT = function(a) {
		ifelse(a == 0, runif(1) * 5, 82 + (runif(1) - .5) * 5)
	  },
	  dredT = function(a) {
		ifelse(a == 0, 0, runif(1) * 7)
	  }
	)


The second list should provide two Rfunctions of axe number, 
leaf position and leaf position counted from top 
(plus leaf stage for Lindex, defined as curent length/final length). 
Returned values should be :

    * **azim** : the azimuth (deg) of the leaf compared to the previous one
    * **Lindex** : the index of the collection to use for leaf shape

This list could be generated by genGeoLeaf or freeGeoLeaf. 
A sample dfinition may be : 

.. code-block:: r

	geoLeaf <- list(
		Azim = function(a,n,ntop) {0 * runif(1)},
		Lindex = function(a,n,ntop,stage) {ntop + 1}
		)

Inputs describing simulation
********************************

Time step is given as a list of date for which a mock-up is wished
position of plants within the plot are given externally from adel to a planter.

Description of Adel's outputs (write authors names)
============================================================


.. _plantgen:

Construction of the input tables (Mariem ABICHOU, Bruno ANDRIEU and Camille CHAMBON)
============================================================================================

ADEL expects inputs characterising the development of plants in the canopy. These 
inputs are described in :ref:`development_input`.
ADEL user who does not have a complete set of inputs may wish to use ADEL anyway,
constructing the missing inputs. That's the aim of the ``plantgen`` package 
(see :mod:`alinea.adel.plantgen`). 

According to the level of completeness of the raw inputs, and given some parameters, 
the ``plantgen.plantgen`` module provides routines to construct 
:ref:`axeT <axeT>`, :ref:`dimT <dimT>` and :ref:`phenT <phenT>`, and some other 
dataframes for debugging purpose (see :mod:`alinea.adel.plantgen.plantgen`).

In the next subsections, we first describe the different levels of completeness 
of the inputs and of the parameters set by the user. 

Then we see how to construct the inputs of ADEL from a Python interpreter, using 
the routine ``gen_adel_input_data(...)``. This routine can be used whatever the 
level of completeness of the raw inputs, adapting the processing automatically 
(see :func:`alinea.adel.plantgen.plantgen.gen_adel_input_data`).

Finally, we see how to construct the inputs of ADEL from the Visualea interface, 
using the following convenience routines: 

    * ``gen_adel_input_data_from_min(...)``,
    * ``gen_adel_input_data_from_short(...)``,  
    * and ``gen_adel_input_data_from_full(...)``.
    
All these routines belong to :mod:`alinea.adel.plantgen.plantgen`. They 
permit to ease the construction from respectively a minimum, short, and full set 
of inputs. 


.. _levels_of_completeness:

The levels of completeness
***************************

The inputs and the parameters needed for the construction are ``dynT_user`` and 
``dimT_user``.
``dynT_user`` and ``dimT_user`` can have different levels of completeness: ``FULL``, 
``SHORT`` and ``MIN``. 
According to their level of completeness, ``dynT_user`` and ``dimT_user`` take 
different types, shapes and/or contents:

.. list-table::
    :widths: 10 25 25
    :header-rows: 1

    * - Level of completeness
      - dynT_user
      - dimT_user
    * - **FULL** 
      - :ref:`dynT_user_FULL`
      - :ref:`dimT_user_FULL`
    * - **SHORT** 
      - :ref:`dynT_user_SHORT`
      - :ref:`dimT_user_SHORT`
    * - **MIN** 
      - :ref:`dynT_user_MIN`
      - :ref:`dimT_user_MIN`
      
.. seealso:: :class:`alinea.adel.plantgen.plantgen.DataCompleteness`

      
.. _construct_inputs_from_interpreter:

Construct the inputs from Python interpreter
***********************************************

``gen_adel_input_data(...)`` is aimed to be used from Python interpreter 
(see :func:`alinea.adel.plantgen.plantgen.gen_adel_input_data`).

In the next subsections, we explain how to define the arguments of ``gen_adel_input_data(...)``.

.. note : in the examples below, the csv tables are supposed to be located in the 
          working directory.
          

* dynT_user : *the leaf dynamic parameters set by the user*

  *dynT_user* can be either a :class:`pandas.DataFrame` or a :class:`dict`, 
  depending on its level of completeness. See :ref:`levels_of_completeness` for more 
  details.

  For example, if *dynT_user_completeness* is ``SHORT``, then the user may import the 
  table: :ref:`dynT_user_SHORT_example <dynT_user_SHORT_example>` example, using :mod:`pandas` as follows::

      import pandas

      dynT_user = pandas.read_csv('dynT_user_SHORT.csv')
    

* dimT_user : *the dimensions of the axes set by the user*

  *dimT_user* is a :class:`pandas.DataFrame`, which content depends on 
  *dynT_user_completeness*. See :ref:`levels_of_completeness` for more details. 
    
  For example, if *dimT_user_completeness* is ``SHORT``, then the user may import the 
  table: :ref:`dimT_user_SHORT_example <dimT_user_SHORT_example>` example, using :mod:`pandas` as follows::
    
      import pandas
    
      dimT_user = pandas.read_csv('dimT_user_SHORT.csv')


* dynT_user_completeness and dimT_user_completeness : *the levels of completeness of dynT_user and dimT_user*

  *dynT_user_completeness* and *dimT_user_completeness* are the levels of completeness 
  of respectively *dynT_user* and *dimT_user* (see :ref:`levels_of_completeness`). 
    
  *dynT_user_completeness* and *dimT_user_completeness* have to be coherent with 
  respectively *dynT_user* and *dimT_user*.
    
  For example, if the levels of completeness of *dynT_user* and *dimT_user* are 
  both ``SHORT``, then *dynT_user_completeness* and *dimT_user_completeness* 
  must be defined as follows::
    
      from alinea.adel.plantgen.plantgen import DataCompleteness
    
      dynT_user_completeness = DataCompleteness.SHORT
      dimT_user_completeness = DataCompleteness.SHORT
    

* *plant_number*, *decide_child_cohort_probabilities*, *MS_leaves_number_probability_distribution*, ...

  The other arguments of the routine are: 
    
      * *plant_number*, the number of plants to be generated,
      * *decide_child_cohort_probabilities*, for each child cohort the probability of emergence of an axis when the parent axis is present,  
      * *MS_leaves_number_probability_distribution*, the probability distribution of the final number of main stem leaves,
      * *TT_bolting*, the date in thermal time at which the bolting starts,
      * *TT_flowering*, the flowering date in thermal time,
      * *final_axes_number*, the final number of axes which have an ear, per square meter,
      * *GL_number*, the GL decimal numbers measured at several thermal times (including the senescence end),
      * *delais_TT_stop_del_axis*, the thermal time between an axis stop growing and its disappearance,
      * and *TT_col_break*, the thermal time when the rate of Haun Stage is changing.
            
  They can be defined as follows::
  
      plant_number = 100
      decide_child_cohort_probabilities = {'3': 0.0, '4': 0.900, 
                              '5': 0.983, '6': 0.817, 
                              '7': 0.117}
      MS_leaves_number_probabilities = {'10': 0.145, 
                                        '11': 0.818, 
                                        '12': 0.036, 
                                        '13': 0.0, 
                                        '14': 0.0}
      TT_bolting = 500
      TT_flowering = 1440
      final_axes_number = 250
      GL_number = {1117.0: 5.6, 1212.1:5.4, 
                   1368.7:4.9, 1686.8:2.4, 
                   1880.0:0.0}
      delais_TT_stop_del_axis = 600
      TT_col_break = 0.0
        
    
  See :func:`alinea.adel.plantgen.plantgen.gen_adel_input_data` for more details.


Launch the construction 
--------------------------

To launch the construction, simply call ``gen_adel_input_data(...)`` 
with the appropriate arguments:: 

    from alinea.adel.plantgen.plantgen import gen_adel_input_data

    (axeT, 
    dimT, 
    phenT, 
    phenT_abs, 
    dimT_abs, 
    dynT, 
    phenT_first,
    HS_GL_SSI_T,
    tilleringT,
    cohortT) = gen_adel_input_data(dynT_user, 
                                   dimT_user, 
                                   plant_number, 
                                   decide_child_cohort_probabilities, 
                                   MS_leaves_number_probability_distribution, 
                                   TT_bolting, 
                                   TT_flowering, 
                                   final_axes_number, 
                                   GL_number, 
                                   delais_TT_stop_del_axis, 
                                   TT_col_break, 
                                   dynT_user_completeness, 
                                   dimT_user_completeness)

The returned values are all :class:`pandas.DataFrame`. 

*axeT*, *dimT* and *phenT* can be converted to csv files and used as ADEL inputs::

    # write axeT, dimT and phenT to csv files in the working directory, replacing
    # missing values by 'NA' and ignoring the indexes (the indexes are the labels of
    # the rows)
    axeT.to_csv('axeT.csv', na_rep='NA', index=False)
    dimT.to_csv('dimT.csv', na_rep='NA', index=False)
    phenT.to_csv('phenT.csv', na_rep='NA', index=False)

See :ref:`axeT <axeT>`, :ref:`dimT <dimT>`, :ref:`phenT <phenT>`, :ref:`phenT_abs <phenT_abs>`, 
:ref:`dimT_abs <dimT_abs>`, :ref:`dynT <dynT>`, :ref:`phenT_first <phenT_first>`, 
:ref:`HS_GL_SSI_T <HS_GL_SSI_T>`, :ref:`tilleringT <tilleringT>`, 
:ref:`cohortT <cohortT>`.


.. _construct_inputs_from_visualea:

Construct the data from Visualea
***********************************

The following routines are convenience routines to construct the inputs of ADEL: 

    * ``gen_adel_input_data_from_min(...)``: construct the inputs of ADEL from 
      :ref:`dynT_user_MIN` and :ref:`dimT_user_MIN`,
    * ``gen_adel_input_data_from_short(...)``: construct the inputs of ADEL from 
      :ref:`dynT_user_SHORT` and :ref:`dimT_user_SHORT`,  
    * and ``gen_adel_input_data_from_full(...)``: construct the inputs of ADEL from 
      :ref:`dynT_user_FULL` and :ref:`dimT_user_FULL`.
    
All these routines belong to :mod:`alinea.adel.plantgen.plantgen`.

These routines are wrapped in the following Visualea nodes:

.. list-table::
    :widths: 10 10 10
    :header-rows: 1

    * - ``plantgen_MIN``
      - ``plantgen_SHORT``
      - ``plantgen_FULL``
    * - .. image:: image/plantgen_MIN.png
      - .. image:: image/plantgen_SHORT.png
      - .. image:: image/plantgen_FULL.png

The following table summarizes the nodes, the routines and the levels of completeness 
of :ref:`dynT <dynT>` and :ref:`dimT <dimT>`:

.. list-table::
    :widths: 15 15 40 30
    :header-rows: 1

    * - Completeness of :ref:`dynT`
      - Completeness of :ref:`dimT <dimT>`
      - Convenience routine
      - Visualea node
    * - **MIN** 
      - **MIN**
      - ``gen_adel_input_data_from_min``
      - ``plantgen_MIN``
    * - **SHORT** 
      - **SHORT**
      - ``gen_adel_input_data_from_short``
      - ``plantgen_SHORT``
    * - **FULL** 
      - **FULL**
      - ``gen_adel_input_data_from_full``
      - ``plantgen_FULL``
 
The following dataflow demonstrates how to use ``plantgen_MIN``, ``plantgen_SHORT``, 
and ``plantgen_FULL`` through Visualea:

.. image:: image/plantgen_dataflow.png

This dataflow is accessible from the Package explorer of Visualea, in 
``alinea.adel.tutorials.plangen``.


Appendices
***********

The appendices contain the description of the tables or dictionaries referred to 
in :ref:`construct_inputs_from_interpreter` and :ref:`construct_inputs_from_visualea`.

.. _dynT_user_FULL:

dynT_user_FULL
---------------

*dynT_user_FULL* is a table which contains the dynamic of the leaves. Actually, 
each type of axis is described by one row. The type of an axis is defined by its 
cohort index and its final number of leaves. Each row contains the following data: 
*N_cohort*, *Nff*, *a_cohort*, *TT_col_0*, *TT_col_nff*, *n0*, *n1* and *n2*.
See :ref:`dynT` for a description of these data.

.. _dynT_user_FULL_example:

Example:

    .. csv-table::
        :file: ./data/dynT_user_FULL.csv
        :header-rows: 1

.. seealso:: :download:`dynT_user_FULL.csv <./data/dynT_user_FULL.csv>`


.. _dynT_user_SHORT:

dynT_user_SHORT
----------------

*dynT_user_SHORT* is a table which contains the dynamic of a subset of the leaves. 
Actually, there is one row for each cohort, each row referring to the most frequent 
axis of the cohort. Each row contains the following data: 
*N_cohort*, *a_cohort*, *TT_col_0*, *TT_col_nff*, *n0*, *n1* and *n2*.
See :ref:`dynT` for a description of these data.

.. _dynT_user_SHORT_example:

Example:

    .. csv-table::
        :file: ./data/dynT_user_SHORT.csv
        :header-rows: 1
        
.. seealso:: :download:`dynT_user_SHORT.csv <./data/dynT_user_SHORT.csv>`


.. _dynT_user_MIN:

dynT_user_MIN
--------------

*dynT_user_MIN* is a dictionary which contains the dynamic of a subset of the 
leaves. This subset is composed by the leaves of the main stem. 
The dictionary contains the following keys: *a_cohort*, *TT_col_0*, *TT_col_nff*, 
*n0*, *n1* and *n2*. 
See :ref:`dynT` for a description of these data.

.. _dynT_user_MIN_example:

Example::

    # first, define TT_col_nff: the thermal time 
    # when Haun Stage is equal to Nff
    TT_col_nff = {'1': 1078, '4': 1148, '5': 1158, 
                  '6': 1168, '7': 1178}
    # then define dynT_user_MIN, which includes TT_col_nff
    dynT_user_MIN = {'a_cohort': 0.0102, 
                     'TT_col_0': -0.771289027, 
                     'TT_col_nff': TT_col_nff, 
                     'n0': 4.871559739, 
                     'n1': 3.24283148, 
                     'n2': 5.8}


.. _dimT_user_FULL:

dimT_user_FULL
----------------

*dimT_user_FULL* is the same as :ref:`dimT <dimT>`. 

.. _dimT_user_FULL_example:

Example:

    .. csv-table::
        :file: ./data/dimT_user_FULL.csv
        :header-rows: 1
       
.. seealso:: :download:`dimT_user_FULL.csv <./data/dimT_user_FULL.csv>`


.. _dimT_user_SHORT:

dimT_user_SHORT
----------------

*dimT_user_SHORT* is a table which contains the dimensions of the organs, for each 
phytomer of the most frequent axis of each cohort. Each row contains the following 
data: *id_axis*, *index_phytomer*, *L_blade*, *W_blade*, *L_sheath*, *W_sheath*, 
*L_internode* and *W_internode*. 
*id_axis* is the index of the cohort to which belongs the current most frequent axis.
See :ref:`dimT <dimT>` for a description of the other data. 

.. _dimT_user_SHORT_example:

Example:

    .. csv-table::
        :file: ./data/dimT_user_SHORT.csv
        :header-rows: 1
        
.. seealso:: :download:`dimT_user_SHORT.csv <./data/dimT_user_SHORT.csv>`


.. _dimT_user_MIN:

dimT_user_MIN
--------------

*dimT_user_MIN* is a table which contains the dimensions of the organs, for each 
phytomer of the most frequent axis of the main stem. Each row contains the following 
data: *index_phytomer*, *L_blade*, *W_blade*, *L_sheath*, *W_sheath*, *L_internode* 
and *W_internode*.
See :ref:`dimT <dimT>` for a description of these data. 

.. _dimT_user_MIN_example:

Example:

    .. csv-table::
        :file: ./data/dimT_user_MIN.csv
        :header-rows: 1

.. seealso:: :download:`dimT_user_MIN.csv <./data/dimT_user_MIN.csv>`


.. _phenT_abs:

phenT_abs
----------

:ref:`phenT_abs` is exactly the same as :ref:`phenT <phenT>`, except that:
    * the positions of the phytomers are not normalized,
    * the dates of developmental events are absolute.
    
:ref:`phenT_abs` is an intermediate dataframe used to construct :ref:`phenT <phenT>`.

.. _phenT_abs_example:

Example:

    .. csv-table::
        :file: ./data/phenT_abs.csv
        :header-rows: 1

.. seealso:: :download:`phenT_abs.csv <./data/phenT_abs.csv>`
      

.. _dimT_abs:

dimT_abs
----------

:ref:`dimT_abs` is exactly the same as :ref:`dimT <dimT>`, except that the positions 
of the phytomers are not normalized.

:ref:`dimT_abs` is an intermediate dataframe used to construct :ref:`dimT <dimT>`.  

.. _dimT_abs_example:

Example:

    .. csv-table::
        :file: ./data/dimT_abs.csv
        :header-rows: 1

.. seealso:: :download:`dimT_abs.csv <./data/dimT_abs.csv>`


.. _dynT:        

dynT
-----

:ref:`dynT` is a table which contains the dynamic of the leaves, for each type 
of axis. The type of an axis is defined by its cohort index and its final number 
of leaves. 
There is one row per type of axis. Each row contains the following data:

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **N_cohort** 
      - the index of the cohort to which belongs **id_axis**
    * - **id_axis** 
      - the current type of axis. This type is the concatenation of **N_cohort** 
        and **Nff**.
    * - **cardinality**
      - the cardinality of the set composed of **id_axis**
    * - **Nff** 
      - the final number of leaves of **id_axis**
    * - **a_cohort** 
      - the rate of Haun Stage vs Thermal time. This is the rate of the 
        first phase in case of bilinear behavior.
    * - **TT_col_0** 
      - the thermal time for Haun Stage equal to 0
    * - **TT_col_break**
      - the thermal time when the rate of phytomers emergence is changing
    * - **TT_col_nff** 
      - the thermal time when Haun Stage is equal to **Nff**
    * - **n0** 
      - number of green leaves at **t0**
    * - **n1** 
      - number of green leaves at **t1**
    * - **n2** 
      - number of green leaves at **TT_col_nff**
    * - **t0**
      - the thermal time at the start of leaf senescence 
    * - **t1**
      - the date in thermal time at which the senescence starts
    * - **hs_t1**
      - the Haun Stage at t1
    * - **a**
      - the coefficient of the 3rd order term of the polynomial describing the 
        dynamics of Green Leaf number after flowering 
    * - **c**
      - the coefficient of the 1st order term of the polynomial describing the 
        dynamics of Green Leaf number after flowering 
    * - **RMSE_gl**
      - the RMSE for the dynamic of green leaf number after estimation of 
        parameter a.

The rows are ordered by cohort index (**N_cohort**), and, within each cohort index, 
by **cardinality**.   


.. _dynT_example:

Example:

    .. csv-table::
        :file: ./data/dynT.csv
        :header-rows: 1

.. seealso:: :download:`dynT.csv <./data/dynT.csv>`
        

.. _phenT_first:

phenT_first
------------

:ref:`phenT_first` is a subset of :ref:`phenT_abs`. Actually, :ref:`phenT_first` 
contains only the rows of :ref:`phenT_abs` which correspond to the first phytomer 
of each axis. These rows have *index_phytomer* equal to 1. 

:ref:`phenT_first` is an intermediate dataframe used to construct :ref:`phenT <phenT>`.

.. _phenT_first_example:

Example:

    .. csv-table::
        :file: ./data/phenT.csv
        :header-rows: 1

.. seealso:: :download:`phenT.csv <./data/phenT.csv>`


.. _HS_GL_SSI_T:

HS_GL_SSI_T
------------

:ref:`HS_GL_SSI_T` describes, for each type of axis, the dynamic of *HS*, *GL* 
and *SSI* when *TT* varies. The type of an axis is defined by its cohort index 
and its final number of leaves.

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **id_axis** 
      - the current type of axis. The type of an axis is defined by its cohort 
        index and its final number of leaves.
    * - **TT** 
      - the thermal time.
    * - **HS** 
      - the Haun Stage.
    * - **GL** 
      - the number of green leaves.
    * - **SSI** 
      - the number of senescent leaves.
      
.. note::

   For each axis, *TT* varies from 0 to :attr:`alinea.adel.plantgen.params.TT_del_Fhaut`. 
   
:ref:`HS_GL_SSI_T` is constructed for debugging purpose.        

.. _HS_GL_SSI_T_example:

Example:

    .. csv-table::
        :file: ./data/HS_GL_SSI_T.csv
        :header-rows: 1
        
.. seealso:: :download:`HS_GL_SSI_T.csv <./data/HS_GL_SSI_T.csv>`

.. note:: this is a shortened version of :ref:`HS_GL_SSI_T`.


.. _tilleringT:

tilleringT
------------

:ref:`tilleringT` describes the dynamic of tillering. It stores the number of axes at 
important dates: the start of growth, the bolting date, and the flowering date.

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **TT** 
      - the date in thermal time.
    * - **NbrAxes** 
      - the number of axes.

:ref:`tilleringT` is constructed for debugging purpose.

.. _tilleringT_example:

Example:

    .. csv-table::
        :file: ./data/tilleringT.csv
        :header-rows: 1

.. seealso:: :download:`tilleringT.csv <./data/tilleringT.csv>`


.. _cohortT:

cohortT
------------

:ref:`cohortT` describes the theoretical and the simulated cardinalities of 
each cohort. It permits the user to validate the simulated cardinalities against 
the theoretical ones. The theoretical cardinalities are calculated from the 
probabilities of emergence of an axis when the parent axis is present (given by 
the user). The simulated cardinalities are calculated for each plant using 
:func:`alinea.adel.plantgen.tools.decide_child_cohorts`.

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **cohort** 
      - the index of the cohort
    * - **theoretical_cardinality** 
      - the theoretical cardinality
    * - **simulated_cardinality** 
      - the simulated cardinality

:ref:`cohortT` is constructed for debugging purpose.

.. _cohortT_example:

Example:

    .. csv-table::
        :file: ./data/cohortT.csv
        :header-rows: 1

.. seealso:: :download:`cohortT.csv <./data/cohortT.csv>`

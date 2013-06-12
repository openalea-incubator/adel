
   
++++++++++++++++++++

.. sidebar:: Summary

    :Topic: *Alinea.Adel Documentation*
    :Release: |release|
    :Date: |today|
    :Authors: Christian FOURNIER, Mariem ABICHOU, Christophe PRADAL, Bruno ANDRIEU and Camille CHAMBON
    :Target: users, developers and administrators
 
.. contents:: **Alinea.Adel Documentation**
   

General introduction 
=====================

ADEL-Wheat (Architectural model of DEvelopment based on L-systems) is designed for
simulating the 3D architectural development of the shoots of wheat plants. The model has been
coupled with a light model adapted to field conditions (Caribu model), in order to allow calculating the
light interception by individual phytoelements during the crop cycle. Applications of such a tool range from the
interpretation of remote-sensing signals, to the estimation of crop light use efficiency, or to the
assistance in the ecophysiological analysis of plant response to light conditions.
The model is based on an analysis of developmental and geometrical similarities that exists among phytomers, to allow a
concise parameterisation. Parameterising the model requires using experimental data to document a set of inputs described below.

Beside these user inputs, Adel also make use of constants and of relations (set and 'readable' in Adel.R, but not documented)
describing coordination of leaves, dynamics of geometry, computation of visibility and progression of senescence as a function of ssi.
See Robert et al (Functional Plant Biology 2008) for description.


Simulations are visually realistic (Fig 1). And the model was shown to correctly reproduce important agronomic features such as kinetics of LAI and of plant height.


.. _adel_input:

Description of Adel's inputs
=============================

Adel requires three kinds of inputs to be provided by the user:

 * inputs characterizing the development of plants in the canopy (topology, developmental rate and size of plants)
 * inputs characterizing the geometry of axes and of leaves
 * inputs characterizing the simulation (time step..) and the plot configuration (number of plants, position)


Units and conventions
- Dimension are expressed in cm
- Thermal time is expressed in °CD

Position of a phytomer on an axis: Most often (as described below) phytomer position in Adel are counted acropetally and are normalized relatively to the total number of phytomers
(leaf 1 relative position is 1/nf and flag leaf relative position is 1). Use of relative positions was chosen to allow sharing data between axes differing by the total number of phytomers
In later versions of Adel the relative phytomer number should be changed by the absolute one. With the convention that 1 is referring to the first true leaf.

.. figure:: image/global_dataflow.png
   :width: 100%

   Global data-flow of Adel.


.. _development_input:

Inputs describing development
********************************

Development is parameterized in a Rlist containing 3 tables (:ref:`axeT <axeT>`, :ref:`dimT <dimT>` and :ref:`phenT <phenT>`). 

These tables have dependencies (cross references). However some may be compatible with others if cross references are ???. This allows for recombination of parameters.

.. _axeT:

axeT : *Master table allowing to organize the information plant per plant*
---------------------------------------------------------------------------

This table is the master table that organizes how each plant is described.
For each plant, the table contains a few explicit parameters that describe the phenology and the number of modules (eg time of emergence, number of axes and number of leaves on axes)
and identifiers that refer to information given in the other tables (:ref:`dimT <dimT>`, :ref:`phenT <phenT>`, :ref:`earT <earT>`).

All plants to be used for the reconstruction must be listed in Axis_Table. If only one plant is given, Adel will clone that plant. 
To have a correct simulation of tiller dynamics at the plot level, a minimum of 30 plants is recommended.

There is one line per axis. Columns are :

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **id_plt**
      - Number (int) identifying the plant to which the axe belongs
    * - **id_cohort_axis**
      - Number (int) identifying the cohort to which the axe belongs
    * - **id_axis**
      - Identifier of the botanical position of the axis on the plant. "MS" refers 
        to the main stem. "T0", "T1", "T2",..., refers to the primary tillers. "T0.0", 
        "T0.1", "T0.2",..., refers to the secondary tillers of the primary tiller "T0". 
        "T0.0.0", "T0.0.1", "T0.0.2",..., refers to the tertiary tillers of the secondary 
        tiller T0.0. See :ref:`botanical_positions`. 
    * - **N_phytomer**
      - The total number of vegetative phytomers formed on the axis.
    * - **HS_final**
      - The Haun Stage at the end of growth of the axis.
    * - **TT_stop_axis**
      - If the axis dyes: thermal time (since crop emergence) of end of growth. If the axis grows up to flowering:  *NA*  
    * - **TT_del_axe**
      - If the axis dyes: thermal time (since crop emergence) of disappearance. If the axis grows up to flowering:  *NA*  
    * - **id_dim**
      - key (int) linking to DimTable. id_dim allows referring to the data that describe the dimensions of the phytomers of the axis
    * - **id_phen**
      - key (int) linking to PhenTable. id_phen allows referring to the data that describe the phenology of the axis
    * - **id_ear**
      - Key (int) linking to EarTable. id_ear allows referring to the data that describe the ear of the axis. 
        For the regressive axes, **id_ear**="NA".
    * - **TT_em_phytomer1**
      - Thermal time (relative to canopy emergence) of tip appearance of the first true leaf (not coleoptile or prophyll)
    * - **TT_col_phytomer1**       
      - Thermal time (relative to canopy emergence) of collar appearance of the first true leaf                              
    * - **TT_sen_phytomer1**
      - Thermal time (relative to canopy emergence) of full senescence of the first true leaf (this is : thermal time when SSI= 1)
    * - **TT_del_phytomer1**
      - Thermal time (relative to canopy emergence) of disappearance of the first true leaf
       

.. _botanical_positions:

.. figure:: ./image/botanical_positions.png
   :width: 100%
   :align: center

   Botanical position of the axis on the plant. 

See :download:`an example of axeT <../../test/data/test_plantgen/min_min/axeT.csv>`.
   

.. _dimT:

dimT : *description of the dimensions of leaf blades, sheath and internodes*
------------------------------------------------------------------------------

The table allows to describe a number of profiles of dimension, each profile being associated to a value of id_dim.
Dimensions of organs can be given for each of the axis mentioned in :ref:`axeT <axeT>` (using different values of id_dim for different axes) or dimensions can be given by categories of axes (axe of a given category should share a same value of id_dim and will have the same profile of dimensions).

Positions on an axis are expressed as relative position (index_rel_phytomer = phytomer rank/N_phytomer);

Use of relative position makes it possible to use a same profile of dimension for axes differing in the final number of phytomers (N_phytomer);
Use of relative position makes it possible to document a profile with only some the phytomers on an axis:
Missing data will be estimated by linear interpolation according to index_rel_phytomer;  

Actual dimension of the blade sheath and internode axis are hence calculated using id_dim and N_phytomer.

There is one line per phytomer documented.

Columns are :

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **id_dim**
      - the identifier referred to in :ref:`axeT <axeT>`. By convention, if the current **id_dim** 
        ends by ``0`` (e.g. **id_dim**=``1110``), then the current line documents 
        the dimensions of a regressive axis. If the current **id_dim** ends by 
        ``1`` (e.g. **id_dim**=``1111``), then the current line documents the 
        dimensions of a non-regressive axis.
    * - **index_rel_phytomer** 
      - The relative phytomer position 
    * - **L_blade**
      - length of the mature blade (cm)
    * - **W_blade**
      - Maximum width of the mature leaf blade (cm)
    * - **L_sheath** 
      - Length of a mature sheath (cm)
    * - **W_sheath** 
      - Diameter of the stem or pseudo stem at the level of sheath (cm)
    * - **L_internode** 
      - Length of an internode (cm)
    * - **W_internode** 
      - Diameter of an internode (cm)
      
See :download:`an example of dimT <../../test/data/test_plantgen/min_min/dimT.csv>`.


.. _phenT:

phenT : *description of phenology of axes*
-----------------------------------------------------------------

PhenT controls the dynamics of leaf appearance, ligulation, senescence and disappearance.
Internal rules of Adel coordinate sheaths and internodes to the blades so that :ref:`phenT <phenT>` controls indirectly the whole dynamics of plant development.

Positions on an axis are expressed as relative positions.

One timing of development has to be documented for each value taken by id_phen in :ref:`axeT <axeT>`; axes sharing a same value of id_phen will share the same timing;
Use of relative position makes it possible to use a same developmental timing for axes differing in the final number of phytomers;
Use of relative position makes it possible to document a developmental timing with a number of value higher than the number of phytomers on an axis:
this is required because the dynamics of SSI shows a complex behavior(see below)

Timing of developmental events on a leaf is given relative to the timing of the event on leaf 1 of the axis;
Actual timing is computed from :ref:`phenT <phenT>` and the data concerning leaf 1 in :ref:`axeT <axeT>`. 

There is one line per value of index_rel_phytomer documented. For a smooth description of the dynamics of SSI from crop emergence to maturity,
approximately 40 values of index_rel_phytomer should be documented (for each value of id_phen).
More over for each value of id_phen, one line should be documented for index_rel_phytomer = 0, so as to allow interpolation.

Columns are :

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **id_phen** 
      - the index referred to in :ref:`axeT <axeT>`
    * - **index_rel_phytomer** 
      - normalized phytomer position, starting from index_rel_phytomer = 0
    * - **dTT_em_phytomer** 
      - Thermal time of the appearance of the tip of leaf out of the whorl made by the older blade; expressed as thermal time since TT_em_phytomer1
    * - **dTT_col_phytomer*
      - Thermal time of the appearance of collar; expressed as thermal time since TT_col_phytomer1
    * - **dTT_sen_phytomer** 
      - Thermal time for which SSI = n (where n is the phytomer rank); expressed as thermal time since TT_sen_phytomer1
    * - **dTT_del_phytomer** 
      - Thermal time after which the leaf blade is destroyed and is not displayed in the 3D mock-up anymore; expressed as thermal time since TT_del_phytomer1

See :download:`an example of phenT <../../test/data/test_plantgen/min_min/phenT.csv>`.


.. _earT:

earT : *description of ear dimension and phenology*
----------------------------------------------------

There is one line per ear type (referred by id_ear in :ref:`axeT <axeT>`)

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **id_ear** 
      - the identifier referred to in :ref:`axeT <axeT>`
    * - **dTT_ap_ear** 
      - Thermal time interval between flag leaf ligulation and ear appearance (appearance of the tip of highest spike, discounting the awn)
    * - **dTT_ap_peduncle** 
      - Thermal time interval between flag leaf ligulation and peduncle appearance (appearance of the base of the ear) 
    * - **TT_z92** 
      - Thermal time (relative to canopy emergence) of the end of grain filling (corresponding on z92 on Zadoks scale)
    * - **L_peduncle** 
      - length of the ear peduncle (cm)
    * - **W_peduncle** 
      - diameter of the ear peduncle (cm)
    * - **L_ear** 
      - length of the ear without awns (cm)
    * - **A_ear** 
      - projected area of ear without awn  (cm2)
    * - **L_spike** 
      - Total length of the spike : from base of the ear to the top of the awns (cm)    


.. _ssi2sen:

ssi2sen : *description of progression of senescence in upper leaf blades as a function of SSI*
-----------------------------------------------------------------------------------------------

Adel considers two categories of phytomers for describing the progression of senescence in leaf blades.

* for lower leaves, the senescence progresses linearly as function of SSI and blades sequentially: the senescence of blade at rank n starts when senescence of blade n-1 has finished. 
  This means that the senesced fraction of leaf n is : 1+SSI -n. It depends only in ssi and there is no need for additional parameters.
* for upper leaves, the progress of senescence is more complex and several leaf blades senesce simultaneously: 
  SSi2senT contains data to calculate the fraction of senesced area of each upper leaves as function of ssi.

The upper leaves correspond approximately to the leaves beard by an elongated internode. 
The number of lower leaves showing a linear progress of senescence is called Nsenlow;
The number of upper leaves showing a complex progress is called Nsenup

All upper leaf blades start to senesce at the same time, that is at :math:`SSI = Nsenlow`; 
Senescence of each upper leaf blade progresses first at a slow rate,identical for all leaves, then at a fast rate.

The parameter used to describe these kinetics are the value of the slow rate (R_sen1), the value of ssi (dssit1) at the onset of fast senescence 
and the value of SSI (dssit2) at full senescence for each upper leaf. 

The table defines the parameter values for the upper leaves.
There is one line per upper leaf and the number of lines of the file must be Nsenup
The values d_SSIt1 and dssit2 are specified in term of difference with the ssi at onset of upper leaves senecence (Nsenlow)

It should be noted that the present description of progress of senescence is over-parameterized, resulting in a constraint between parameters value.
This comes from the fact that at any time the sum of the rate of progress of senescence for all leaves should be one. 
Complying with this constraint is not straightforward. So a user that do not know precisely the value of parameters in his experiment should probably use the default values to ensure a consistent behavior.


.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **N_senup**
      - Number of leaves that show two phases during senescence (the value is repeated for all lines!)
    * - **R_sen1**
      - Rate of progress of senescence during phase 1 (the value is repeated for all lines !)
    * - **dssit1**
      - (SSI when the leaf blade starts phase 2) - Nsenlow)
    * - **dssit2**  
      - (SSI when the leaf blade is 100% senesced - Nsenlow)



Inputs describing geometry
*****************************

Input are required to define the geometry of leaves (normalized 2D shape, midrib curvature and azimuth) and the geometry of stems (inclination, azimuth)

Normalized 2D shapes are leaf width variations with distance to the base of the leaf, both axes being normalized so that max values is 1.

Normalized 2D shapes and midrib curvature are stored as collections and Adel will draw and individual leaf by scaling a 2D shape plus taking a midrib curvature from these collections. 

The inclination of axes is defined by two parameters DredT and Tillerinc.
DredT represents the horizontal distance between the main stem and a tiller at flowering.
Tillerinc represents the angle of insertion of a tiller at flowering.
When a tiller grows, it starts with angle of 3° compared to the vertical. Then, during the period of extension of the lower internode, insertion angle increases up to the value Tillerinc.
It will keep this value until the top of the stem reaches the distance DredT from the main stem. When this is reach, 
the two upper visible nodes rotate so that the top of the tillers remains at distance DredT. Any internode that elongate
later is vertical. Note that when sheath disappear, new node become visible and will become involved in the process.

genGoeaxe (see below) includes a parameter to randomly tilt the main stem of a small value around the vertical. When the main stem is tilted, all the plant follows


How to specify the geometric data
---------------------------------

The collections for 2D leaf shape and for leaf curvature should be specified as one list of lists of matrices for 2D shape and one list of matrices for midrib curvature.

* the first level in the list is for collection index 
* the second level is for matrix index.

See alea for more information.

Besides these collections, R functions should be provided as inputs. A first list of function is for defining the axis geometry;
A second list of functions is for selecting shapes in the collections mentioned above.

The first list should provide 3 R functions of axis number (0 = main stem) that return:
    * **azT** : the azimuth(deg) of the first leaf of the axis with reference to the azimuth of the parent leaf
    * **incT** : the inclination (deg) of the base of the tiller compared with main stem
    * **dredT** : the distance (at maturity) between tiller and main stem

These functions can be generated by the predefined *genGeoAxe* node or be freely user-defined in a *freeGeoAxe* node.

In genGeoAxe 
The azimuth of a tiller stem is the same as that of the axilling main stem leaf. 
The azimuth of the first leaf of a primary tiller is with an angle of 75° relatively to that of the axilling main stem leaf.
For secondary tillers, the azimuth of the first leaf is also with a fixed angle relatively to that of the parent tiller.


A sample code of "geoAxe" function is:                                              

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


The second list should provide two Rfunctions for drawing in the collections of leaf shape

Inputs have to be axis number, leaf position, leaf position counted from top, and leaf stage, defined as current length/final length. 
Returned values have to be :

    * **azim** : the azimuth (deg) of the leaf compared to the previous one
    * **Lindex** : the index of the collection to use for leaf curvature

These functions can be generated by the predefined genGeoLeaf node or be freely user-defined in a *freeGeoLeaf* node:

A sample code for a "geoLeaf" function is be : 

.. code-block:: r

    geoLeaf <- list(
        Azim = function(a,n,ntop) {0 * runif(1)},
        Lindex = function(a,n,ntop,stage) {ntop + 1}
        

Inputs describing simulation
********************************

Time step is given as a list of values of thermal times for which a mock-up is to be produced.
Positions of plants within the plot are given externally from adel to a planter.


Some revelations about deep ADEL
=================================

They are in Adel a number of secret options that are very important and for this reason have never been documented before

Coordination in organ extension
*******************************

The thermal time of leaf tip appearance and leaf collar appearance given in Phen_T are used to calculate a number of features;
- the leaf extension (blade + sheath) is simulated as starting 0,4 phyllochron between tip appearance, and having a constant rate (cm.°C-1.J-1) for a duration of 2 phyllochrons
- The model calculate the length of the hidden part of a leaf (whorl length) : at tip emergence, this hidden length is the blade length; 
at collar emergence this hidden length is taken as the length of sheath n-1; Between it is approximated by linear interpolation. 
This is used to calculate the length of the visible part of the leaf in the post processing treatments. Note that this calculation is not fully accurate because sheath n-1 stop growing before collar n emerges


The leaf extension is simulated as consisting sequentially of the blade extension, followed by the sheath extension. 

The internode extension is simulated as following sequentially the sheath extension, and taking place at a constant rate, for a duration of 1/(stemleaf) phyllochron
It is known that in grass, internode fast extension start at collar emergence. However there is no such calculation of collar emergence in the model: 
it expected that the synchronization with collar emergence will be reasonably well approximated by the synchronization implemented with the end of leaf extension.


The parameters for these coordinations are defined in AdelRunOption, which remained to be documented


Senescence of sheath and internodes
***********************************

The senescence of sheath n is simulated as being synchronous with the senescence of blade n+2
The disappearance of sheath n is simulated as synchronous with disappearance of blade n+1

There is no senescence implemented for internodes : they stay green.
For ear and peduncle : to be documented

On regressing tillers, individual leaf senescence is simulated from SSI with the same pattern as on non-regressing tillers.


Disappearance of dead tillers
*****************************

A dead tiller can be programmed to disappear some time after it stops growing. 
Only the blades and sheaths, not the internodes, disappear. This will be changed in further version, so that internode also disappear
When this happens, it has priority over the process of disappearance following leaf senescence. 


Description of Adel's outputs
==============================


.. _plantgen:

Construction of the input tables 
=================================

Authors: Mariem ABICHOU, Bruno ANDRIEU and Camille CHAMBON

ADEL requires inputs characterizing the development of plants as described 
in :ref:`development_input`.

The :mod:`plantgen <alinea.adel.plantgen>` package allows the user who does not have 
a complete set of data to estimate the missing inputs. 
Inside this package, the :mod:`plantgen <alinea.adel.plantgen.plantgen>` module 
provides routines to construct :ref:`axeT <axeT>`, :ref:`dimT <dimT>` and :ref:`phenT <phenT>`. 
It provides also some other tables for debugging purpose.

We have considered three possible levels of completeness of data, denote as MIN, 
SHORT, and FULL. In the next subsections, we 

* describe the levels of completeness of the data and of the parameters set 
  by the user,
* describe how to construct the inputs of ADEL from a Python interpreter, 
  using the routine :func:`gen_adel_input_data <alinea.adel.plantgen.plantgen.gen_adel_input_data>`. 
  This routine can be used whatever the level of completeness of the raw inputs, 
  adapting the processing automatically,
* describe how to construct the inputs of ADEL from the Visualea interface, 
  using one of the following routines:
  
  * :func:`gen_adel_input_data_from_min <alinea.adel.plantgen.plantgen.gen_adel_input_data_from_min>`
  * :func:`gen_adel_input_data_from_short <alinea.adel.plantgen.plantgen.gen_adel_input_data_from_short>`
  * :func:`gen_adel_input_data_from_full <alinea.adel.plantgen.plantgen.gen_adel_input_data_from_full>`.
        
All routines belong to :mod:`plantgen <alinea.adel.plantgen.plantgen>`.
All routines produce the same output tables: 

* :ref:`axeT <axeT>`
* :ref:`dimT <dimT>`
* :ref:`phenT <phenT>`
* :ref:`phenT_abs <phenT_abs>`: the equivalent of :ref:`phenT <phenT>`, but 
  with absolute dates and absolute positions.
* :ref:`dimT_abs <dimT_abs>`: the equivalent of :ref:`dimT <dimT>`, but with 
  absolute positions.
* :ref:`dynT <dynT>`: the dynamic of the leaves for each type of axis. 
* :ref:`phenT_first <phenT_first>`: a subset of :ref:`phenT_abs <phenT_abs>`, 
  containing only the lines of :ref:`phenT_abs` which correspond to the first 
  phytomer of each axis.
* :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`: the dynamic of *HS*, *GL* and *SSI* when 
  *TT* varies, for each type of axis. 
* :ref:`tilleringT <tilleringT>`: the dynamic of tillering.
* :ref:`cohortT <cohortT>`: the theoretical and the simulated cardinalities of 
  each cohort.

.. _levels_of_completeness:

The levels of completeness
***************************

The information provided to generate Adel input must be provided in two tables: 
``dynT_user`` and ``dimT_user``. ``dynT_user`` and ``dimT_user`` can  have 
different  levels  of  completeness:  ``FULL``,  ``SHORT`` and  ``MIN``.  
According  to  their  level  of completeness, ``dynT_user`` and ``dimT_user`` 
take different types, shapes and/or contents.

The table below list the specific designation in :func:`plantgen <alinea.adel.plantgen>`
for ``dynT_user``  and ``dimT_user`` for each level of completeness:

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
      
.. seealso:: the documentation of :class:`DataCompleteness <alinea.adel.plantgen.plantgen.DataCompleteness>`
             in the :ref:`adel_reference`.                       
      
.. _construct_inputs_from_interpreter:

Construction of Adel input tables using the Python interpreter
***************************************************************

:func:`gen_adel_input_data <alinea.adel.plantgen.plantgen.gen_adel_input_data>` 
is aimed to be used from Python interpreter.

First we explain the arguments of :func:`gen_adel_input_data <alinea.adel.plantgen.plantgen.gen_adel_input_data>` 
that the user has to define. Second we present a complete code example to use 
:func:`gen_adel_input_data <alinea.adel.plantgen.plantgen.gen_adel_input_data>` 
from a Python interpreter.          

The arguments to define by the user
-------------------------------------

The arguments to define are:

* dynT_user : *the leaf dynamic parameters set by the user*

  *dynT_user* can be either a :class:`pandas.DataFrame` or a :class:`dict`, 
  depending on the argument :ref:`*dynT_user_completeness* <levels_of_completeness>`. 

* dimT_user : *the dimensions of the axes set by the user*

  *dimT_user* is a :class:`pandas.DataFrame`, which content depends on 
  :ref:`*dimT_user_completeness* <levels_of_completeness>`.

* dynT_user_completeness and dimT_user_completeness : *the levels of completeness of dynT_user and dimT_user*

  :ref:`*dynT_user_completeness* <levels_of_completeness>` and :ref:`*dimT_user_completeness* <levels_of_completeness>` have to be consistent with 
  respectively *dynT_user* and *dimT_user*.

* *plant_number*, *decide_child_cohort_probabilities*, *MS_leaves_number_probability_distribution*, ...

  The other arguments of the routine are: 
    
  * *plant_number*, the number of plants to be generated,
  * *decide_child_cohort_probabilities*, for each cohort the probability of 
    emergence of an axis when the parent axis is present, 
  * *MS_leaves_number_probability_distribution*, the probability distribution 
    of the final number of main stem leaves,
  * *TT_bolting*, the thermal time at which the bolting starts,
  * *final_axes_density*, the final number of axes which have an ear, per square meter,
  * *GL_number*, the thermal times of GL measurements and corresponding values of green leaves number, 
  * *delais_TT_stop_del_axis*, the thermal time between an axis stop growing and its disappearance,
  * and *TT_col_break*, the thermal time when the rate of progress Haun Stage vs thermal time is changing. 
    If phyllochron is constant, then *TT_col_break* is null.
  
Code example
-------------
  
Now let's see a complete code example to use 
:func:`gen_adel_input_data <alinea.adel.plantgen.plantgen.gen_adel_input_data>` 
from a Python interpreter::
    
    # define the levels of completeness
    from alinea.adel.plantgen.plantgen import DataCompleteness
    dynT_user_completeness = DataCompleteness.SHORT
    dimT_user_completeness = DataCompleteness.SHORT
    
    # import the pandas library. In this example, pandas is used to read and 
    # write the tables.
    import pandas

    # read the dynT_user_SHORT table. "dynT_user_SHORT.csv" must be in the working directory. 
    # "dynT_user_SHORT.csv" must be coherent with dynT_user_completeness.
    dynT_user = pandas.read_csv('dynT_user_SHORT.csv')
        
    # read the dimT_user_SHORT table. "dimT_user_SHORT.csv" must be in the working directory.
    # "dimT_user_SHORT.csv" must be coherent with dimT_user_completeness.
    dimT_user = pandas.read_csv('dimT_user_SHORT.csv')    
    
    # define the other arguments
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
    final_axes_density = 250
    GL_number = {1117.0: 5.6, 1212.1:5.4, 
                 1368.7:4.9, 1686.8:2.4, 
                 1880.0:0.0}
    delais_TT_stop_del_axis = 600
    TT_col_break = 0.0
    
    # launch the construction
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
                                   final_axes_density, 
                                   GL_number, 
                                   delais_TT_stop_del_axis, 
                                   TT_col_break, 
                                   dynT_user_completeness, 
                                   dimT_user_completeness)

    # write axeT, dimT and phenT to csv files in the working directory, replacing
    # missing values by 'NA' and ignoring the indexes (the indexes are the labels of
    # the lines). 
    axeT.to_csv('axeT.csv', na_rep='NA', index=False)
    dimT.to_csv('dimT.csv', na_rep='NA', index=False)
    phenT.to_csv('phenT.csv', na_rep='NA', index=False)
    
    # "axeT.csv", "dimT.csv" and "phenT.csv" are now ready to be used by Adel.
    
    
.. _construct_inputs_from_visualea:

Construction of Adel input tables using Visualea
*************************************************

The following routines allow to construct the inputs of ADEL: 

* :func:`gen_adel_input_data_from_min <alinea.adel.plantgen.plantgen.gen_adel_input_data_from_min>`: 
  construct the inputs of ADEL from :ref:`dynT_user_MIN` and :ref:`dimT_user_MIN`,
* :func:`gen_adel_input_data_from_short <alinea.adel.plantgen.plantgen.gen_adel_input_data_from_short>`: 
  construct the inputs of ADEL from :ref:`dynT_user_SHORT` and :ref:`dimT_user_SHORT`,  
* and :func:`gen_adel_input_data_from_full <alinea.adel.plantgen.plantgen.gen_adel_input_data_from_full>`: 
  construct the inputs of ADEL from :ref:`dynT_user_FULL` and :ref:`dimT_user_FULL`.
    
All these routines belong to :mod:`alinea.adel.plantgen.plantgen`.

These routines are wrapped in the following Visualea nodes:

.. list-table::
    :widths: 10 10 10
    :header-rows: 1

    * - ``plantgen_MIN``
      - ``plantgen_SHORT``
      - ``plantgen_FULL``
    * - .. image:: image/plantgen_MIN_node.png
      - .. image:: image/plantgen_SHORT_node.png
      - .. image:: image/plantgen_FULL_node.png
    * - .. image:: image/plantgen_MIN_widget.png
      - .. image:: image/plantgen_SHORT_widget.png
      - .. image:: image/plantgen_FULL_widget.png

The following table summarizes the nodes, the routines and the levels of completeness 
of :ref:`dynT <dynT>` and :ref:`dimT <dimT>`:

.. list-table::
    :widths: 15 50 20
    :header-rows: 1

    * - Level of completeness
      - Convenience routine
      - Visualea node
    * - **MIN** 
      - :func:`gen_adel_input_data_from_min <alinea.adel.plantgen.plantgen.gen_adel_input_data_from_min>`
      - ``plantgen_MIN``
    * - **SHORT** 
      - :func:`gen_adel_input_data_from_short <alinea.adel.plantgen.plantgen.gen_adel_input_data_from_short>`
      - ``plantgen_SHORT``
    * - **FULL** 
      - :func:`gen_adel_input_data_from_full <alinea.adel.plantgen.plantgen.gen_adel_input_data_from_full>`
      - ``plantgen_FULL``
 
The following dataflow demonstrates how to use ``plantgen_MIN``, ``plantgen_SHORT``, 
and ``plantgen_FULL`` through Visualea:

.. image:: image/plantgen_dataflow.png

The user must select existing data nodes to set the input and ouput tables.

The following data-flow also demonstrates how to use ``plantgen_MIN`` through 
Visualea:

.. image:: image/plantgen_MIN_csv_dataflow.png
  
In this case the user must give the paths of csv files for inputs and outputs. 
Attention: the paths set in the example will not work on your computer. You have 
to adapt them to your needs. This example is more straightful because you don't 
have to create output data nodes before running, but it is also less portable.     

These dataflows are accessible from the Package explorer of Visualea, in 
``alinea.adel.tutorials.plangen``.


Appendices
***********

The appendices contain the description of the following data:

* :ref:`dynT_user_FULL <dynT_user_FULL>`: the dynamic of the Haun stage for 
  *at least* the most frequent non-regressive axis of each cohort.
* :ref:`dynT_user_SHORT <dynT_user_SHORT>`: the dynamic of the Haun stage for 
  *exactly* the most frequent non-regressive axis of each cohort.
* :ref:`dynT_user_MIN <dynT_user_MIN>`: the dynamic of the Haun stage for 
  the main stem.
* :ref:`dimT_user_FULL <dimT_user_FULL>`: the dimensions of *at least* the 
  most frequent non-regressive axis of each cohort.
* :ref:`dimT_user_SHORT <dimT_user_SHORT>`: the dimensions of *exactly* the 
  most frequent non-regressive axis of each cohort.
* :ref:`dimT_user_MIN <dimT_user_MIN>`: the dimensions of the main stem. 
* :ref:`phenT_abs <phenT_abs>`: the equivalent of :ref:`phenT <phenT>`, but 
  with absolute dates and absolute positions.
* :ref:`dimT_abs <dimT_abs>`: the equivalent of :ref:`dimT <dimT>`, but with 
  absolute positions.
* :ref:`dynT <dynT>`: the dynamic of the leaves for each type of non-regressive 
  axis. 
* :ref:`phenT_first <phenT_first>`: a subset of :ref:`phenT_abs <phenT_abs>`, 
  containing only the lines of :ref:`phenT_abs` which correspond to the first 
  phytomer of each non-regressive axis.    
* :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`: the dynamic of *HS*, *GL* and *SSI* when 
  *TT* varies, for each type of non-regressive axis. 
* :ref:`tilleringT <tilleringT>`: the dynamic of tillering.
* :ref:`cohortT <cohortT>`: the theoretical and the simulated cardinalities of 
  each cohort.
    
.. _dynT_user_FULL:

dynT_user_FULL
---------------

:ref:`dynT_user_FULL` is a table which describes the dynamic of the Haun stage for each 
type of non-regressive axis. The type of an axis is defined by its cohort index 
and its final number of leaves. One type of axis is described by one line, which 
contains the following parameters *N_cohort*, *Nff*, *a_cohort*, *TT_col_0*, 
*TT_col_nff*, *n0*, *n1* and *n2*. See :ref:`dynT` for the definition of these 
parameters.

See :download:`an example of dynT_user_FULL <../../test/data/test_plantgen/full_full/dynT_user.csv>`.


.. _dynT_user_SHORT:

dynT_user_SHORT
----------------

:ref:`dynT_user_SHORT` is a table which describes the dynamic Haun stage for each cohort. 
Values are supposed to represent a non-regressive axis having the most frequent 
number of leaves for the cohort. One line refers to one cohort and contains the 
following parameters: *N_cohort*, *a_cohort*, *TT_col_0*, *TT_col_nff*, *n0*, *n1* 
and *n2*. See :ref:`dynT` for a description of these parameters.

See :download:`an example of dynT_user_SHORT <../../test/data/test_plantgen/short_short/dynT_user.csv>`.


.. _dynT_user_MIN:

dynT_user_MIN
--------------

:ref:`dynT_user_MIN` is a dictionary which describes the dynamic of the Haun stage for 
the main stem. The dictionary contains the following keys: *a_cohort*, *TT_col_0*, 
*TT_col_nff*, *n0*, *n1* and *n2*. See :ref:`dynT` for a description of these 
parameters.

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

:ref:`dimT_user_FULL` is the same as :ref:`dimT_abs`, except that :ref:`dimT_user_FULL` 
does not contain the dimensions of the regressive axes. Thus, each **id_dim** of 
:ref:`dimT_user_FULL` refers implicitly to the corresponding non-regressive axis. 
For example, **id_dim**=``111`` in :ref:`dimT_user_FULL` corresponds to **id_dim**=``1111`` 
in :ref:`dimT_abs`.   

See :download:`an example of dimT_user_FULL <../../test/data/test_plantgen/full_full/dimT_user.csv>`.


.. _dimT_user_SHORT:

dimT_user_SHORT
----------------

:ref:`dimT_user_SHORT <dimT_user_SHORT>` is a table which contains one profile of dimensions of the organs 
for each cohort. Values represent a non-regressive axis having the most frequent 
leaf number for that cohort. Each line contains the following data: *id_axis*, *index_phytomer*, 
*L_blade*, *W_blade*, *L_sheath*, *W_sheath*, *L_internode* and *W_internode*. 
*id_axis* is the index of the cohort. See :ref:`dimT_abs` for a description 
of the other data.

See :download:`an example of dimT_user_SHORT <../../test/data/test_plantgen/short_short/dimT_user.csv>`.
        

.. _dimT_user_MIN:

dimT_user_MIN
--------------

:ref:`dimT_user_MIN <dimT_user_MIN>` is a table which contains the dimensions of the organs, for each 
phytomer of the main stem. Values are given only for a main stem having the 
most frequent number of phytomers. Each line contains the following data: 
*index_phytomer*, *L_blade*, *W_blade*, *L_sheath*, *W_sheath*, *L_internode* 
and *W_internode*. See :ref:`dimT_abs` for a description of these data.

See :download:`an example of dimT_user_MIN <../../test/data/test_plantgen/min_min/dimT_user.csv>`.


.. _phenT_abs:

phenT_abs
----------

:ref:`phenT_abs` is an intermediate table used to construct :ref:`phenT <phenT>`. 
This table is not an input of Adel. Thus the user normally needn't it. This table 
can be useful for debugging.

:ref:`phenT_abs` is the same as :ref:`phenT <phenT>`, except that:
    * the positions of the phytomers are not normalized,
    * the dates of developmental events are absolute.

See :download:`an example of phenT_abs <../../test/data/test_plantgen/min_min/phenT_abs.csv>`.
        

.. _dimT_abs:

dimT_abs
----------

:ref:`dimT_abs` is an intermediate table used to construct :ref:`dimT <dimT>`. 
This table is not an input of Adel. Thus the user normally needn't it. This table 
can be useful for debugging.

:ref:`dimT_abs` is the same as :ref:`dimT <dimT>`, except that the positions 
of the phytomers are not normalized.

See :download:`an example of dimT_abs <../../test/data/test_plantgen/min_min/dimT_abs.csv>`.


.. _dynT:

dynT
-----

:ref:`dynT` is an intermediate table used to construct :ref:`phenT_abs <phenT_abs>`. 
This table is not an input of Adel. Thus the user normally needn't it. This table 
can be useful for debugging.

:ref:`dynT` is a table which contains the dynamic of the leaves, for each type 
of non-regressive axis. The type of an axis is defined by its cohort index and 
its final number of leaves. 
There is one line per type of axis. Each line contains the following data:

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **N_cohort** 
      - the index of the cohort to which belongs *id_axis*
    * - **id_axis** 
      - the current type of axis. This type is the concatenation of *N_cohort* 
        and *Nff*.
    * - **cardinality**
      - the cardinality of the set composed of *id_axis*
    * - **Nff** 
      - the final number of leaves of *id_axis*
    * - **a_cohort** 
      - the rate of Haun Stage vs Thermal time. This is the rate of the 
        first phase in case of bilinear behavior.
    * - **TT_col_0** 
      - the thermal time for Haun Stage equal to 0
    * - **TT_col_break**
      - the thermal time when the rate of phytomers emergence is changing
    * - **TT_col_nff** 
      - the thermal time when Haun Stage is equal to *Nff*
    * - **n0** 
      - number of green leaves at *t0*
    * - **n1** 
      - number of green leaves at *t1*
    * - **n2** 
      - number of green leaves at *TT_col_nff*
    * - **t0**
      - the thermal time at the start of leaf senescence 
    * - **t1**
      - the thermal time at which the senescence starts
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

The lines are ordered by cohort index (i.e. **N_cohort**), and, within each cohort 
index, by **cardinality**.   

See :download:`an example of dynT <../../test/data/test_plantgen/min_min/dynT.csv>`.
        

.. _phenT_first:

phenT_first
------------

:ref:`phenT_first` is an intermediate table used to construct :ref:`phenT <phenT>` and 
:ref:`axeT <axeT>`. This table is not an input of Adel. Thus the user normally 
needn't it. This table can be useful for debugging.

:ref:`phenT_first` contains only the lines of :ref:`phenT_abs` which correspond to 
the first phytomer of each non-regressive axis. These lines have *index_phytomer* 
equal to 1.

See :download:`an example of phenT_first <../../test/data/test_plantgen/min_min/phenT_first.csv>`.


.. _HS_GL_SSI_T:

HS_GL_SSI_T
------------

:ref:`HS_GL_SSI_T` is constructed for debugging purpose.    

:ref:`HS_GL_SSI_T` describes, for each type of non-regressive axis, the dynamic 
of *HS*, *GL* and *SSI* when *TT* varies. The type of an axis is defined by its 
cohort index and its final number of leaves.

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

   For each axis, *TT* varies from 0 to :attr:`alinea.adel.plantgen.params.TT_DEL_FHAUT`.     

See :download:`an example of HS_GL_SSI_T <../../test/data/test_plantgen/min_min/HS_GL_SSI_T.csv>`.


.. _tilleringT:

tilleringT
------------

:ref:`tilleringT` is constructed for debugging purpose.

:ref:`tilleringT` describes the dynamic of tillering. It stores the number of axes at 
important dates: the start of growth, the thermal time of the bolting, and the thermal 
time of the flowering.

.. list-table::
    :widths: 10 50
    :header-rows: 1

    * - Column
      - Description
    * - **TT** 
      - the thermal time.
    * - **NbrAxes** 
      - the number of axes.

See :download:`an example of tilleringT <../../test/data/test_plantgen/min_min/tilleringT.csv>`.


.. _cohortT:

cohortT
------------

:ref:`cohortT` is constructed for debugging purpose.

:ref:`cohortT` describes the theoretical and the simulated cardinalities of 
each cohort. It permits the user to validate the simulated cardinalities against 
the theoretical ones. The theoretical cardinalities are calculated from the 
probabilities of emergence of an axis when the parent axis is present. These 
probabilities are given by the user. The simulated cardinalities are calculated 
for each plant using :func:`alinea.adel.plantgen.tools.decide_child_cohorts`.

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

See :download:`an example of cohortT <../../test/data/test_plantgen/min_min/cohortT.csv>`.

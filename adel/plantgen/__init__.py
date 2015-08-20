# -*- python -*-
#
#       Adel.PlantGen
#
#       Copyright 2012-2014 INRIA - CIRAD - INRA
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
"""
:mod:`plantgen <alinea.adel.plantgen>` permits to generate the :class:`pandas.DataFrame` which contain the 
plant data to be used as input for generating plot with ADEL.

Available submodules are:

* :mod:`plantgen_interface <alinea.adel.plantgen.plantgen_interface>`:
    Front-end for the generation of the :class:`pandas.DataFrame` which contain the plant data 
    expected by ADEL. Uses the routines of the modules :mod:`plantgen_core <alinea.adel.plantgen.plantgen_core>` 
    and :mod:`tools <alinea.adel.plantgen.tools>`, and uses the parameters of the 
    module :mod:`params <alinea.adel.plantgen.params>`.
* :mod:`plantgen_core <alinea.adel.plantgen.plantgen_core>`: 
    Routines defining the main steps of the process: :func:`init_axes <alinea.adel.plantgen.plantgen_core.init_axes>`, 
    :func:`phenology_functions <alinea.adel.plantgen.plantgen_core.phenology_functions>`, 
    :func:`plants_structure <alinea.adel.plantgen.plantgen_core.plants_structure>`, 
    :func:`organs_dimensions <alinea.adel.plantgen.plantgen_core.organs_dimensions>`, 
    :func:`axes_phenology <alinea.adel.plantgen.plantgen_core.axes_phenology>`, 
    :func:`init_axes <alinea.adel.plantgen.plantgen_core.phenology_functions>`, 
    :func:`init_axes <alinea.adel.plantgen.plantgen_core.phenology_functions>`.
* :mod:`tools <alinea.adel.plantgen.tools>`: 
    Generic routines used in the :mod:`plantgen <alinea.adel.plantgen>` package. 
    These routines can also be used by other packages.
* :mod:`alinea.adel.plantgen.params`: 
    The constant parameters used in :mod:`plantgen <alinea.adel.plantgen>`.
    
Authors: M. Abichou, B. Andrieu, C. Chambon
"""
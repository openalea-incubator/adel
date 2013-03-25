"""
:mod:`alinea.adel.plantgen` permits to generate the :class:`pandas.DataFrame` which contain the 
plant data to be used as input for generating plot with ADEL.

Available submodules are:

* :mod:`alinea.adel.plantgen.plantgen`:
    Front-end for the generation of the :class:`pandas.DataFrame` which contain the plant data 
    expected by ADEL. Uses the routines of the modules listed below.
* :mod:`alinea.adel.plantgen.axeT`: 
    Routines for the generation of the :ref:`axeT <axeT>` :class:`pandas.DataFrame`.
* :mod:`alinea.adel.plantgen.dynT`: 
    Routines for the generation of the :ref:`dynT <dynT>` :class:`pandas.DataFrame`.
* :mod:`alinea.adel.plantgen.dimT`: 
    Routines for the generation of the :ref:`dimT <dimT>` :class:`pandas.DataFrame`.
* :mod:`alinea.adel.plantgen.phenT`: 
    Routines for the generation of the :ref:`phenT <phenT>` :class:`pandas.DataFrame`.
* :mod:`alinea.adel.plantgen.params`: 
    The constant parameters used in :mod:`alinea.adel.plantgen`.
* :mod:`alinea.adel.plantgen.tools`: 
    Generic routines used in the :mod:`alinea.adel.plantgen` package. These routines can also be 
    used by other packages.
    
Authors: Mariem Abichou, Camille Chambon, Bruno Andrieu
"""
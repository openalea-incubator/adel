============
Alinea.Adel
============

**Authors** : C. Fournier, C. Pradal, B. Andrieu

**Contributors** : 
  * M Abichou (wheat parameteristion, documentation, plantgen), 
  * C Chambon (plantgen, documentation, povray, fit, io widgets)

**Institutes** : INRA, CIRAD

**Status** : R and Python package (+ L-system)

**License** : Cecill-C

About
------

Description
============

Alinea.Adel (Architectural model of DEvelopment based on L-systems) allows
to simulate the 3D architectural development of the shoot of gramineaous plant. 




Content
========

The package hosts generic data structure and simulation tools for gramineaous plants(Fournier & Pradal, unpublished),
the Adel-Maize (Fournier & Andrieu, 1998), Adel-Wheat (Fournier et al. 2003) models, 
together with the wheat parameterisation model of Abichou et al. (2013) and the plastic leaf model of Fournier & Pradal (2012)


Installation
=============

  python setup.py install
  
Requirements
============

* OpenAlea.Deploy
* OpenAlea.Mtg
* OpenAlea.core
* OpenAlea.visualea
* OpenAlea.PlantGl
* NumPy
* Scipy
* MatplotLib
* Pandas
* Rpy2

Installation with conda or mamba
================================

Create an environment (either use mamba or conda):
  
  mamba create -n adel -c conda-forge -c openalea3 alinea.adel -y
  

Activate the environment:

  conda activate adel

Installation from source
========================

Create an environment 

  mamba create -n adel -c openalea3 -c conda-forge alinea.adel --only-deps -y

Install adel

.. code-block:: console

   git clone https://github.com/openalea-incubator/adel.git 
   cd adel
   python setup.py develop


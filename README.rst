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

Installation with conda
=======================

Create an environment:

  conda create -n adel -c openalea openalea.plantgl boost=1.66

Activate the environment:

  [source] activate adel

Install the different packages

  conda install -c conda-forge rpy2 

  conda install -c anaconda -c openalea openalea.mtg openalea.visualea openalea.core numpy scipy matplotlib pandas


Install adel

  python setup.py develop

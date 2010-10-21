====== Caribu ======

**Authors** :  Michael Chelle for models, Christian Fournier and Christophe Pradal for Python and Visualea interfaces 

**Institutes** : INRA / INRIA 

**Status** : Python Package 

**License** : INRA Lisence agreement

**URL** : http://openalea.gforge.inria.fr

===== About =====

=== Description ===

Caribu is a set of tolls for mdelling light distribution on 3D canopies. It is based on the nested radiosity model developped by M. Chelle in 1998.

=== Content ===

Caribu is a python script that combines four models : periodise, s2v, MCSail and Canestra, that respectively allow for scene infinitisation, scene conversion to 1D layers, light flux computations within layers and nested radiosity computation. It operates on a special scene object called Caribuscene, composed of a list of triangles with optical properties representing the plants, a set of direction and intensities (called light sources) representing the sky and a pattern delimitating the scene, used for infinitisation.
This model is completed with a set of utilities for visualisation (using PlantGL), import of caribuscene from files or MTG, and tools for building ligth sources
from meteorological data

===== Installation =====


* Download the egg and use Visualea Installer -> install other egg
* Alternatively, if you're codeveloping the package, dowload the Caribu directory from openaleapkg and run <code> python setup.py develop </code> from wihin the directory.

=== Download ===

* egg is available from Christian.Fournier@supagro.inra.fr
* sources are on the openaleapkg svn server

=== Requirements ===

* Scons 
* SconsX
* OpenAlea.Deploy
* PlantGl (for visualisation)


=== Utilisation ===

* A tutorial dataflow is availabable in Visualea, together with example data
* A pdf documenting the original models is avalable in the doc directory of the package



===== Snapshot =====

A snapshot of Caribu in action under visualea is available here



	





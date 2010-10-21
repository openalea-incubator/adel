#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import os, sys
from setuptools import setup, find_packages
pj = os.path.join

name = 'Alinea.Caribu'
version= '6.0.6'
namespace = 'alinea'
pkg_root_dir = 'src'


description= 'Python/Visualea interface to Caribu Light model'


long_description= ''' Caribu allows to compute light distribution on a set of triangles. 
It combines a projection algorythm for direct lighting and nested radiosity algorythm for multiple rediffusion. 
It can simulate an infinite canopy arround your scene. 
The package is based on four C binaries (Canestra, s2v, MCsail and periodise) 
which are combined in a python script (caribu). 
Visualisation of outputs is available through PlantGL. '''
author= '''Michael Chelle (Canestra,s2v,MCSail,Periodise), 
Christian Fournier (Python packaging and Visualea nodes), 
Christophe Pradal (interface with PlantGL)'''
author_email= 'chelle@grignon.inra.fr,Christian.Fournier@supagro.inra.fr, christophe.pradal@cirad.fr'
url= ''
license= 'INRA License agreement (see pdf in doc)' 


# dependencies 
setup_requires = ['openalea.deploy']
if ("win32" in sys.platform):
    install_requires = ['VPlants.PlantGL']
else:
    install_requires = []
dependency_links = ['http://openalea.gforge.inria.fr/pi']

# Scons build directory
build_prefix= "build-scons"

#retrieving packages
pkgs = [ pkg for pkg in find_packages(pkg_root_dir) if namespace not in pkg]
top_pkgs = [pkg for pkg in pkgs if  len(pkg.split('.')) < 2]
packages = [ namespace + "." + pkg for pkg in pkgs]
package_dir = dict( [('',pkg_root_dir)] + [(namespace + "." + pkg, pkg_root_dir + "/" + pkg) for pkg in top_pkgs] )
wralea_entry_points = ['%s = %s'%(pkg,namespace + '.' + pkg) for pkg in top_pkgs]

# Call to setup
setup(
    name=name,
    version=version,
    description=description,
    long_description=long_description,
    author=author,
    author_email=author_email,
    url=url,
    license=license,
    packages= packages,	
    package_dir= package_dir,
    namespace_packages = [namespace],
    create_namespaces = True,
    zip_safe= False,
    setup_requires = setup_requires,
    install_requires = install_requires,
    dependency_links = dependency_links,                  
    # Include data
    include_package_data = True,
    package_data = {'' : ['*.can', '*.R', '*.8', '*.opt', '*.light', '*.csv', '*.png','*.pyd', '*.so', '*.dylib']},
    # Binary construction and install
    scons_scripts = ['SConstruct'],
    bin_dirs = {'bin':  build_prefix + '/bin'},
    # extensions
    entry_points = { 'wralea':  wralea_entry_points },
   )



    

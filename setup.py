#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


import os, sys
from setuptools import setup, find_packages
pj = os.path.join


# Setup script

name = 'adel'
namespace = 'alinea'
pkg_name = 'alinea.adel'
src_rep = 'adel'

version= '0.8.1'

description= 'ADEL' 
long_description= ''' '''

author= 'Chrisian Fournier, Christophe Pradal'
author_email= ''
url= ''

license= 'INRA License agreement' 

if("win32" in sys.platform):
    install_requires = []
    setup_requires = install_requires + []
else:
    install_requires = []
    setup_requires = []
    
packages =[pkg_name]+['%s.%s'%(pkg_name,x) for x in find_packages(src_rep)]

# Main setup
setup(
    name="Alinea.Adel",
    version=version,
    description=description,
    long_description=long_description,
    author=author,
    author_email=author_email,
    url=url,
    license=license,
    
    namespace_packages = ["alinea"],
    create_namespaces = True,

    py_modules = [],
    # pure python  packages
    packages= packages,
    # python packages directory
    package_dir= {pkg_name : src_rep},

                   
    # Add package platform libraries if any
    include_package_data=True,
    package_data = {'' : ['*.RData', '*.R', '*.8', '*.h', '*.str','*.txt', '*.l', '*.map', '*.csv', '*.png'],},

    # Add package platform libraries if any
    zip_safe = False,


    # Scripts
    entry_points = { 'wralea': [ 'adel= alinea.adel',] },
 
    # Dependencies
    setup_requires = setup_requires + ['openalea.deploy'],
    install_requires = install_requires,
    dependency_links = ['http://openalea.gforge.inria.fr/pi'],
   )



    

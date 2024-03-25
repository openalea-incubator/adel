#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import os, sys 
from setuptools import setup, find_namespace_packages
#from openalea.deploy.metainfo import read_metainfo
pj = os.path.join

version = '2.0.0'
name = 'alinea.adel'

description= '3D plant simulation of graminae crops'
long_description= 'The Adel package characterise 3D plant development for graminae crops.'

authors= 'Christian Fournier, Christophe Pradal'
authors_email = 'christian fournier at inrae fr'

url = 'https://github.com/openalea-incubator/adel'
license = 'Cecill-C'

# Main setup
setup(
    name=name,
    version=version,
    description=description,
    long_description=long_description,
    author=authors,
    author_email=authors_email,
    url=url,
    license=license,
    #py_modules = [],
    # pure python  packages
    packages= find_namespace_packages(where='src', 
                                      include=['alinea.*']),
    # python packages directory
    package_dir= {'': 'src'},


    # Namespace packages creation by deploy
    # Add package platform libraries if any
    include_package_data=True,
    package_data = {'' : ['*.RData', '*.R', '*.8', '*.h', '*.str','*.txt', '*.l', '*.map', '*.csv', '*.png'],},
    share_dirs = {pj(*('alinea', 'adel', 'data')): pj(*('src', 'alinea','adel', 'data')),
                  pj(*('alinea', 'adel', 'echap_leaf_data')): pj(*('src', 'alinea','adel', 'echap_leaf_data'))},

    # Add package platform libraries if any
    zip_safe = False,

    # Scripts
    entry_points = { 'wralea': [ 'adel= alinea.adel',] },
 
    # Dependencies (other are listed in doc to avoid setputools/pip/conda possible conflicts in automatic installs)
    setup_requires = ['openalea.deploy'],
   )



    

#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


import os, sys
from setuptools import setup, find_packages
from openalea.deploy.metainfo import read_metainfo
pj = os.path.join


# Reads the metainfo file
metadata = read_metainfo('metainfo.ini', verbose=True)
for key,value in metadata.iteritems():
    exec("%s = '%s'" % (key, value))


if("win32" in sys.platform):
    install_requires = ['pandas <= 0.7.3']
    setup_requires = install_requires + []
else:
    install_requires = ['pandas <= 0.7.3']
    setup_requires = []
    
pkg_name = 'alinea.adel'
src_rep = 'adel'
packages =[pkg_name]#+['%s.%s'%(pkg_name,x) for x in find_packages(src_rep)]

# Main setup
setup(
    name=name,
    version=version,
    description=description,
    long_description=long_description,
    authors=authors,
    authors_email=authors_email,
    url=url,
    license=license,
    
    namespace_packages = [namespace],
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



    

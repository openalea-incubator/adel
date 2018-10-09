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

    
# #retrieving packages
# pkg_root_dir = '.'
# pkgs = [ pkg for pkg in find_packages(pkg_root_dir) if namespace not in pkg]
# top_pkgs = [pkg for pkg in pkgs if  len(pkg.split('.')) < 2]
# packages = [ namespace + "." + pkg for pkg in pkgs]
# package_dir = dict( [('',pkg_root_dir)] + [(namespace + "." + pkg, pkg_root_dir + "/" + pkg) for pkg in top_pkgs] )
# wralea_entry_points = ['%s = %s'%(pkg,namespace + '.' + pkg) for pkg in top_pkgs]

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
    packages= find_packages('src'),
    # python packages directory
    package_dir= {'': 'src'},

                   
    # Add package platform libraries if any
    include_package_data=True,
    package_data = {'' : ['*.RData', '*.R', '*.8', '*.h', '*.str','*.txt', '*.l', '*.map', '*.csv', '*.png'],},
    share_dirs = {os.path.join(*('alinea', 'adel', 'data')): os.path.join(*('src', 'alinea','adel', 'data'))},

    # Add package platform libraries if any
    zip_safe = False,


    # Scripts
    entry_points = { 'wralea': [ 'adel= alinea.adel',] },
 
    # Dependencies (other are listed in doc to avoid setputools/pip/conda possible conflicts in automatic installs)
    setup_requires = [],
    install_requires = [],
    dependency_links = [],
   )



    

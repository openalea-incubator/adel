# -*- coding: utf-8 -*-
__revision__ = "$Id: $"

import sys
import os

from setuptools import setup, find_packages
from openalea.deploy.metainfo import read_metainfo

# Reads the metainfo file
metadata = read_metainfo('metainfo.ini', verbose=True)
for key,value in metadata.iteritems():
    exec("%s = '%s'" % (key, value))


pkg_root_dir = 'src'
pkgs = [ pkg for pkg in find_packages(pkg_root_dir) if namespace not in pkg]
top_pkgs = [pkg for pkg in pkgs if  len(pkg.split('.')) < 2]
packages = [ namespace + "." + pkg for pkg in pkgs]
package_dir = dict( [('',pkg_root_dir)] + [(namespace + "." + pkg, pkg_root_dir + "/" + pkg) for pkg in top_pkgs] )

setup_requires = ['openalea.deploy']
install_requires = []
dependency_links = ['http://openalea.gforge.inria.fr/pi']

setup(
    name=name,
    version=version,
    description=description,
    long_description=long_description,
    author=authors,
    author_email=authors_email,
    url=url,
    license=license,
    keywords = '',	
    # package installation
    packages= packages,	
    package_dir= package_dir,
    # Namespace packages creation by deploy
    namespace_packages = [namespace],
    create_namespaces = True,

    zip_safe= False,
    # Dependencies
    setup_requires = setup_requires,
    install_requires = install_requires,
    dependency_links = dependency_links,

    include_package_data = True,

    # Declare scripts and wralea as entry_points (extensions) of your package 
    entry_points = { 'wralea': [ 'prospect= alinea.prospect_wralea',]},
    )



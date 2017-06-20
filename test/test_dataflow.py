"""adel dataflow Tests"""

__license__ = "Cecill-C"
__revision__ = " $Id: $"

from openalea.core.alea import run
from openalea.core.pkgmanager import PackageManager
from random import random, randint
""" A unique PackageManager is created for all test of dataflow """
pm = PackageManager()
pm.init(verbose=False)


def test_adelr1():
    """ Test AdelR MonoRun """
    res = run(('alinea.adel.tutorials', 'AdelR MonoRun'),
        inputs={}, pm=pm, vtx_id=13)



# to be repaired
# def test_adelr2():
#     """ Test AdelR MonoRun """
#     res = run(('alinea.adel.tutorials', 'AdelR MonoRun'),
#         inputs={}, pm=pm, vtx_id=39)
#
#
# def test_arvalis():
#     """ Test Leaf db Explorer """
#     res = run(('alinea.adel.tutorials', 'Leaf db Explorer'),
#         inputs={}, pm=pm)
#     assert res[0] == []

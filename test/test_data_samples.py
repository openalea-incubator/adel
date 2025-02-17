from alinea.adel.data_samples import *


def test_leaf_db():
    db = leaves_db()
    sr = srdb()
    xy = xydb()
    leaves = wheat_leaf_db()
    return db, sr, xy, leaves


def test_devT():
    dev = devT()
    return dev


def test_adel():
    two_metamer = adel_two_metamers()
    two_metamer_2 = adel_two_metamers(leaf_sectors=2)
    one_leaf = adel_one_leaf()
    one_sect = adel_one_leaf_element()
    return two_metamer, two_metamer_2, one_leaf, one_sect


def test_stand():
    g, domain_area, domain, convUnit = adel_two_metamers_stand()
    return g, domain_area, domain, convUnit

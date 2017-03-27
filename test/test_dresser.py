from alinea.adel.dresser import blade_dimension, stem_dimension, ear_dimension, \
    dimension_table, AdelDress


def test_dimension():
    blades = blade_dimension()
    stem = stem_dimension()
    ear = ear_dimension()
    dim = dimension_table()


def test_dresser():
    adel = AdelDress()
    g = adel.canopy()
    s = adel.scene(g)
    assert len(s) == 6
    stats = adel.midrib_statistics(g)
    assert all(stats['insertion_height'].values ==(40, 50, 60))


from openalea.mtg import MTG
from alinea.adel.mtg_editions import add_plant, add_axe, add_vegetative_metamer, find_plants, \
    find_metamers, find_label, insert_elements


def test_add_plant():
    g = MTG()
    vid = add_plant(g)
    plant = g.node(vid)
    assert plant.label =='plant1'
    vid = add_plant(g)
    plant = g.node(vid)
    assert plant.label =='plant2'
    vid = add_plant(g, 5)
    plant = g.node(vid)
    assert plant.label =='plant5'
    vid = add_plant(g)
    plant = g.node(vid)
    assert plant.label =='plant6'
    vid = add_plant(g, 1)
    plant = g.node(vid)
    assert plant.label =='plant1'
    assert len(find_plants(g)) == 4

def test_add_metamer():
    g = MTG()
    add_plant(g)
    vid_plant, vid_axe, metamers = find_metamers(g)
    assert len(metamers) == 1
    vid = add_vegetative_metamer(g)
    vid_plant, vid_axe, metamers = find_metamers(g)
    metamer = g.node(vid)
    assert len(metamers) == 2
    assert metamer.label == 'metamer1'

def test_add_axe():
    g = MTG()
    labels = g.property('label')
    p1 = add_plant(g)
    p2 = add_plant(g)
    vid_axe = add_axe(g, 'T0')
    assert g.complex(vid_axe) == p1
    assert labels[g.parent(vid_axe)] == 'MS'
    vid_axe = add_axe(g, 'T1.3', plant_number=2)
    assert g.complex(vid_axe) == p2
    assert labels[g.parent(vid_axe)] == 'T1'
    vid_metamer0 = g.component_roots_at_scale(vid_axe, 3)[0]
    assert labels[g.parent(vid_metamer0)] == "metamer3"


def test_insert_elements():
    g = MTG()
    labels = g.property('label')
    add_plant(g)
    elts = [{'label':'elt1'},{'label':'elt2'}]
    collar = find_label('collar', g)[0]
    insert_elements(g, collar, elts)
    elt1 = find_label('elt1', g)[0]
    assert labels[g.parent(elt1)] == 'baseElement'
    assert labels[g.children(elt1)[0]] == 'elt2'
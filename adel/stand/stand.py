from numpy import linspace, cos, sin
from itertools import *
from math import pi, sqrt
from random import random,sample
from numpy.random import vonmises
from openalea.core.path import path


def regular(nb_plants, nb_rank, dx, dy):
    
    nx = int(nb_plants/nb_rank)
    ny = nb_rank
    domain = ((0,0),(nx*dx, ny*dy))
    return [(i*dx+dx/2., j*dy+dy/2., 0.) for j in xrange(ny) for i in xrange(nx)], domain

def agronomicplot(length, width, sowing_density, plant_density, inter_row, noise = 0,convunit=100,center_scene = True):
    """ Returns the number of plants, the positions, the domain and the domain area of a micro-plot specified with agronomical variables
    length (m) is plot dimension along row direction
    width (m) is plot dimension perpendicular to row direction
    sowing density is the density of seeds sawn
    plant_density is the density of plants that are present (after loss due to bad emergence, early death...)
    inter_row (m) is for the  distance between rows
    noise (%), indicates the precision of the sowing for the inter plant spacing
    convunit is the conversion factor from meter to scene unit
    center_scene allows to center the position arround origin. If False, the scene is in the x+,y+ sector, the origin being at the lower left corner of the domain
    
    Rows are parrallel to x-axis
    Length and Width are adjusted to produce a canopy centered in its domain and compliant with infinitisation
    """
    from operator import itemgetter
    
    inter_plant = 1. / inter_row / sowing_density
    nrow = max(1, int(float(width) / inter_row))
    dx = inter_plant * convunit
    dy = inter_row * convunit
    plant_per_row = max(1,int(float(length) / inter_plant))
    nplants = nrow * plant_per_row
    positions, domain = regular(nplants, nrow, dx, dy)
    n_emerged = int(nplants * plant_density / sowing_density)
    positions = sample(positions, n_emerged)
    # sorting by ranks
    positions = sorted(positions,key= itemgetter(1,0))
    domain_area = abs(domain[1][0] - domain[0][0]) / convunit * abs(domain[1][1] - domain[0][1]) / convunit
    if center_scene:
        xc = float(domain[1][0] + domain[0][0]) / 2
        yc = float(domain[1][1] + domain[0][1]) / 2
        positions = [(x - xc, y - yc, z) for x,y,z in positions]
        domain = ((domain[0][0] - xc,domain[0][1] - yc),(domain[1][0] - xc,domain[1][1] - yc))
    
    return n_emerged, positions, domain, domain_area

def agronomicplotwithdistributions(length, width, sowing_density, plant_density, inter_row, mu=0.0, kappa=3.0, noise = 0,convunit=100):
    """
    Returns the
        - number of plants
        - the positions
        - azimuthal orientation of each plant
        - the domain
        - the simulated density
        of a micro-plot specified with agronomical variables
       
    Inputs:
        - length (m) is plot dimension along row direction
        - width (m) is plot dimension perpendicular to row direction
        - sowing density is the density of seeds sawn
        - plant_density is the density of plants that are present (after loss due to bad emergence, early death...)
        - inter_row (m) is for the  distance between rows
        - mu, kappa are the von Mises parameters of the azimuthal distribution
        - noise (%), indicates the precision of the sowing for the inter plant spacing
        - unit (m or cm) is for the unit of the position and domain
    """
    inter_plant = 1. / inter_row / sowing_density
    nrow = max(1, int(float(width) / inter_row))
    plant_per_row = max(1,int(float(length) / inter_plant))
    nplants = nrow * plant_per_row
    positions, domain = regular(nplants, nrow, inter_plant * convunit, inter_row * convunit)
    n_emerged = int(nplants * plant_density / sowing_density)
    positions = sample(positions, n_emerged)
    density = int(n_emerged / ( abs(domain[1][0] - domain[0][0]) / convunit * abs(domain[1][1] - domain[0][1]) / convunit))
    azimuths = vonmises(mu, kappa, nplants)
    return n_emerged, positions, azimuths, domain, density
    

def regularband(nb_plants, nb_rank, dx, dy):
    
    nx = int(nb_plants/nb_rank)
    ny = nb_rank
    db=sample(xrange(100),nb_plants)
    domain = ((0,0),(nx*dx, ny*dy))
    return [(i*dx+dx/2., j*dy+dy/2.+(db[i+j]-50)/100.*dy/3., 0.) for j in xrange(nb_rank) for i in xrange(nx)], domain

def concentric(nb_plants, distance_plant):
    nb_circle = int(sqrt(nb_plants))
    dp = distance_plant
    points = []
    for i in range(nb_circle):
        n = 2*i+1
        theta = linspace(0,pi,n) 
        radius = i * dp
        x = cos(theta)*radius
        y = sin(theta)*radius
        points.extend([(x[i],y[i],0.) for i in range(len(x)) ])
        
    return points

def uniform(density, nb_plants):
    return []

def planter(scene, distribution):
    """
    Returns a CanestraScene by positionning all elements.
    """
    from openalea.plantgl.all import Translation, Vector3, AxisRotation, Transform4
    from alinea.adel import symbol

    new_scene = symbol.CanScene()
    def transfo(plant,pt):
        r = AxisRotation((0,0,1), random()*pi).getMatrix()
        return plant.transform2(Transform4(Translation(pt).getMatrix()*r))
    l=[new_scene.add_plant(p) for p in imap(transfo, cycle(scene.itervalues()),iter(distribution))]
    return new_scene

def sample_selection(points, gap_fraction):
    """
    Choose a sample from a list of points.
    The gap fraction indicates the proportion of points to remove.
    """
    if gap_fraction == 0:
        return points
    n = len(points)
    k = int(n * (1. - gap_fraction / 100.))
    return sample(points, k)

def clumping_selection(points, nb_clumps, radius_mu, radius_sigma):
    """
    Select a set of points based on andomly choosen groups 
    with a radius following a normal distribution.
    """
    n = len(points)
    centers = sample(range(n), nb_clumps)

    for index in centers:
        radius = normalvariate(radius_mu, radius_sigma)
        r2 = radius**2
        # filter to select all the poit in a ball of radius r2
        # Add these points to the result

def sample_regular_gaps(points, pattern = [0,1]):
    """ Sample points along pattern.
    Pattern is replicated to become as long as points and then used as a filter (0= gap) on points
    Returns positions of plants and gaps 
    """
    
    if not isinstance(points, list):
        points = [points]
    
    p = pattern
    length = len(points)
    p = p * (length / len(p)) + [p[i] for i in range(length % len(p))]
    
    #selection = compress(points,p) only python >= 2.7
    return [point for point,i in izip(points,p) if i],[point for point,i in izip(points,p) if not i]


def post_processing(adel_ouput_path='', plant_number=0, domain_area=0.0, result_path=''):
    import pandas
    assert adel_ouput_path != '' and plant_number != 0 and domain_area != 0.0 \
    and result_path != ''
    
    adel_ouput_path = path(adel_ouput_path)
    adel_ouput_df = pandas.read_csv(adel_ouput_path)
    result_path = path(result_path)
    
    if result_path.exists():
        result_df = pandas.read_csv(result_path)
    else:
        columns = ['Filename', 
                   'aire du plot', 
                   'nombre de plantes', 
                   'temps thermique', 
                   'LAI totale', 
                   'LAI vert', 
                   'PAI total', 
                   'PAI vert', 
                   'nombre axes par m2']
        result_df = pandas.DataFrame(columns=columns)
    
    filename = adel_ouput_path.basename()
    thermal_time = adel_ouput_df['date'][0]
    intermediate_result_1 = domain_area * 10000
    tot_LAI = adel_ouput_df['Slv'].sum() / intermediate_result_1
    green_LAI = (adel_ouput_df['Slv'].sum() - adel_ouput_df['SLsen'].sum()) / intermediate_result_1
    intermediate_result_2 = adel_ouput_df['Slv'].sum() + adel_ouput_df[['SGv', 'SEv']].sum().sum() / 2.0
    tot_PAI = intermediate_result_2 / intermediate_result_1
    intermediate_result_3 = adel_ouput_df['SLsen'].sum() + adel_ouput_df[['SGsen', 'SEsen']].sum().sum() / 2.0
    green_PAI = (intermediate_result_2 - intermediate_result_3) / intermediate_result_1
    axes_density = adel_ouput_df[adel_ouput_df['numphy'] == 1].index.size / domain_area
    
    data = [[filename, domain_area, plant_number, thermal_time, tot_LAI, green_LAI, tot_PAI, green_PAI, axes_density]]
    new_index = [result_df.index.size]
    new_df = pandas.DataFrame(data, index=new_index, columns=result_df.columns)
    
    result_df = pandas.concat([result_df, new_df])
    result_df.to_csv(result_path, na_rep='NA', index=False)
    
    return str(result_path) 

    
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

    print "mariem"
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


def post_processing(adel_output_path='', plant_number=0, domain_area=0.0, postprocessing_results_path=''):
    '''
    Apply post processing on the ADEL output.
    
    For one date, the variables calculated are : TODO: list and explain the 
    variables.
    The results of the post processing are added at the end of a csv file, 
    creating it if necessary.
    Furthermore, for each call to the current function, intermediate post 
    processing results are saved in a new temporary file.
        
    :Parameters:
    
        - `adel_output_path` (str) -  
          path to the csv file describing ADEL output. This file contains data 
          for one date. 
          
        - `plant_number` (int) -  
          the number of plants simulated by ADEL
          
        - `plant_number` (int) -  
          the area of the domain on which ADEL simulation has been done.
          
        - `postprocessing_results_path` (str) -  
          path to the csv file describing the results of the post processing 
          for each date. 
    
    :Returns:
        Two paths: one for the post processing results and another one for the
        intermediate results.
    
    :Returns Type:
        tuple of 2 str

    .. warning:: adel_output_path and postprocessing_results_path must be 
                 non-empty, plant_number and domain_area must be non-null.
    
    '''
    import pandas
    assert adel_output_path != '' and plant_number != 0 and domain_area != 0.0 \
    and postprocessing_results_path != ''
    
    adel_output_path = path(adel_output_path)
    adel_output_df = pandas.read_csv(adel_output_path)
    # Construct the intermediate table
    grouped = adel_output_df.groupby(['plant', 'axe'], as_index=False)
    intermediate_df = pandas.DataFrame(columns=['date', 'refplant_id', 'plant','axe_id', 'axe', 'Slv', 'SLsen', 'SGv', 'SGsen', 'SEv', 'SEsen','SLgreen', 'SGgreen', 'SEgreen', 'NFF', 'HS', 'SSI', 'GreenLeaf'])
    
    date = adel_output_df['date'][0]
    for name, group in grouped:
        refplant_id = group['refplant_id'][group.index[0]]
        plant = name[0]
        axe_id = group['axe_id'][group.index[0]]
        axe = name[1]
        Slv = group['Slv'].sum()
        SLsen = group['SLsen'].sum()
        SGv = group['SGv'].sum()
        SGsen = group['SGsen'].sum()
        SEv = group['SEv'].sum()
        SEsen = group['SEsen'].sum()
        SLgreen = Slv - SLsen
        SGgreen = SGv - SGsen
        SEgreen = SEv - SEsen
        NFF = group['L_shape'][group['L_shape'] != 0.0].count()

	index_of_last_non_null_Ll = group.index[NFF - 1]
        indexes_of_all_non_null_Ll = group.index[:index_of_last_non_null_Ll+1]
        HS = (group['Lv'][indexes_of_all_non_null_Ll] / \
             (group['L_shape'][indexes_of_all_non_null_Ll]).astype(float)) \
	     .sum()
        SSI = (group['Lsen'][indexes_of_all_non_null_Ll] / \
               (group['L_shape'][indexes_of_all_non_null_Ll]).astype(float)) \
              .sum()
        GreenLeaf = HS - SSI
        new_intermediate_data = [[date, refplant_id, plant, axe_id, axe, Slv,SLsen, SGv, SGsen, SEv, SEsen, SLgreen, SGgreen, SEgreen, NFF, HS, SSI, GreenLeaf]]
        new_intermediate_df = pandas.DataFrame(new_intermediate_data, columns=intermediate_df.columns)
        intermediate_df = pandas.concat([intermediate_df, new_intermediate_df], ignore_index=True)
    
    import tempfile
    intermediate_file_suffix = '-%d.csv' % int(date)
    intermediate_results_path = path(tempfile.mktemp(intermediate_file_suffix))
    intermediate_df.to_csv(intermediate_results_path, na_rep='NA', index=False)
    
    # Construct the post processing results
    postprocessing_results_path_ = path(postprocessing_results_path)
    
    if postprocessing_results_path_.exists():
        postprocessing_results_df = pandas.read_csv(postprocessing_results_path_)
    else:
        columns = ['Filename', 'aire du plot', 'Nbr.plant/plot', 'ThermalTime', 
                   'LAI totale', 'LAI vert', 'PAI total', 'PAI vert', 
                   'Nbr.axe.tot/m2', 'Nbr.axe.actif/m2', 'HS.axe0', 'HS.axe1', 
                   'HS.axe2', 'HS.axe3', 'HS.axe4', 'HS.axe5', 'SSI.axe0',    
                   'SSI.axe1', 'SSI.axe2', 'SSI.axe3', 'SSI.axe4', 'SSI.axe5']
        postprocessing_results_df = pandas.DataFrame(columns=columns)
    
    Filename = intermediate_results_path.basename()
    area_in_cm = domain_area * 10000.0
    tot_LAI = intermediate_df['Slv'].sum() / area_in_cm 
    green_LAI = intermediate_df['SLgreen'].sum() / area_in_cm 
    tot_PAI = (intermediate_df['Slv'].sum() \
               + intermediate_df[['SGv', 'SEv']].sum().sum() / 2.0) \
               / area_in_cm
    green_PAI = (intermediate_df['SLgreen'].sum() \
                 + intermediate_df[['SGgreen', 'SEgreen']].sum().sum() / 2.0) \
                 / area_in_cm
    axes_density = intermediate_df.index.size / float(domain_area)                
    active_axes_density_df = intermediate_df[intermediate_df['HS'] > 0.5]
    active_axes_density_df = active_axes_density_df[active_axes_density_df['HS'] < active_axes_density_df['NFF']]
    active_axes_density_df = active_axes_density_df[active_axes_density_df['SLgreen'] > 0.0]
    active_axes_density = active_axes_density_df.index.size / float(domain_area)
    
    new_postprocessing_results_data = [[Filename, domain_area, plant_number, date, tot_LAI, green_LAI, tot_PAI, green_PAI, axes_density, active_axes_density]]
    
    for column_name in ['HS', 'SSI']: 
        for axe_idx in range(6):
            current_mean = intermediate_df[intermediate_df['axe'] == axe_idx][column_name].mean()
            new_postprocessing_results_data[0].append(current_mean)    
        
    new_postprocessing_results_df = pandas.DataFrame(new_postprocessing_results_data,columns=postprocessing_results_df.columns)
    postprocessing_results_df = pandas.concat([postprocessing_results_df, new_postprocessing_results_df], ignore_index=True)
    postprocessing_results_df.to_csv(postprocessing_results_path_,na_rep='NA', index=False)
    
    return (postprocessing_results_path, str(intermediate_results_path))

    

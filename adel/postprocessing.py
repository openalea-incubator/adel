# -*- python -*-
#
#       adel.povray
#
#       Copyright 2006-2012 INRIA - CIRAD - INRA
#
#       File author(s): Camille Chambon <camille.chambon@grignon.inra.fr>
#                       Christian.Fournier <christian.fournier@supagro.inra.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################
import numpy as np
import pandas as pd
import openalea.plantgl.all as pgl

#analysis of Ground cover

def domain3D(domain2D, scene):
    t=pgl.Tesselator()
    bbc = pgl.BBoxComputer(t)
    bbc.process(scene)
    bbox = bbc.result

    z_base = bbox.getZMin()
    z_top = bbox.getZMax()
    domain3D = (domain2D[0] + (z_base,), domain2D[1] + (z_top,))
    return domain3D
        
        
def stand_box(domain):
    '''
    
    domain: 3D bounding box of the stand
    '''
    # list of points
    z_base = domain[0][2]
    z_top = domain[1][2]
    sides_points = [domain[0], # coordinates of bottom right corner
                  (domain[1][0], domain[0][1], z_base),
                  (domain[1][0], domain[1][1], z_base), # coordinates of bottom left corner
                  (domain[0][0], domain[1][1], z_base),
                  (domain[0][0], domain[0][1], z_top), # coordinates of top right corner
                  (domain[1][0], domain[0][1], z_top),
                  (domain[1][0], domain[1][1], z_top),    # coordinates of top left corner
                  (domain[0][0], domain[1][1], z_top)]
                  
    bottom_points = [domain[0], # coordinates of bottom right corner
                  (domain[1][0], domain[0][1], z_base),
                  (domain[1][0], domain[1][1], z_base), # coordinates of bottom left corner
                  (domain[0][0], domain[1][1], z_base)]

    # list of indices to make the quads of the sides from the points
    side_indices = [(0, 1, 5, 4), #
               (1, 2, 6, 5), # indices for 
               (2, 3, 7, 6), # side faces
               (3, 0, 4, 7)] #         
     
    # list of indices to make the quads of the bottom from the points
    bottom_indices = [(0, 1, 2, 3)] # indices for bottom face

    # list of colors
    side_color = pgl.Color3(0, 0, 0)
    bottom_color = pgl.Color3(255, 255, 255)

    # construction of the geometry for the sides
    side_box = pgl.QuadSet(sides_points, side_indices)
    # construction of the geometry for the bottom
    bottom_box = pgl.QuadSet(bottom_points, bottom_indices)
                           
    # create 2 shapes: 1 with side_color, 1 with bottom_color 
    sides_shape = pgl.Shape(side_box, pgl.Material(side_color))
    bottom_shape = pgl.Shape(bottom_box, pgl.Material(bottom_color))
    
    scene = pgl.Scene()
    scene.add(sides_shape)
    scene.add(bottom_shape)
    
    return scene

def color_count(image,RGBcolor = (0,255,0)):
    import cv2
    BGRcolor = np.array(RGBcolor[::-1])
    res = cv2.inRange(image, BGRcolor, BGRcolor)
    return res.sum() / 255.
    
def replicate_scene(scene, domain, replicate=None):
    """ replicate the scene around the domain
    """
    if replicate is None:
        newscene = scene
        newdomain = domain
    else:
        dx = abs(domain[0][0] - domain[1][0])
        dy = abs(domain[0][1] - domain[1][1])
        
        newscene = pgl.Scene()
        for tx in (-dx,0,dx):
            for ty in (-dy,0,dy): 
                rep = scene.deepcopy()
                for sh in rep:
                    sh.geometry = pgl.Translated(tx,ty,0,sh.geometry)
                newscene += rep
        newdomain = ( (domain[0][0] - dx, domain[0][1] - dy),(domain[1][0] + dx, domain[1][1] + dy) )
        if replicate == 1:
            replicate = None
        else:
            replicate -= 1
        newscene, newdomain = replicate_scene(newscene, newdomain, replicate)
    return newscene, newdomain
    
def color_ground_cover(scene, domain, colors_def = {'green':[0, 255, 0], 'senescent' : [255, 0, 0]}, camera = {'type':'perspective', 'distance':200., 'fov':50., 'azimuth':0, 'zenith':0.}, image_width = 4288, image_height = 2848, getImages=False, replicate=None):
    """
    Compute ground_cover fraction over domain for each color declared in colors
    
    depends on cv2 and povray
    
    """
    import cv2
    from alinea.adel.povray.povray import PovRay
    
    newscene, domain = replicate_scene(scene, domain, replicate=replicate)
    
    camera['xc'] = 0.5 * (domain[0][0] + domain[1][0])
    camera['yc'] = 0.5 * (domain[0][1] + domain[1][1])
    pov = PovRay(camera=camera, image_width=image_width, image_height=image_height)
    pov.render(newscene)
    im = pov.get_image(cv2.imread)
    
    d3D = domain3D(domain, newscene)
    scene_box = pgl.Scene()
    scene_box.add(stand_box(d3D))
    pov.render(scene_box)
    box = pov.get_image(cv2.imread)
    
    mask = np.uint8(box[:,:,0] / box.max() * 255)
    total_domain = mask.sum() / 255.
    masked = cv2.bitwise_and(im,im,mask=mask)
    
    res = {k:color_count(masked,v) / total_domain for k,v in colors_def.iteritems()}
    
    if not getImages:
        im = None
        box = None
    return res, im, box
        
def ground_cover(g, domain, camera = {'type':'perspective', 'distance':200., 'fov':50., 'azimuth':0, 'zenith':0.}, image_width = 4288, image_height = 2848, getImages=False, replicate=None):
    """
    compute ground cover based on is_green/not_green property of g
    """
    from alinea.adel.mtg_interpreter import plot3d
    
    colors_def = {'green':[0, 255, 0], 'senescent' : [255, 0, 0]}
    greeness = g.property('is_green')
    colors = {k:colors_def['green'] if greeness[k] else colors_def['senescent'] for k in greeness}
    
    scene = plot3d(g, colors=colors)
    
    gc, im, box = color_ground_cover(scene, domain, colors_def = colors_def, camera = camera, image_width = image_width, image_height = image_height, getImages=getImages, replicate=replicate)
    
    return gc, im, box


# Analysis of axis and LAI dynamics



def aggregate_adel_output(adel_output_df, by=['TT', 'plant', 'axe_id']):
    '''
    Aggregate the output of Adel.
    
    :Parameters:
    
        - `adel_output_df` (pandas.DataFrame) - the output of Adel (return by method adel.get_exposed_areas(g)).
        
        - `by` (list) - the name of the columns to group adel_output_df on.
    
    :Returns:
        A table with the following columns: 'TT', 'plant', 'axe_id', 'refplant_id', 
        'Slv', 'Slvsen', 'SGv', 'SGvsen', 'SEv', 'SEvsen', 'Slvgreen', 'SGvgreen', 
        and 'SEvgreen'.
        
    :Returns Type:
        pandas.DataFrame

    '''

    aggregations_to_apply_list = [('TT', lambda x: x[x.first_valid_index()]),
                                  ('plant', lambda x: x[x.first_valid_index()]),
                                  ('axe_id', lambda x: x[x.first_valid_index()]),
                                  ('refplant_id', lambda x: x[x.first_valid_index()]),
                                  ('Slv', 'sum'),
                                  ('Slvsen', 'sum'),
                                  ('SGv', 'sum'),
                                  ('SGvsen', 'sum'),
                                  ('SEv', 'sum'),
                                  ('SEvsen', 'sum'),
                                  ('Slvgreen', 'sum'),
                                  ('SGvgreen', 'sum'),
                                  ('SEvgreen', 'sum')]
                             
    column_names = np.array(aggregations_to_apply_list)[:, 0]
    aggregations_to_apply_dict = dict(aggregations_to_apply_list)
    aggregations_to_apply_dict = dict([(key, value) for key, value in aggregations_to_apply_dict.iteritems() if key not in by])
                             
    grouped = adel_output_df.groupby(by, as_index=False)
    adel_output_aggregated_df = grouped.aggregate(aggregations_to_apply_dict)
    
    return adel_output_aggregated_df.reindex_axis(column_names, axis=1)
            
    
def phenology(adel_output_df):
    '''
    Calculate the phenology from the output of Adel.
    
    :Parameters:
    
        - `adel_output_df` (pandas.DataFrame) - the output of Adel.
    
    :Returns:
        A table with the following columns: 'TT', 'plant', 'axe_id', 'NFF', 'HS', 
        'SSI', 'GreenLeaf', 'NFL', 'NFV', 'has_ear', 'd_base-lastcol'. 
        
    :Returns Type:
        pandas.DataFrame
 
    '''
    
    phenology_df = pd.DataFrame(columns=['TT', 'plant', 'axe_id', 'NFF',
                                         'HS', 'SSI', 'GreenLeaf', 'NFL',
                                         'NFV', 'has_ear', 'd_base-lastcol', 
                                         'HS_final'])
    
    for name, group in adel_output_df.groupby(['TT', 'plant', 'axe_id'], as_index=False):
        TT, plant, axe_id = name
        L_shape = group['L_shape']
        indexes_of_vegetative_phytomers = group[group['numphy'] <= group['nff']].index
        NFF = indexes_of_vegetative_phytomers.size
        NFF = group['nff'].values[0]
        # HS
        indexes_of_all_non_null_Lv = group['Lv'][group['Lv'] != 0].index
        NFV = len(indexes_of_all_non_null_Lv)
        Lv_non_null_series = group['Lv'][indexes_of_vegetative_phytomers][indexes_of_all_non_null_Lv]
        if len(Lv_non_null_series) == 0:
            HS = 0.0
            NFL = 0.0
            d_base_lastcol = np.nan
        else:
            Lv_equal_L_shape_series = Lv_non_null_series[Lv_non_null_series == L_shape[indexes_of_all_non_null_Lv]]
            if Lv_equal_L_shape_series.index.size == 0:
                HS_indexes = Lv_non_null_series.index
                NFL = 0.0
                d_base_lastcol = np.nan
            else:
                index_of_last_Lv_equal_L_shape = Lv_equal_L_shape_series.index[-1]
                index_of_first_Lv_not_equal_L_shape = index_of_last_Lv_equal_L_shape + 1
                index_of_last_Lv_non_null = Lv_non_null_series.index[-1]
                HS_indexes = range(index_of_first_Lv_not_equal_L_shape, index_of_last_Lv_non_null + 1)
                NFL = group['numphy'][index_of_last_Lv_equal_L_shape]
                d_base_lastcol = group['d_basecol'][index_of_last_Lv_equal_L_shape]
            HS = NFL + \
                 (Lv_non_null_series[HS_indexes] / \
                  L_shape[HS_indexes].astype(float)) \
                 .sum()
        indexes_of_all_non_null_Lsen = group['Lvsen'][group['Lvsen'] != 0].index
        Lsen_non_null_series = group['Lvsen'][indexes_of_vegetative_phytomers][indexes_of_all_non_null_Lsen]
        Lsen_equal_L_shape_series = Lsen_non_null_series[Lsen_non_null_series == L_shape[indexes_of_all_non_null_Lsen]]
        # SSI
        SSI = pd.Series(0.0, index=group.index)
        if indexes_of_all_non_null_Lsen.size != 0:
            nonzero_Slvsen_indexes = group.index[group.Slvsen.nonzero()]
            last_nonzero_Slvsen_index = nonzero_Slvsen_indexes[-1]
            for nonzero_Slvsen_index in nonzero_Slvsen_indexes:
                n = group['numphy'][nonzero_Slvsen_index]
                Slvsen_values = group['Slvsen'].loc[nonzero_Slvsen_index:last_nonzero_Slvsen_index]
                S_shape_values = group['S_shape'].loc[nonzero_Slvsen_index:last_nonzero_Slvsen_index]
                SSI[nonzero_Slvsen_index] = n + (Slvsen_values / S_shape_values).sum()
                
        indexes_of_all_null_Lshape = L_shape[L_shape == 0.0].index
        HS_final = group['HS_final'][group.first_valid_index()]
        has_ear = int(group[group['HS_final'] >= group['nff']].any().any())
              
        GreenLeaf = HS - SSI
       
        new_phenology_df = pd.DataFrame(index=group.index, 
                                        columns=phenology_df.columns)
        
        new_phenology_df['TT'] = TT
        new_phenology_df['plant'] = plant
        new_phenology_df['axe_id'] = axe_id
        new_phenology_df['NFF'] = NFF
        new_phenology_df['HS'] = HS
        new_phenology_df['SSI'] = SSI
        new_phenology_df['GreenLeaf'] = GreenLeaf
        new_phenology_df['NFL'] = NFL
        new_phenology_df['NFV'] = NFV
        new_phenology_df['has_ear'] = has_ear
        new_phenology_df['d_base-lastcol'] = d_base_lastcol
        new_phenology_df['HS_final'] = HS_final
        
        phenology_df = pd.concat([phenology_df, new_phenology_df], 
                                 ignore_index=True)
    
    return phenology_df
                                                

def axis_statistics(adel_output_df, domain_area, convUnit=0.01):
    '''
    Calculate statistics on the axes from adel output. 
    Calls aggregate_adel_output and phenology.
    
    :Parameters:
    
        - `adel_output_df` (pandas.DataFrame) - the output of Adel.
        
        - `domain_area` (float) -  ground surface area occupied by the plants.
        
        - `convUnit` (float) -  
          Factor to convert the length unit of adel output to meter. Default '0.01'.
    
    :Returns:
        A table with the following columns: 'ThermalTime', 'axe_id', 'NFF', 'HS', 
        'SSI', 'LAI totale', 'LAI vert', 'PAI total', 'PAI vert', 'has_ear', 
        'd_base-lastcol', 'axes_cardinality', 'growing_axes_cardinality', 
        'active_axes_cardinality'
        
    :Returns Type:
        pandas.DataFrame
 
    '''
    aggregated_adel_output_df = aggregate_adel_output(adel_output_df)
    phenology_df = phenology(adel_output_df)
    intermediate_df = pd.merge(aggregated_adel_output_df, phenology_df, on=['TT', 'plant', 'axe_id'])
    
    area_in_cm = domain_area * 1.0 / convUnit ** 2

    axis_statistics_df = \
        pd.DataFrame(columns=['ThermalTime', 'axe_id', 'NFF', 'HS', 'SSI', 'LAI totale', 'LAI vert', 
                                  'PAI total', 'PAI vert', 'has_ear', 'd_base-lastcol', 'axes_cardinality', 
                                  'growing_axes_cardinality', 'active_axes_cardinality', 'axis_order'])
    for (ThermalTime, axe_id, NFF, has_ear), group in intermediate_df.groupby(['TT', 'axe_id', 'NFF', 'has_ear'], as_index=False):
        HS = group['HS'].mean()
        SSI = group['SSI'].mean()
        tot_LAI = group['Slv'].sum() / area_in_cm
        green_LAI = group['Slvgreen'].sum() / area_in_cm
        tot_PAI = (group['Slv'] + (group['SGv'] + group['SEv']) / 2.0 ).sum() / area_in_cm
        green_PAI = (group['Slvgreen'] + (group['SGvgreen'] + group['SEvgreen']) / 2.0 ).sum() / area_in_cm
        d_base_lastcol = group['d_base-lastcol'].mean()
        axes_cardinality = len(group)
        growing_axes_cardinality_df = group[(group['HS'] < group['HS_final']) & (group['HS'] > 0.5) & (group['Slv'] > 0)]
        #growing_axes_cardinality_df = growing_axes_cardinality_df[growing_axes_cardinality_df['HS'] > 0.5]
        growing_axes_cardinality = len(growing_axes_cardinality_df)
        
        # Calculate the number of active axes. An axis is active if:
        # - any of its metamers is growing,
        # - XOR one of its phytomers is a ear (i.e. HS_final == NFF)
        growing_indexes = growing_axes_cardinality_df.index
        not_growing_indexes = group.index - growing_indexes
        is_ear_indexes = not_growing_indexes.intersection(group[group['HS_final'] == group['NFF']].index)
        active_axes_indexes = growing_indexes + is_ear_indexes
        active_axes_cardinality = len(active_axes_indexes)
        if axe_id == 'MS':
            axis_order = 0
        else:
            axis_order = axe_id.count('.') + 1
        
        new_axis_statistics_data = [[ThermalTime, axe_id, NFF, HS, SSI, tot_LAI, 
                                     green_LAI, tot_PAI, green_PAI, has_ear, 
                                     d_base_lastcol, axes_cardinality, 
                                     growing_axes_cardinality, 
                                     active_axes_cardinality, axis_order]]
        new_axis_statistics_df = pd.DataFrame(new_axis_statistics_data, 
                                                  columns=axis_statistics_df.columns)
        axis_statistics_df = pd.concat([axis_statistics_df, 
                                            new_axis_statistics_df], 
                                           ignore_index=True)
    
    return axis_statistics_df
    
    
def plot_statistics(axis_statistics_df, plant_number, domain_area):
    '''
    Calculate the statistics on the plot from adel output.
    
    :Parameters:
    
        - `axis_statistics_df` (pandas.DataFrame) - the statistics on the axes.
        
        - `plant_number` (float) - the number of plants.
        
        - `domain_area` (float) -  ground surface area occupied by the plants.
    
    :Returns:
        A table with the following columns: 'aire du plot', 'Nbr.plant.perplot', 
        'ThermalTime', 'LAI_tot', 'LAI_vert', 'PAI_tot', 'PAI_vert', 'Nbr.axe.tot.m2', 
        'Nbr.axe.actif.m2.old', 'number_of_active_axes_per_m2'
        
    :Returns Type:
        pandas.DataFrame
    '''
    plot_statistics_df = \
        pd.DataFrame(columns=['aire du plot', 'Nbr.plant.perplot', 'ThermalTime', 
                              'LAI_tot', 'LAI_vert', 'PAI_tot', 'PAI_vert', 
                              'Nbr.axe.tot.m2', 'Nbr.axe.actif.m2.old', 'number_of_active_axes_per_m2'])
                                  
    for ThermalTime, group in axis_statistics_df.groupby('ThermalTime', as_index=False):
                                  
        tot_LAI = group['LAI totale'].sum() 
        green_LAI = group['LAI vert'].sum()
        tot_PAI = group['PAI total'].sum()
        green_PAI = group['PAI vert'].sum()
        axes_density = group['axes_cardinality'].sum() / float(domain_area) 
        growing_axes_density = group['growing_axes_cardinality'].sum() / float(domain_area)
        active_axes_density = group['active_axes_cardinality'].sum() / float(domain_area)
                                      
        new_plot_statistics_data = [[domain_area, plant_number, ThermalTime, 
                                     tot_LAI, green_LAI, tot_PAI, green_PAI, 
                                     axes_density, growing_axes_density, 
                                     active_axes_density]]
        new_plot_statistics_df = \
        pd.DataFrame(new_plot_statistics_data, 
                         columns=plot_statistics_df.columns)
        plot_statistics_df = \
            pd.concat([plot_statistics_df, new_plot_statistics_df], 
                          ignore_index=True)
                          
    return plot_statistics_df
    
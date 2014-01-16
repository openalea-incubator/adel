import numpy as np
import pandas as pd

def aggregate_adel_output(adel_output_df, by=['TT', 'plant', 'axe_id']):
    '''
    Aggregate the output of Adel.
    
    :Parameters:
    
        - `adel_output_df` (pandas.DataFrame) - the output of Adel.
        
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
        if indexes_of_all_non_null_Lsen.size == 0:
            SSI = 0.0
        elif Lsen_equal_L_shape_series.index.size == 0:
            SSI = (Lsen_non_null_series[indexes_of_all_non_null_Lsen] / \
                   L_shape[indexes_of_all_non_null_Lsen].astype(float)) \
                  .sum()
        else:
            index_of_last_Lsen_equal_L_shape = Lsen_equal_L_shape_series.index[-1]
            NFS = group['numphy'][index_of_last_Lsen_equal_L_shape]
            index_of_first_Lsen_not_equal_L_shape = index_of_last_Lsen_equal_L_shape + 1
            index_of_last_Lsen_non_null = indexes_of_all_non_null_Lsen[-1]
            SSI_indexes = range(index_of_first_Lsen_not_equal_L_shape, index_of_last_Lsen_non_null + 1)
            SSI = NFS + \
                 (Lsen_non_null_series[SSI_indexes] / \
                  L_shape[SSI_indexes].astype(float)) \
                 .sum()
        indexes_of_all_null_Lshape = L_shape[L_shape == 0.0].index
        HS_final = group['HS_final'][group.first_valid_index()]
        has_ear = int(group[group['HS_final'] >= group['nff']].any().any())
              
        GreenLeaf = HS - SSI
        new_phenology_data = [[TT, plant, axe_id, NFF, HS, SSI, GreenLeaf, NFL,
                              NFV, has_ear, d_base_lastcol, HS_final]]
       
        new_phenology_df = pd.DataFrame(new_phenology_data, 
                                        columns=phenology_df.columns)
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
                                  'growing_axes_cardinality', 'active_axes_cardinality'])
    for (ThermalTime, axe_id, NFF, has_ear), group in intermediate_df.groupby(['TT', 'axe_id', 'NFF', 'has_ear'], as_index=False):
        HS = group['HS'].mean()
        SSI = group['SSI'].mean()
        tot_LAI = group['Slv'].sum() / area_in_cm
        green_LAI = group['Slvgreen'].sum() / area_in_cm
        tot_PAI = (group['Slv'] + (group['SGv'] + group['SEv']) / 2.0 ).sum() / area_in_cm
        green_PAI = (group['Slvgreen'] + (group['SGvgreen'] + group['SEvgreen']) / 2.0 ).sum() / area_in_cm
        d_base_lastcol = group['d_base-lastcol'].mean()
        axes_cardinality = len(group)
        growing_axes_cardinality_df = group[group['HS'] < group['HS_final']]
        growing_axes_cardinality_df = growing_axes_cardinality_df[growing_axes_cardinality_df['HS'] > 0.5]
        growing_axes_cardinality = len(growing_axes_cardinality_df)
        
        # Calculate the number of active axes. An axis is active if:
        # - any of its metamers is growing,
        # - XOR one of its phytomers is a ear (i.e. HS_final == NFF)
        growing_indexes = growing_axes_cardinality_df.index
        not_growing_indexes = group.index - growing_indexes
        is_ear_indexes = not_growing_indexes.intersection(group[group['HS_final'] == group['NFF']].index)
        active_axes_indexes = growing_indexes + is_ear_indexes
        active_axes_cardinality = len(active_axes_indexes)
        
        new_axis_statistics_data = [[ThermalTime, axe_id, NFF, HS, SSI, tot_LAI, 
                                     green_LAI, tot_PAI, green_PAI, has_ear, 
                                     d_base_lastcol, axes_cardinality, 
                                     growing_axes_cardinality, active_axes_cardinality]]
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
                                  
        tot_LAI = axis_statistics_df['LAI totale'].sum() 
        green_LAI = axis_statistics_df['LAI vert'].sum()
        tot_PAI = axis_statistics_df['PAI total'].sum()
        green_PAI = axis_statistics_df['PAI vert'].sum()
        axes_density = axis_statistics_df['axes_cardinality'].sum() / float(domain_area) 
        growing_axes_density = axis_statistics_df['growing_axes_cardinality'].sum() / float(domain_area)
        active_axes_density = axis_statistics_df['active_axes_cardinality'].sum() / float(domain_area)
                                      
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
    
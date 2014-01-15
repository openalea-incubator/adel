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
                                         'NFV', 'has_ear', 'd_base-lastcol'])
    
    for name, group in adel_output_df.groupby(['TT', 'plant', 'axe_id'], as_index=False):
        TT, plant, axe_id = name
        L_shape = group['L_shape']
        # Utiliser nff (dispo dans adel output)  et remplacer 'indexes_of_all_non_null_Lshape' par 'indexes_of_vegetative_phytomers' = phytomer avec numphy <= nff
        indexes_of_all_non_null_Lshape = L_shape[L_shape != 0.0].index 
        NFF = indexes_of_all_non_null_Lshape.size
        # HS
        indexes_of_all_non_null_Lv = group['Lv'][group['Lv'] != 0].index
        NFV = len(indexes_of_all_non_null_Lv)
        Lv_non_null_series = group['Lv'][indexes_of_all_non_null_Lshape][indexes_of_all_non_null_Lv]
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
        Lsen_non_null_series = group['Lvsen'][indexes_of_all_non_null_Lshape][indexes_of_all_non_null_Lsen]
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
        # 15/01 (christian): a mon avis enlever 'has_ear' et garder juste  la colone 'HS_final', disponible dans nouveau adel output . 
        # si tu veux garder la definition est :
        # has_ear = group[group['HS_final'] >= group['NFF']] ; rename has_ear to survivors
        has_ear = int(group[['El','Ed']].ix[indexes_of_all_null_Lshape].any().any())
              
        GreenLeaf = HS - SSI
        new_phenology_data = [[TT, plant, axe_id, NFF, HS, SSI, GreenLeaf, NFL,
                              NFV, has_ear, d_base_lastcol]]
       
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
        # TODO: use HS and HS_final
        # growing_axes_cardinality_df = group['HS'] > 0.5 and (group['HS'] < group['HS_final'])
        growing_axes_cardinality_df = group[group['HS'] > 0.5]
        growing_axes_cardinality_df = \
        growing_axes_cardinality_df[growing_axes_cardinality_df['HS'] <= growing_axes_cardinality_df['NFF']]
        growing_axes_cardinality_df = growing_axes_cardinality_df[growing_axes_cardinality_df['Slvgreen'] > 0.0]
        growing_axes_cardinality = len(growing_axes_cardinality_df)
        # calculate the number of active axes without ear
        # TODO: use HS_final and nff to make selection :
        # active = axes with nff <= HS_final OR (axes with nff > HSfinal AND HS < HS_final)
        # - (any of its metamers) is growing
        # - OR one of its is a ear.
        # active_axes = growing_axes_cardinality_df or (not growing_axes_cardinality_df and group['HS_final'] == group['NFF'])
        HS_positive = group[group['HS'] > 0]
        has_ear_0 = HS_positive[HS_positive['has_ear'] == 0]
        active_axes_without_ear_df = has_ear_0[has_ear_0['HS'] < has_ear_0['NFF']]
        number_of_active_axes_without_ear = len(active_axes_without_ear_df)
        # calculate the number of active axes with ear
        active_axes_with_ear_df = HS_positive[HS_positive['has_ear'] == 1]
        number_of_active_axes_with_ear = len(active_axes_with_ear_df)
        # calculate of the number of active axes per square meter
        active_axes_cardinality = number_of_active_axes_with_ear + number_of_active_axes_without_ear
        
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
    
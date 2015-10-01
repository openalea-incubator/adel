def AdelRunOptions(leafDuration, fracLeaf, stemDuration, dHS_col, dHS_en, senescence_leaf_shrink, epsillon,HSstart_inclination_tiller, rate_inclination_tiller, drop_empty):
    '''    Node interface for Adel General Option dict
    '''
    adel_general_option_dict = {'leafDuration':leafDuration,
		    'fracLeaf' : fracLeaf,
		    'stemDuration' : stemDuration,
            'dHS_col': dHS_col,
            'dHS_en':dHS_en,
		    'senescence_leaf_shrink' : senescence_leaf_shrink,
		    'epsillon' : 1e-6,
            'HSstart_inclination_tiller': HSstart_inclination_tiller,
            'rate_inclination_tiller': rate_inclination_tiller,
            'drop_empty':drop_empty} 
    return adel_general_option_dict,

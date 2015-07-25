def AdelRunOptions(startleaf, endleaf, stemleaf, senescence_leaf_shrink, epsillon,HSstart_inclination_tiller, rate_inclination_tiller, drop_empty):
    '''    Node interface for Adel General Option dict
    '''
    adel_general_option_dict = {'startLeaf':startleaf,
		    'endLeaf' : endleaf,
            'endLeaf1' : endleaf,
		    'stemLeaf' : stemleaf,
		    'senescence_leaf_shrink' : senescence_leaf_shrink,
		    'epsillon' : 1e-6,
            'HSstart_inclination_tiller': HSstart_inclination_tiller,
            'rate_inclination_tiller': rate_inclination_tiller,
            'drop_empty':drop_empty} 
    return adel_general_option_dict,

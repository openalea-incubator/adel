def AdelRunOptions(startleaf, endleaf, stemleaf, senescence_leaf_shrink, epsillon):
    '''    Node interface for Adel General Option dict
    '''
    adel_general_option_dict = {'startLeaf':startleaf,
		    'endLeaf' : endleaf,
            'endLeaf1' : endleaf,
		    'stemLeaf' : stemleaf,
		    'senescence_leaf_shrink' : senescence_leaf_shrink,
		    'epsillon' : 1e-6} 
    return adel_general_option_dict,

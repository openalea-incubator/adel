def vcaribu_adaptor(caribuscene, direct, nz, dz, ds):
    '''    Adaptor to make vcaribu compatible with caribuScene and Scattering options nodes
    '''
    scene = None
    lightsources = None
    pattern = None
    if caribuscene.hasScene:	
        scene = caribuscene.scene
    if caribuscene.hasSources:
        lightsources = caribuscene.sources
    opticals = caribuscene.PO
    if caribuscene.hasPattern:
        pattern = caribuscene.pattern
    optiondict = {'1st':direct,'Nz':nz,'Hc':dz,'Ds':ds,'wavelength':caribuscene.wavelength}; 
    # write the node code here.

    # return outputs
    return scene, lightsources, opticals, pattern, optiondict

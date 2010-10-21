def vcaribu(canopy, lightsource, optics, pattern, options):
    """
    Visualea interface to Caribu class
    Caribu allows nested radiosity illumination on a 3D scene.
    """

    from alinea.caribu.caribu  import Caribu
    
    sim = Caribu(resdir = None, resfile = None)#no output on disk
    # --canfile 
    sim.scene = canopy
    # --optics 
    sim.opticals = optics
    #--skyfile 
    sim.sky = lightsource               
    #--pattern 
    sim.pattern = pattern
    #--options (if different from caribu defaults)
    if options is not None:
        #--scatter
        if '1st' in options.keys():
            sim.direct= options['1st']
        #--nb_layers
        if 'Nz' in options.keys():
            sim.nb_layers =  options['Nz'] 
        #--can_height
        if 'Hc' in options.keys():
            sim.can_height =  options['Hc'] 
        #--sphere_diameter
        if 'Ds' in options.keys():
            sim.sphere_diameter =  options['Ds']
        #--debug mode (if True, prevent removal of tempdir)
        if 'debug' in options.keys():
            sim.my_dbg = options['debug']
        #--names of optical properties (usefull if opticals are given as strings
        if 'wavelength' in options.keys():
            sim.optnames = options['wavelength']
    status = str(sim)
    sim.run()
    irradiances=sim.nrj

    # return outputs
    return irradiances,status

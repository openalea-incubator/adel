from alinea.nema import *

def nemafile(parameters, drivingvariables, state0dd, outdir,steps,substeps):
    '''    run nema with files inputs
    '''
 
    # write the node code here.
    DParamPlant=parameters+"/ParamPlant.txt";
    DParamEnv=parameters+"/ParamEnv.txt";
    DParamGrain=parameters+"/ParamGrain.txt";
    DParamLamina=parameters+"/ParamLamina.txt";
    DParamSheath=parameters+"/ParamSheath.txt";
    DParamInternode=parameters+"/ParamInternode.txt";
    DParamPeduncle=parameters+"/ParamPeduncle.txt";
    DParamChaff=parameters+"/ParamChaff.txt";
    DParamRoot=parameters+"/ParamRoot.txt";

    #path of DrivingVariables 

    # H0, H1, H2 refers to treatments H0, H3 and H15 (paper definition) 
    DNsoilH0=drivingvariables+"/NsoilH2.txt"; 
    DTimeSoil=drivingvariables+"/TimeSoil.txt";
    Dmeteo=drivingvariables+"/meteo(old2).txt";

    #path of State0dd
    DState0Grain=state0dd+"/State0Grain.txt";
    DState0Lamina=state0dd+"/State0Lamina.txt";
    DState0Sheath=state0dd+"/State0Sheath.txt";
    DState0Internode=state0dd+"/State0Internode.txt";
    DState0Peduncle=state0dd+"/State0Peduncle.txt";
    DState0Chaff=state0dd+"/State0Chaff.txt";

    grain = Grain(DParamGrain, DState0Grain, outdir)
    root = Root(DParamRoot, outdir)
    environment = Environment(DParamEnv,DNsoilH0,DTimeSoil,Dmeteo)
    lamina = Entity(DParamLamina,DState0Lamina,4, outdir)
    sheath = Entity(DParamSheath,DState0Sheath,4, outdir)
    internode = Entity(DParamInternode,DState0Internode,4, outdir)
    peduncle = Entity(DParamPeduncle,DState0Peduncle,1, outdir)
    chaff = Entity(DParamChaff,DState0Chaff,1, outdir)

    plant = Plant(grain, root, environment, lamina, sheath, internode, peduncle, chaff)

    nSTEPS = steps
    step = 1
    nbsubsteps = substeps

    plant.InitialPlant(nbsubsteps)

    for i in range(nSTEPS):
        plant.ActuPlant(step, nbsubsteps)
    print 'Simulation %d'%(i+1)

    plant.Display(outdir, nbsubsteps)
    print 'Nema: Simulation succeed'

    # return outputs
    return outdir

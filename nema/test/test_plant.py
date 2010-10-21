from alinea.nema import *
from numpy import *
from pylab import *

inputdir = ''
outputdir = ''

DParamPlant=inputdir+"parameters/ParamPlant.txt";
DParamEnv=inputdir+"parameters/ParamEnv.txt";
DParamGrain=inputdir+"parameters/ParamGrain.txt";
DParamLamina=inputdir+"parameters/ParamLamina.txt";
DParamSheath=inputdir+"parameters/ParamSheath.txt";
DParamInternode=inputdir+"parameters/ParamInternode.txt";
DParamPeduncle=inputdir+"parameters/ParamPeduncle.txt";
DParamChaff=inputdir+"parameters/ParamChaff.txt";
DParamRoot=inputdir+"parameters/ParamRoot.txt";

#path of DrivingVariables 

# H0, H1, H2 refers to treatments H0, H3 and H15 (paper definition) 
DNsoilH0=inputdir+"DrivingVariables/NsoilH2.txt"; 
DTimeSoil=inputdir+"DrivingVariables/TimeSoil.txt";
Dmeteo=inputdir+"DrivingVariables/meteo(old2).txt";

#path of State0dd
DState0Grain=inputdir+"State0dd/State0Grain.txt";
DState0Lamina=inputdir+"State0dd/State0Lamina.txt";
DState0Sheath=inputdir+"State0dd/State0Sheath.txt";
DState0Internode=inputdir+"State0dd/State0Internode.txt";
DState0Peduncle=inputdir+"State0dd/State0Peduncle.txt";
DState0Chaff=inputdir+"State0dd/State0Chaff.txt";

#DState0Root=inputdir+"State0dd/State0Root.txt"; #ELMER

grain = Grain(DParamGrain, DState0Grain, outputdir)

root = Root(DParamRoot, outputdir)
#root = Root(DParamRoot, DStateRoot, outputdir)

environment = Environment(DParamEnv,DNsoilH0,DTimeSoil,Dmeteo)
lamina = Entity(DParamLamina,DState0Lamina,4, outputdir)
sheath = Entity(DParamSheath,DState0Sheath,4, outputdir)
internode = Entity(DParamInternode,DState0Internode,4, outputdir)
peduncle = Entity(DParamPeduncle,DState0Peduncle,1, outputdir)
chaff = Entity(DParamChaff,DState0Chaff,1, outputdir)

plant = Plant(grain, root, environment, lamina, sheath, internode, peduncle, chaff)

nSTEPS = 49
step = 1
nbsubsteps = 4

plant.InitialPlant(nbsubsteps)

for i in range(nSTEPS):
    plant.ActuPlant(step, nbsubsteps)
    print 'Simulation %d'%(i+1)

plant.Display(outputdir, nbsubsteps)
grain = plant.GetGrain()
root = plant.GetRoot()
lamina = plant.GetLamina()

print 'Simulation succeed'


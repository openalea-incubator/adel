from alinea.nema import *
from numpy import *
from pylab import *

"""TODO:
Readers for the files to extract parameters.
Readers for the initial values.
Build entities with these initial values (using addEntity).
Rewrote the methode Affich in Python  (named Display).
"""
inputdir = ''
outputdir = ''

default_par= { 

    'Grain' : {'TTexp': 700.0,
               'TTinit': 0.0,
               'alpha': 2.0,
               'beta': 2.0,
               'ggrain': 0.0055,
               'p': 6.0,
               'sink': 7.0,
               'tgrain': 250.0},
    
    'Chaff' : {'DegDMRate': 0.008,
               'DegRate': 0.008,
               'InsertionAngle': 1.57,#adel 
               'SynthRate': 3.0e-05,
               'TTexp': 1100.0,
               'alpha': 2.0,
               'beta': 2.0,
               'k1': 0.0018,
               'k2': 864000,
               'p': 0.1,
               'pdeath': 0.4,
               'peff': 5.0e-06,
               'pm1': 0.8,
               'pm2': 6.0,
               'sink': 2.0},
    
    'Internode' : {'DegDMRate': 0.008,
                   'DegRate': 0.008,
                   'InsertionAngle': 1.57,# adel : average for all internode
                   'SynthRate': 3.0e-05,
                   'TTexp': 920.0,
                   'alpha': 2.0,
                   'beta': 2.0,
                   'k1': 0.0018,
                   'k2': 864000.0,
                   'p': 0.1,
                   'pdeath': 0.4,
                   'peff': 5.0e-06,
                   'pm1': 0.8,
                   'pm2': 6.0,
                   'sink': 1.0},
 
    'Lamina' : {'DegDMRate': 0.008,
                'DegRate': 0.008,
                'InsertionAngle': 0.52,#adel average
                'SynthRate': 0.00015,
                'TTexp': 1090.0,
                'alpha': 2.0,
                'beta': 2.0,
                'k1': 0.0018,
                'k2': 864000.0,
                'p': 0.1,
                'pdeath': 0.4,
                'peff': 5.0e-06,
                'pm1': 0.8,
                'pm2': 6.0,
                'sink': 1.0},
    
    'Peduncle' : {'DegDMRate': 0.008,
                  'DegRate': 0.008,
                  'InsertionAngle': 1.57,#adel
                  'SynthRate': 3.0e-05,
                  'TTexp': 830.0,
                  'alpha': 2.0,
                  'beta': 2.0,
                  'k1': 0.0018,
                  'k2': 864000.0,
                  'p': 0.1,
                  'pdeath': 0.4,
                  'peff': 5.0e-06,
                  'pm1': 0.8,
                  'pm2': 6.0,
                  'sink': 2.0},
    
    'Root' : {'degDMcoef': 0.008,
              'degcoef': 0.008,
              'NsoilMin': 0.0,
              'QonDmin': 0.0013,
              'TTexp': 1700.0,
              'TTinit': -1000.0,
              'Umax': 0.0001,
              'alpha': 2.0,
              'beta': 2.0,
              'beta1': 130.0,
              'beta2': 2300.0,
              'kroot1': 2.5,
              'kroot2': 5.0e-06,
              'p': 0.1,
              'sink': 1.0},
    
    'Sheath' : {'DegDMRate': 0.008,
                'DegRate': 0.008,
                'InsertionAngle': 1.57,#adel
                'SynthRate': 0.00015,
                'TTexp': 940.0,
                'alpha': 2.0,
                'beta': 2.0,
                'k1': 0.0018,
                'k2': 864000.0,
                'p': 0.1,
                'pdeath': 0.4,
                'peff': 5.0e-06,
                'pm1': 0.8,
                'pm2': 6.0,
                'sink': 1.0}
    
    }
#reference is relative to dimensions of the mean plant measured by Jessica
default_init= {
    
    'Grain' : {'Ngrain': 0.0024, #adel if we choose a model Ngrain = Ngrainref/reference_volume_of_ear * adel_volume_of_ear
               'DMgrain': 0.12},#adel if we choose a model Ngrain = Ngrainref/reference_volume_of_ear * adel_volume_of_ear

    
    'Root' : {'DMrem': 0.056,# scale with whole plante volume ????????
              'DMstruct': 0.504,
              'Nrem': 0.0011,
              'Nstruct': 0.01},

    'Chaff' : {'Area': [0.00075], #adel
           'DMrem': [0.05], # reference value / reference_volume_of_ear * adel_volume_of_ear = 0.05/(0.00075*0.09) * adel_volume_of_ear
           'DMstruct': [0.21],# reference value / reference_volume_of_ear * adel_volume_of_ear
           'Length': [0.09],#adel
           'Nphs': [0.00368],# reference value / reference_area_of_ear * adel_volume_of_ear
           'Nstruct': [0.00107],# reference value / reference_area_of_ear * adel_volume_of_ear
           'TTinit': [-400]},#adel

    'Internode' : {'Area' : [0.0001,0.00025,0.0004,0.0015],#idem chaff
               'DMrem': [0.0093,0.033,0.039,0.041],
               'DMstruct': [0.0428,0.1541,0.1797,0.188],
               'Length': [0.05,0.086,0.128,0.186],
               'Nphs': [0.0,0.0,0.00028,0.00123],
               'Nstruct': [0.00005,0.00014,0.00033,0.00085],
               'TTinit': [-490.0,-400.0,-310.0,-220.0]},

    'Lamina' : {'Area' : [0.0016,0.00228,0.0034,0.00346],#idem chaff except that all is scaled with area
            'DMrem': [0.0,0.02,0.04,0.04],
            'DMstruct': [0.04,0.05,0.09,0.14],
            'Length': [0.182,0.211,0.2265,0.1735],
            'Nphs': [0.00009,0.0012,0.0029,0.0053],
            'Nstruct': [0.00038,0.00053,0.00083,0.00102],
            'TTinit': [-660.0,-570.0,-480.0,-390.0]},

    'Peduncle' : {'Area' : [0.0024],#idem chaff
              'DMrem': [0.0555],
              'DMstruct': [0.2568],
              'Length': [0.219],
              'Nphs': [0.00253],
              'Nstruct': [0.0013],
              'TTinit': [-130.0]},

    'Sheath' : {'Area' : [0.0002,0.0004,0.0005,0.0006],#idem chaff except that all is scaled with area
            'DMrem': [0.0019,0.0093,0.0148,0.0222],
            'DMstruct': [0.0086,0.0428,0.0685,0.1027],
            'Length': [0.11,0.125,0.14,0.145],
            'Nphs': [0.000026,0.00018,0.00066,0.0018],
            'Nstruct': [0.00012,0.00011,0.00021,0.00068],
            'TTinit': [-510.0,-420.0,-330.0,-240.0]}

}

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

def read_state0(filename):
    f = open(filename)
    def first_float(l):
        l = l.strip()
        return float(l.split()[0]) if l else None
    states = [first_float(l) for l in f if first_float(l) is not None]
    
    f.close()
    return states

def update_entity(entity, state0):
    states = read_state0(state0)
    n = len(states)/7
    for i in range(n):
        Nstruct,Nph,Length,Area,TTinit,DMrem,DMstruct = states[7*i:7*(i+1)]
        entity.addEntity(Nstruct,Nph,Length,Area,TTinit,DMrem,DMstruct)

def setInitGrain(grain, state0dict):
	grain.N0 = state0dict['Ngrain']
	grain.DM0 = state0dict['DMgrain']
	
def setInitRoot(root, state0dict):
	root.DMrem0 = state0dict['DMrem']
	root.DMstruct0 = state0dict['DMstruct']
	root.Nrem0 = state0dict['Nrem']
	root.Nstruct0 = state0dict['Nstruct']
	
	
def setInitEntity(entity, state0dict):
    n = len(state0dict['Area'])
    for v in state0dict.itervalues():
        assert len(v) == n
    entity.Areas0 = state0dict['Area']
    entity.DMrems0 = state0dict['DMrem']
    entity.DMstruct0 = state0dict['DMstruct']
    entity.Lengths0 = state0dict['Length']
    entity.Nphs0 = state0dict['Nphs']
    entity.Nstructs0 = state0dict['Nstruct']
    entity.TTinits0 = state0dict['TTinit']

#grain = Grain(DParamGrain, DState0Grain, outputdir)
grain = Grain()
grain.params = default_par['Grain']
setInitGrain(grain,default_init['Grain'])

#root = Root(DParamRoot, outputdir)
root = Root()
root.params = default_par['Root']
setInitRoot(root,default_init['Root'])


#chaff = Chaff (DParamChaff, outputdir)
chaff = Entity()
chaff.params = default_par['Chaff']
setInitEntity(chaff,default_init['Chaff'])

#internode = Internode (DParamInternode, outputdir)
internode = Entity()
internode.params = default_par['Internode']
setInitEntity(internode,default_init['Internode'])

#lamina = Lamina (DParamLamina, outputdir)
lamina = Entity()
lamina.params = default_par['Lamina']
setInitEntity(lamina,default_init['Lamina'])

#peduncle = Peduncle (DParamPeduncle, outputdir)
peduncle = Entity()
peduncle.params = default_par['Peduncle']
setInitEntity(peduncle,default_init['Peduncle'])

#sheath = sheath (DParamsheath, outputdir)
sheath = Entity()
sheath.params = default_par['Sheath']
setInitEntity(sheath,default_init['Sheath'])


 
#environment = Environment(DParamEnv,DNsoilH0,DTimeSoil,Dmeteo)
environment = Environment()

#lamina = Entity(DParamLamina,DState0Lamina,0, outputdir)
#update_entity(lamina, DState0Lamina)

#sheath = Entity(DParamSheath,DState0Sheath,0, outputdir)
#update_entity(sheath, DState0Sheath)

#internode = Entity(DParamInternode,DState0Internode,0, outputdir)
#update_entity(internode, DState0Internode)

#peduncle = Entity(DParamPeduncle,DState0Peduncle,0, outputdir)
#update_entity(peduncle, DState0Peduncle)

#chaff = Entity(DParamChaff,DState0Chaff,0, outputdir)
#update_entity(chaff, DState0Chaff)


plant = Plant(grain, root, environment, lamina, sheath, internode, peduncle, chaff)

nSTEPS = 49
step = 1
dTT = 10
nbsubsteps = 4
Nmob = 0.0024
Nsoil = 1
PARs = [1] * 4
PAR = [0]
PARlamina = PARs
PARsheath = PARs
PARinternode = PARs
PARpeduncle = PAR
PARchaff = PAR
PARss = [PARlamina,PARsheath,PARinternode,PARpeduncle,PARchaff]


#plant.InitialPlantExt(nbsubsteps, Nmob,PARss, Nsoil)
#for i in range(nSTEPS):
#    plant.ActuPlantExt(step, nbsubsteps, dTT, PARss)

#print 'Simulation succeed'

#print 'Nb lamina ', len(lamina)

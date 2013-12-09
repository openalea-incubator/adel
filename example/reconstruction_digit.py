""" Tutorial Reconstructing canopy from digitised data """

import pandas
import numpy

#input
L_ear = 7.2
L_ped = 14.8
W_ped = 0.2
A_ear = 10.28
L_spike = 8.8
df_dim = pandas.DataFrame({ 'id_plt' : [1] * 9,
                        'ntop' : [6,5,4,3,2,1, 0, -1, -2],
                        'L_blade' : [0] * 3 + [17.07,12.1,8.87] + 3 * [0],
                        'W_blade' : [0] * 3 + [1,1.23,1.01] + 3 * [0],
                        'L_sheath' : [0] * 3 + [11.21,12.55,13.59] + 3 * [0],
                        'W_sheath' : [0] * 3 + [0.39,.49,.36] + [W_ped, A_ear / L_ear, A_ear / L_ear],
                        'L_internode' : [0.55,4.97,8.28,13.06,13.36,18.32] + [L_ped, L_ear, L_spike - L_ear], 
                        'W_internode' : [0] * 3 + [0.34, 0.44, 0.31] + [W_ped, A_ear / L_ear, A_ear / L_ear] # 
                    })
        
#geometry (to be replaced by calling fit_leaves on your data)
#from alinea.adel.data_samples import leaves_db 
#db = leaves_db()
from alinea.adel.wheat import extract_wheat
db = extract_wheat.extract_leaf_info('laminaCurv.RData', 'lamina2D.RData')
from alinea.adel.fit import fit
db, discard = fit.fit_leaves(db, 7)
       
# compute insertion height of leaves
df_dim['h_insertion'] = df_dim['L_internode'].cumsum() + df_dim['L_sheath']
#keep only uper metamer
df_dim = df_dim[df_dim.ntop <= 3]
# create pseudo-stem element (spacer between leaves)
df_dim['pstem'] = numpy.diff([0] + df_dim['h_insertion'].tolist())

# construction of canopy table
df = df_dim[['id_plt', 'ntop', 'L_blade', 'W_blade', 'pstem', 'W_sheath']]
df.rename(columns={'id_plt':'plant', 'L_blade':'Ll', 'W_blade':'Lw_shape', 'pstem': 'El', 'W_sheath': 'Ed'}, inplace=True)
# add mandatory topological info
df['axe_id'] = 'MS'
df['ms_insertion'] = 0
df['numphy'] = df_dim.ntop.max() + 1 - df_dim.ntop # odre from base to top
# add missing mandatory data  (does like adel)                  
df['Laz'] = [180 + (numpy.random.random() - 0.5) * 30 for i in range(len(df))] #leaf azimuth
df['LcType'] = df['ntop'] # selector for first level in leaf db
df['Lc_index'] = 1 # selector for second level in leaf_db (ranging 1:max_nb_leaf_per_level)
# fill other columns
df['Lv'] = df['Ll']
df['Lr'] = 0
df['Lsen'] = 0
df['L_shape'] = df['Ll']
df['Linc'] = 1
df['Gv'] = 0
df['Gsen'] = 0
df['Ginc'] = 0
df['Ev'] = df['El']
df['Esen'] = 0
df['Einc'] = 0

#create mtg
from alinea.adel.stand.stand import agronomicplot
from alinea.adel.newmtg import mtg_factory, adel_metamer, adel_label
from alinea.adel.mtg_interpreter import mtg_interpreter, plot3d
from openalea.plantgl.all import Viewer

# create a canopy with the same number of plants as in the canopy table (play with width, length)
nplants, positions, domain, domain_area, convUnit = agronomicplot(length=0.1, 
                                                            width=0.2, 
                                                            sowing_density=150, 
                                                            plant_density=150,
                                                            inter_row=0.12)
assert(nplants == df.plant.max())
stand = [(pos,0) for pos in positions]
g = mtg_factory(df.to_dict('list'), adel_metamer, leaf_db = db, stand = stand )
#add geometry
g = mtg_interpreter(g)
# plot
scene = plot3d(g)
Viewer.display(scene)

#call caribu
from alinea.caribu.caribu_star import caribu_star
geom = g.property('geometry')
star, exposed_area = caribu_star(geom, directions = 16, domain = domain, convUnit = convUnit)#cf caribu_star doc for output interpretation
res = pandas.DataFrame([(adel_label(g,vid), star[vid], exposed_area[vid]) for vid in star])


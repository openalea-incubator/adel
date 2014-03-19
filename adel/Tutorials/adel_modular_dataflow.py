geometric_data_directory_path = '../data'
botanic_data_directory_path = '../data'
pov_dirpath = '.'
import tempfile
count_pixels_result_filepath = tempfile.mktemp(prefix='maxwell_ref_outputpovray_', suffix='.csv')
caribu_postprocessing_directory_path = tempfile.mkdtemp(prefix='_maxwell_ref_outputcaribu_postprocessing')
import os
LAI_filepath = os.path.join(caribu_postprocessing_directory_path, 'maxwell_ref_caribuoutput_.csv')
global_postprocessing_filepath = os.path.join(caribu_postprocessing_directory_path, 'maxwell_ref_postprocessing_global__RUN2_.csv')
peraxis_postprocessing_filepath = os.path.join(caribu_postprocessing_directory_path, 'maxwell_ref_postprocessing_peraxis_RUN2_.csv')
measurement_filepath = '../data/Maxwell_2011_d220_N1_LAI-PAI_mesure.csv'

# functions implemented in __my package__. To move to alinea.adel ??? 
def AND(xin, bool):
        '''    This function comines two vectors (given as list) one as float and boolean
        '''
        # write the node code here.
        bool=[x-2 for x in bool];
        bool=[abs(x) for x in bool];
        xout=[x*y for x,y in zip(bool,xin)]
        # return outputs
        return xout

def Write_LAI(Tsum, Gc, CAI, CAIg, LAI, LAIg, LAI_filepath):
        '''
        returns a string representation of the content of the table
        '''
        fout = open(LAI_filepath, "a")
        fout.write(str(Tsum))
        fout.write(' , ')
        fout.write(str(Gc))
        fout.write(' , ')
        fout.write(str(CAI))
        fout.write(' , ')
        fout.write(str(CAIg))
        fout.write(' , ')
        fout.write(str(LAI))
        fout.write(' , ')
        fout.write(str(LAIg))
        fout.write('\n')
        fout.close()
        return LAI_filepath

def BOUCLE(in1, nbr_simulation):
    b=[]
    j=0
    while j < len(in1):
        i=0
        while i < nbr_simulation:
            b.append (in1[j])
            i = i + 1
        j = j + 1
    out1 = b; 
    return out1

# select_adel_geometric_data
laminaCur_filename = 'Maxwell_2011_d220_N1_laminaCurv_global.RData'
lamina2D_filename = 'Maxwell_2011_d220_N1_lamina2D_global.RData'
laminaCur_filepath = os.path.join(geometric_data_directory_path, laminaCur_filename)
lamina2D_filepath = os.path.join(geometric_data_directory_path, lamina2D_filename)

# extract_leaves
import openalea
from alinea.adel.wheat import extract_wheat
leaves = extract_wheat.extract_leaf_info(laminaCur_filepath, lamina2D_filepath)

# fit leaves
from alinea.adel.fit import fit
leaves, discard = fit.fit_leaves(leaves, 7)

# getBotanicData, composed of select_adel_geometric_data and devCsv
# select_adel_geometric_data
axeT_filename = 'Maxwell_2011_d220_N1_axeT.csv'
dimT_filename = 'Maxwell_2011_d220_N1_dimT.csv'
phenT_filename = 'Maxwell_2011_d220_N1_phenT.csv'
earT_filename = 'Maxwell_2011_d220_N1_earT.csv'
ssi2sen_filename = 'ssi2sen.csv'

# devCsv
from alinea.adel import AdelR
development_parameters = AdelR.devCsv(os.path.join(botanic_data_directory_path, axeT_filename), 
                                      os.path.join(botanic_data_directory_path, dimT_filename), 
                                      os.path.join(botanic_data_directory_path, phenT_filename), 
                                      os.path.join(botanic_data_directory_path, earT_filename), 
                                      os.path.join(botanic_data_directory_path, ssi2sen_filename))

# geoAxe
axes_geometry = AdelR.genGeoAxe(75.00, 0.0, 0.0, 0.0, 60.0, 5.0, 7.0)

# geoLeaf
leaf_geometry = AdelR.genGeoLeaf(4, 60.0, 10.0)

# agronomic plot
inter_row = 17.5 / 100.0
plot_width = 4.0 * inter_row
from alinea.adel.stand import stand
number_of_plants, positions, domain, domain_area, conv_coeff = stand.agronomicplot(0.5, plot_width, 220.0, 200.0, inter_row, 0.0, 100, True)

# genMTG_3out, composed of symbols, setAdel, RunAdel, Rdataframe, genString, 
# CanMTGInterpreter, CanMTGPlanter
# return stand_MTG, rdataframe and Lstrings
# symbols
from alinea.adel.geometry import geometry
symbols_dict = geometry.symbols(leaves, 0)
# setAdel
simulation_parameters = AdelR.setAdel(development_parameters, leaf_geometry, axes_geometry, number_of_plants, 0, laminaCur_filepath, lamina2D_filepath, 'sequence')

# ThermalTimeInterval
#thermal_times = range(0, 2000, 10)
thermal_times = range(0, 10, 10)
import numpy as np
repeated_thermal_times = np.repeat(thermal_times, 1)

# loop on repeated thermal times
for thermal_time in repeated_thermal_times:
    # RunAdel
    global_parameters = {'endLeaf': 1.6, 'endLeaf1': 1.6, 'epsillon': 1e-06, 'startLeaf': -0.4, 'senescence_leaf_shrink': 0.5, 'stemLeaf': 1.2}
    Lstrings = AdelR.RunAdel(thermal_time, simulation_parameters, global_parameters)
    # Rdataframe
    from alinea.adel.io import io
    rdataframe = io.dataframe(Lstrings)
    # genString
    lsystem_string = AdelR.genString(rdataframe)
    # CanMTGInterpreter
    from alinea.adel.wheat import CanMTGInterpreter
    can_MTG = CanMTGInterpreter.CanMTGInterpreter(symbols_dict, lsystem_string)
    # CanMTGPlanter
    from alinea.adel.stand import CanMTGPlanter
    stand_MTG = CanMTGPlanter.CanMTGPlanter(can_MTG, positions)
    
    # colors
    colors = [[0, 255, 0], [255, 0, 0], [0, 255, 0], [255, 0, 0], [0, 255, 0], [255, 0, 0], [0, 255, 0], [255, 0, 0], [0, 255, 0], [255, 0, 0], [0, 255, 0], [255, 0, 0]]
    # col_item
    from alinea.adel.povray import povray
    lambda_function = povray.col_item(None, colors)
    # apply_property
    from alinea.adel import mtg
    new_values = mtg.apply_property(stand_MTG, 'tissue_type', lambda_function)
    # to_plantgl
    plantgl_scene = mtg.to_plantgl(stand_MTG, (0,180,0), (0, 130, 0), (170,85,0), new_values, True)
    # plot3D
    from openalea.plantgl.wralea.visualization import viewernode
    viewernode.Plot3D(plantgl_scene)
    
    # povray
    pov_filepath = os.path.join(pov_dirpath, '%s%s%s' % ('scene_', thermal_time, '.pov'))
    povray_image_filepath, stand_box_image_filepath = povray.povray(plantgl_scene, pov_filename, 200.0, 50, 4288, 2848, domain, 0, 0, 'perspective', False, 'povray')
    
    # genFilepath
    from alinea.adel.povray import post_processing
    post_processing.count_pixels(povray_image_filepath, stand_box_image_filepath, colors, count_pixels_result_filepath)
    
    # light string
    from alinea.adel.caribu import lightString
    directional_vector = lightString.lightString(1.0, 0.0, 46.0)
    
    # par4.opt
    from openalea.core import pkgmanager
    pm = pkgmanager.PackageManager()
    pm.init(verbose=False)
    from openalea.core.system import systemnodes
    opt_filepath = systemnodes.get_data('par4.opt', 'alinea.caribu.data').values().pop()
    
    # read
    from openalea.file import files
    opt_str = files.FileRead()(opt_filepath)
    
    # CaribuScene
    from alinea.caribu import CaribuScene_nodes
    caribu_scene = CaribuScene_nodes.newCaribuScene(stand_MTG, directional_vector, domain, opt_str)
    
    # addsoil_2
    caribu_scene = CaribuScene_nodes.addSoil(caribu_scene, 2.0)
    
    # Caribu
    caribu_scene, energy = CaribuScene_nodes.runCaribu(caribu_scene, True, {'SphereDiameter': 0.5, 'Nz': 5, 'Zmax': 2, 'keepFF': False}, True)
    
    # LIE
    from alinea.caribu import selectOutput
    canestra_output_1 = selectOutput.selectOutput(energy, 'Opt')[0]
    canestra_output_2 = selectOutput.selectOutput(energy, 'Einc')[0]
    from alinea.caribu import filterby
    filtered = filterby.filterby(canestra_output_1, canestra_output_2, lambda x: x == 0.0)
    total_soil = sum(filtered)
    total_incident = CaribuScene_nodes.getIncidentEnergy(caribu_scene)[2]
    efficience = (total_incident - total_soil) / float(total_incident)
    
    # LAIg
    canestra_output_1 = selectOutput.selectOutput(energy, 'Opt')[0]
    canestra_output_2 = selectOutput.selectOutput(energy, 'Area')[0]
    soil_filtered = filterby.filterby(canestra_output_1, canestra_output_2, lambda x: x == 0.0)
    soil_area = sum(soil_filtered)
    
    canestra_output_1 = selectOutput.selectOutput(energy, 'Opt')[0]
    canestra_output_2 = selectOutput.selectOutput(energy, 'Area')[0]
    canopy_triangles_filtered = filterby.filterby(canestra_output_1, canestra_output_2, lambda x: x > 0.0)
    canopy_triangles_area = sum(canopy_triangles_filtered)
    
    canestra_output_1 = selectOutput.selectOutput(energy, 'Opt')[0]
    canestra_output_2 = selectOutput.selectOutput(energy, 'Area')[0]
    green_triangles_filtered = filterby.filterby(canestra_output_1, canestra_output_2, lambda x: x == 1.0)
    green_triangles_area = sum(green_triangles_filtered)
    
    canestra_output_1 = selectOutput.selectOutput(energy, 'Opak')[0]
    canestra_output_2 = selectOutput.selectOutput(energy, 'Area')[0]
    opak_leaves_filtered = filterby.filterby(canestra_output_1, canestra_output_2, lambda x: x > 0.0)
    opak_leaves_area = sum(opak_leaves_filtered)
    
    canestra_output_1 = selectOutput.selectOutput(energy, 'Opak')[0]
    canestra_output_2 = selectOutput.selectOutput(energy, 'Opt')[0]
    opak_opt_filtered = filterby.filterby(canestra_output_1, canestra_output_2, lambda x: x > 0.0)
    
    LAI = opak_leaves_area / float(soil_area)
    
    LAIg = sum(AND(opak_leaves_filtered, opak_opt_filtered)) / float(soil_area)
    
    CAI = (canopy_triangles_area / float(soil_area) - LAI) / 2 + LAI
    
    CAIg = (green_triangles_area / float(soil_area) - LAIg) / 2 + LAIg
    
    # Write_LAI
    LAI_filepath = Write_LAI(total_incident, efficience, CAI, CAIg, LAI, LAIg, LAI_filepath)
    
    # canL2canS
    rdataframe = AdelR.canL2canS(rdataframe, lamina2D_filepath, 0.5)
    
    # update
    Lstrings.update(rdataframe)
    
    # Write_Table
    from alinea.adel.wheat import Write_Table
    canl2canS_filepath = os.path.join(caribu_postprocessing_directory_path, '%s%s%s' % ('maxwell_ref_canl2canS_RUN2_', thermal_time, '.csv'))
    adel_output_filepath = Write_Table.Write_Table(Lstrings, canl2canS_filepath, sep=',')
    
    # post_processing
    (global_postprocessing_filepath, 
     peraxis_postprocessing_filepath, 
     intermediate_postprocessing_filepath) = stand.post_processing(adel_output_filepath, 
                                                                   number_of_plants, 
                                                                   domain_area, 
                                                                   global_postprocessing_filepath, 
                                                                   peraxis_postprocessing_filepath)
     
    # plot_graph, composed of select_variable, PyLabLine2D, PyLabPlot and PyLabLegend.
    # return axes
    # select_variable, composed of CsvAsDict
    global_postprocessing_dict = csvAsDict(global_postprocessing_filepath)
    global_postprocessing_ThermalTime = global_postprocessing_dict['ThermalTime']
    global_postprocessing_LAI_tot = global_postprocessing_dict['LAI_tot']
    measurement_dict = csvAsDict(measurement_filepath)
    measurement_ThermalTime = measurement_dict['ThermalTime']
    measurement_LAI_tot = measurement_dict['LAI_tot']
    # PyLabLine2D
    from openalea.pylab_plotting_wralea import py_pylab as plotting_py_pylab
    global_postprocessing_pll2D = plotting_py_pylab.PyLabLine2D()
    global_postprocessing_pll2D.set_input('xdata', global_postprocessing_ThermalTime)
    global_postprocessing_pll2D.set_input('ydata', global_postprocessing_LAI_tot)
    global_postprocessing_pll2D.set_input('linestyle', 'dashed')
    global_postprocessing_pll2D.set_input('color', 'black')
    global_postprocessing_pll2D.set_input('marker', 'circle')
    global_postprocessing_pll2D.set_input('markersize', 8)
    global_postprocessing_pll2D.set_input('markeredgewidth', 0.0)
    global_postprocessing_pll2D.set_input('markeredgecolor', 'None')
    global_postprocessing_pll2D.set_input('linewidth', 1.0)
    global_postprocessing_pll2D.set_input('fillstyle', 'full')
    global_postprocessing_pll2D.set_input('label', 'simulation')
    global_postprocessing_pll2D.set_input('alpha', 1.0)
    global_postprocessing_line2d = global_postprocessing_pll2D(inputs=None)
    
    measurement_pll2D = plotting_py_pylab.PyLabLine2D()
    measurement_pll2D.set_input('xdata', measurement_ThermalTime)
    measurement_pll2D.set_input('ydata', measurement_LAI_tot)
    measurement_pll2D.set_input('linestyle', 'dashed')
    measurement_pll2D.set_input('color', 'white')
    measurement_pll2D.set_input('marker', 'circle')
    measurement_pll2D.set_input('markersize', 8)
    measurement_pll2D.set_input('markeredgewidth', 0.0)
    measurement_pll2D.set_input('markeredgecolor', 'None')
    measurement_pll2D.set_input('linewidth', 1.0)
    measurement_pll2D.set_input('fillstyle', 'full')
    measurement_pll2D.set_input('label', 'mesure')
    measurement_pll2D.set_input('alpha', 1.0)
    measurement_line2d = measurement_pll2D(inputs=None)
    # PyLabPlot
    plp = plotting_py_pylab.PyLabPlot()
    plp.set_input('x', [global_postprocessing_line2d, measurement_line2d])
    plp.set_input('marker', 'circle')
    plp.set_input('markersize', 10.0)
    plp.set_input('linestyle', 'solid')
    plp.set_input('color', 'blue')
    plp.set_input('scalex', True)
    plp.set_input('scaley', True)
    plp.set_input('kwargs', {'label': None})
    plp.set_input('figure', 1)
    axes = plp(inputs=None)
    # PyLabLegend
    from openalea.pylab_decorators_wralea import py_pylab as decorators_py_pylab
    pllegend = decorators_py_pylab.PyLabLegend()
    pllegend.set_input('axes', axes)
    pllegend.set_input('shadow', True)
    pllegend.set_input('location', 'upper right')
    pllegend.set_input('numpoints', 1)
    pllegend.set_input('markerscale', 0.8)
    pllegend.set_input('fancybox', True)
    pllegend.set_input('ncol', 1)
    pllegend.set_input('mode', 'None')
    pllegend.set_input('title', 'LAI_tot')
    axes = pllegend(inputs=None)







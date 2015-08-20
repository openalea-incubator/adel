# -*- python -*-
#
#       Adel.PlantGen
#
#       Copyright 2012-2014 INRIA - CIRAD - INRA
#
#       File author(s): Camille Chambon <camille.chambon@grignon.inra.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################
'''
Routines to plot the outputs of :mod:`alinea.adel.plantgen`. 

Authors: M. Abichou, B. Andrieu, C. Chambon
'''

import os

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import numpy as np
import pandas as pd

def plot_HS_GL_SSI_T(HS_GL_SSI_T, id_phen_to_plot=None, plots_dirpath=None):
    '''
    Plot the HS, GL and SSI of `id_phen_to_plot`.
    
    :Parameters:
    
        - `HS_GL_SSI_T` (:class:`pd.DataFrame`) - the table HS_GL_SSI_T.
        - `id_phen_to_plot` (:class:`list`) - the list of id_phen to plot. If None (the default), then plot all the id_phen. 
        - `plots_dirpath` (:class:`str`) - the path of the directory to save the graphs in.  
          If `None` (the default), do not save the graphs but display them.
          
    :Examples:
        
        # plot HS, GL and SSI of id_phen 1101, 1111 and 4081
        import pandas as pd
        HS_GL_SSI_T = pd.read_csv('HS_GL_SSI_T.csv')
        from alinea.adel.plantgen import graphs
        graphs.plot_HS_GL_SSI_T(HS_GL_SSI_T, id_phen_to_plot=[4081, 1111, 1101])
        
    '''
    if id_phen_to_plot is None:
        HS_GL_SSI_T_to_plot = HS_GL_SSI_T
    else:
        HS_GL_SSI_T_to_plot = HS_GL_SSI_T[HS_GL_SSI_T['id_phen'].isin(id_phen_to_plot)]
    
    for id_phen, group in HS_GL_SSI_T_to_plot.groupby('id_phen'):
        group.index = group.TT
        group = group.drop(['TT', 'id_phen'], axis=1)
        current_axis = group.plot()
        current_axis.set_xlabel('Thermal time')
        current_axis.set_ylabel('Decimal leaf number')
        current_axis.legend(prop={'size':10}, framealpha=0.5)
        current_axis.set_title('id_phen = {}'.format(id_phen))
        xmin, xmax = current_axis.get_xlim()
        x_margin = (xmax - xmin) / 100.0
        current_axis.set_xlim(xmin - x_margin, xmax + x_margin)
        ymin, ymax = current_axis.get_ylim()
        y_margin = (ymax - ymin) / 100.0
        current_axis.set_ylim(ymin - y_margin, ymax + y_margin)
        if plots_dirpath is None:
            plt.show()
        else:
            plt.savefig(os.path.join(plots_dirpath, 'HS_GL_SSI_T_{}.png'.format(id_phen)), dpi=200, format='PNG')
            plt.close()


def plot_dimT(dimT, MS_id_dim=None, relative_index_phytomer=False, dimensions_to_plot=None, id_dim_to_plot=None, id_cohort_to_plot=None, plot_non_regressive_tillers=True, plot_regressive_tillers=True, plots_dirpath=None):
    '''
    Plot the dimensions in `dimT` according to filters `dimensions_to_plot`, `id_dim_to_plot`, `id_cohort_to_plot`, `plot_non_regressive_tillers` and `plot_regressive_tillers`.
    
    :Parameters:
    
        - `dimT` (:class:`pd.DataFrame`) - the table dimT.
        - `MS_id_dim` (:class:`int`) - the id_dim of the main stem. The line corresponding to this id_dim is thicker.
        - `relative_index_phytomer` (:class:`bool`) - if True: display the relative index of the phytomers in abscissa. 
          If False (the default), display the absolute index of the phytomers in abscissa.
        - `dimensions_to_plot` (:class:`list`) - the list of dimensions to plot. If None (the default), then plot all the dimensions.
        - `id_dim_to_plot` (:class:`list`) - the list of id_dim to plot. If None (the default), then plot all the id_dim.
        - `id_cohort_to_plot` (:class:`list`) - the list of id_cohort to plot. If None (the default), then plot all the id_cohort.
        - `plot_non_regressive_tillers` (:class:`bool`) - whether to plot the non regressive tillers or not. Non regressive tillers have id_dim ending by '1'. Default is to plot the non regressive tillers.
        - `plot_regressive_tillers` (:class:`bool`) - whether to plot the regressive tillers or not. Regressive tillers have id_dim ending by '0'. Default is to plot the regressive tillers. 
        - `plots_dirpath` (:class:`str`) - the path of the directory to save the graphs in.  
          If `None`, do not save the graphs but display them.
          
    :Examples:
        
        # plot L_blade
        import pandas as pd
        dimT = pd.read_csv('dimT.csv')
        from alinea.adel.plantgen import graphs
        graphs.plot_dimT(dimT, dimensions_to_plot=['L_blade'])
        
    '''
    
    DIM_T_KEY = ['id_dim', 'index_phytomer']
    
    MARKERS = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != ' ':
                MARKERS.append(m)
        except TypeError:
            pass
        
    MARKERS.extend([r'$\lambda$',
                    r'$\bowtie$',
                    r'$\circlearrowleft$',
                    r'$\clubsuit$',
                    r'$\checkmark$'])
    
    COLORS = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
    
    if dimensions_to_plot is None:
        dimT_to_plot = dimT
        dimensions_to_plot = dimT.columns.difference(DIM_T_KEY)
    else:
        dimT_to_plot = dimT[DIM_T_KEY + dimensions_to_plot]
    
    if id_dim_to_plot is None:
        dimT_to_plot = dimT_to_plot
    else:
        dimT_to_plot = dimT_to_plot[dimT_to_plot.id_dim.isin(id_dim_to_plot)]
        
    if id_cohort_to_plot is not None:
        dimT_to_plot = dimT_to_plot[dimT_to_plot.id_dim.astype(str).str[0].astype(int).isin(id_cohort_to_plot)]
    
    if not plot_non_regressive_tillers:
        dimT_to_plot = dimT_to_plot[dimT_to_plot.id_dim.astype(str).str[-1].astype(int) != 1]
        
    if not plot_regressive_tillers:
        dimT_to_plot = dimT_to_plot[dimT_to_plot.id_dim.astype(str).str[-1].astype(int) != 0]
        
    xlabel = 'index_phytomer'
    if relative_index_phytomer:
        xlabel = 'relative_' + xlabel
    
    for dimension in dimensions_to_plot:
        plt.figure()
        current_axis = plt.subplot(111)
        dimension_to_plot = dimT_to_plot[DIM_T_KEY + [dimension]]
        grouped = dimension_to_plot.groupby('id_dim')
        axis_num = 0
        for id_dim, group in grouped:
            axis_num += 1
            if relative_index_phytomer:
                index_phytomer_to_plot = group.index_phytomer / group.index_phytomer[group.last_valid_index()]
            else:
                index_phytomer_to_plot = group.index_phytomer
            if id_dim == MS_id_dim:
                linewidth = 4
                markersize = 10
            else:
                linewidth = 1
                markersize = 7
                
            color = COLORS[axis_num % len(COLORS)]
            marker = MARKERS[axis_num % len(MARKERS)]
            current_axis.plot(index_phytomer_to_plot, group[dimension], marker=marker, color=color, label=id_dim, linewidth=linewidth, markersize=markersize)
        
        current_axis.set_xlabel(xlabel)
        current_axis.set_ylabel(dimension)
        current_axis.legend(prop={'size':10}, framealpha=0.5)
        current_axis.set_title(dimension)
        xmin, xmax = current_axis.get_xlim()
        x_margin = (xmax - xmin) / 100.0
        current_axis.set_xlim(xmin - x_margin, xmax + x_margin)
        ymin, ymax = current_axis.get_ylim()
        y_margin = (ymax - ymin) / 100.0
        current_axis.set_ylim(ymin - y_margin, ymax + y_margin)
        
        if plots_dirpath is None:
            plt.show()
        else:
            plt.savefig(os.path.join(plots_dirpath, '{}.png'.format(dimension)), dpi=200, format='PNG')
            plt.close()


def plot_tillering_dynamic(axeT, plants_density, TT_step=10, plots_dirpath=None):
    '''
    Plot the dynamic of tillering, i.e. the evolution of the density of active axes when TT varies.
    
    A non regressive axis is active at TT if TT >= TT_em_phytomer1. 
    A regressive axis is active at TT if TT_em_phytomer1 <= TT < TT_stop_axis.  
    
    :Parameters:
    
        - `axeT` (:class:`pd.DataFrame`) - the table axeT.
        - `plants_density` (:class:`int`) - the number of plants per square meter.
        - `TT_step` (:class:`int`) - the thermal time step of the plot. Default is 10.
        - `plots_dirpath` (:class:`str`) - the path of the directory to save the graphs in.  
          If `None`, do not save the graphs but display them.
          
    :Examples:
        
        # plot tillering dynamic
        import pandas as pd
        axeT = pd.read_csv('axeT.csv')
        from alinea.adel.plantgen import graphs
        graphs.plot_tillering_dynamic(axeT, plants_density=250)
        
    '''
    
    number_of_plants = axeT.id_plt.nunique()
    domain_area = number_of_plants / float(plants_density)
    
    xmin = axeT.TT_em_phytomer1.min() - TT_step
    xmax = axeT.TT_stop_axis.max() + TT_step
    if TT_step >= 10:
        xmin = round(xmin, -1)
        xmax = round(xmax, -1)
    TT_grid = np.arange(xmin, xmax + TT_step, TT_step)
    
    axeT_non_regressive = axeT[axeT.id_ear.notnull()]
    null_id_ear_index = axeT.index.difference(axeT_non_regressive.index)
    axeT_regressive = axeT.loc[null_id_ear_index, :]
    
    densities_of_active_axes_per_square_meter = []
    densities_of_active_axes_per_plant = []
    
    for TT in TT_grid:
        number_of_non_regressive_active_axes = len(axeT_non_regressive[TT >= axeT_non_regressive.TT_em_phytomer1])
        number_of_regressive_active_axes = len(axeT_regressive[(TT >= axeT_regressive.TT_em_phytomer1) & (TT < axeT_regressive.TT_stop_axis)])
        total_number_of_active_axes = number_of_non_regressive_active_axes + number_of_regressive_active_axes
        densities_of_active_axes_per_square_meter.append(total_number_of_active_axes / domain_area)
        densities_of_active_axes_per_plant.append(total_number_of_active_axes / number_of_plants)
    
    for densities_of_active_axes, type_of_density, file_suffix in ((densities_of_active_axes_per_square_meter, 'per square meter', 'per_square_meter'),
                                                                   (densities_of_active_axes_per_plant, 'per plant', 'per_plant')):
        plt.figure()
        current_axis = plt.subplot(111)
        densities_of_active_axes_series = pd.Series(densities_of_active_axes, index=TT_grid)
        current_axis = densities_of_active_axes_series.plot()
        current_axis.set_xlabel('Thermal time')
        current_axis.set_ylabel('Density of active axes ' + type_of_density)
        current_axis.set_title('Tillering dynamic')
        xmin, xmax = current_axis.get_xlim()
        x_margin = (xmax - xmin) / 100.0
        current_axis.set_xlim(xmin - x_margin, xmax + x_margin)
        ymin, ymax = current_axis.get_ylim()
        y_margin = (ymax - ymin) / 100.0
        current_axis.set_ylim(ymin - y_margin, ymax + y_margin)
        if plots_dirpath is None:
            plt.show()
        else:
            plt.savefig(os.path.join(plots_dirpath, '{}_{}.png'.format('tillering_dynamic', file_suffix)), dpi=200, format='PNG')
            plt.close()
    

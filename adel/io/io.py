# Standard import

import cPickle as Pickle

import numpy as np
from rpy2 import robjects
from rpy2.robjects.numpy2ri import numpy2ri

from alinea.adel.AdelR import RlistAsDict,readRData,saveRData,csvAsDict,dataframeAsdict,canL2canS
from alinea.adel.mtg import (mtg_factory, duplicate, thermal_time, 
                             apply_property, to_plantgl, to_canestra)
from openalea.mtg.io import lpy2mtg, mtg2lpy

from openalea.core.path import path

def dataframe(d):
    """ convert a dict of numbers to an RDataframe  """
    df = {}
    if d is None:
        return robjects.r('as.null()')
    else:
        for k, v in d.iteritems():
            df[k] = numpy2ri(np.array(v))
    dataf = robjects.r['data.frame'](**df)
    return dataf


def load_leaf_data(fn):
    """  Load leaf data obtained by measurement. """ 
    leaves = {}
    try:
        f = open(fn)
        leaves = Pickle.load(f)
    finally:
        f.close()

    return leaves,


def csv2pandasDataframe(csv_filepath, index_col=None, na_values=None, parse_dates=False):
    '''Create a dataframe from a CSV (comma-separated) file.
    
    :Parameters:
        - `csv_filepath` : The filepath of the csv file to import.
        - `index_col` : Column to use as the row labels of the dataframe. If a sequence is given, a MultiIndex is used.
        - `na_values` : List of additional strings to recognize as NA/NaN.
        - `parse_dates` : Attempt to parse dates in the index column(s).
    
    :Types:
        - `csv_filepath` : str
        - `index_col` : int or sequence
        - `na_values` : list-like
        - `parse_dates` : bool
        
    :return: A pandas.DataFrame instance which represents the csv file.
    :rtype: pandas.DataFrame
    
    '''
    import pandas
    return pandas.read_csv(path(csv_filepath), index_col=index_col, na_values=na_values, parse_dates=parse_dates)


def pandasDataframe2csv(dataframe, csv_filepath, na_rep='', index=True, index_label=None):
    '''Write a DataFrame into a comma-separated values (csv) file.
    
    :Parameters:
        - `dataframe` : The DataFrame to write.
        - `csv_filepath` : The file path where the Dataframe is written.
        - `na_rep` : Missing data replacement.
        - `index_label` : Column label for index column(s) if desired. If None is given, and header and index are True, then the index names are used. A sequence should be given if the DataFrame uses MultiIndex.
        
    :Types:
        - `dataframe` : pandas.DataFrame
        - `csv_filepath` : str
        - `na_rep` : str
        - `index_label` : str or sequence

    :return: The file path where the Dataframe is written.
    :rtype: str
    
    '''
    dataframe.to_csv(path(csv_filepath), na_rep=na_rep, index=index, index_label=index_label)
    return csv_filepath   


from openalea.core import Node

class select_data_file(Node):
    '''Given a directory path and the pattern of data filenames, this node \
permits to select data files. For example, if the user chooses 'axeT*.csv' \
for the "axe" tables pattern, then the user can choose any axe table which \
belongs to the given directory and which satisfies the pattern 'axeT*.csv'. 
'''
    def __init__(self, inputs, outputs):
        Node.__init__(self, inputs, outputs)
        self._concatenated_paths = None

    def __call__(self, inputs):
        directory_path = inputs[0]
        nb_of_inputs = self.get_nb_input()
        nb_of_data = (nb_of_inputs - 1) / 2
        filenames = inputs[nb_of_data+1:nb_of_inputs]
        self._concatenated_paths = []
        if directory_path is not None:
            directory_path = path(directory_path)
            for filename in filenames:
                if filename is not None and filename != '':
                    concatenated_path = directory_path/filename
                else:
                    concatenated_path = None
                self._concatenated_paths.append(concatenated_path)
            return self._concatenated_paths
        else:
            return self._concatenated_paths


from PyQt4.QtGui import *
from PyQt4.QtCore import *

from openalea.core.observer import lock_notify
from openalea.visualea.node_widget import NodeWidget

class SelectDataFile(NodeWidget, QDialog):
    def __init__(self, node, parent):

        QDialog.__init__(self, parent)
        NodeWidget.__init__(self, node)

        self.gridlayout = QGridLayout(self)
        self.gridlayout.setMargin(3)
        self.gridlayout.setSpacing(5)

        self.directory_lineedit_label = QLabel('1. Set the directory', self)
        self.gridlayout.addWidget(self.directory_lineedit_label, 0, 0)
        
        self.directory_lineedit = QLineEdit(self)
        self.gridlayout.addWidget(self.directory_lineedit, 0, 1, 1, 3)
        self.connect(self.directory_lineedit, 
                     SIGNAL("textChanged(QString)"), 
                     self.directory_changed)

        self.selectdir_pushbutton = QPushButton('...', self)
        self.gridlayout.addWidget(self.selectdir_pushbutton, 0, 4, 1, 1)
        self.connect(self.selectdir_pushbutton, 
                     SIGNAL("clicked(bool)"), 
                     self.selectdir)

        self.selectdatafiles_label = QLabel('2. Set the patterns and select the \
data files', self)
        self.gridlayout.addWidget(self.selectdatafiles_label, 1, 0, 1, 4)

        self.nb_of_inputs = node.get_nb_input()
        self.nb_of_data = (self.nb_of_inputs - 1) / 2
        
        self.method_mapping = {}
        self.data_mapping = {}
        
        for i in range(self.nb_of_data):
            pattern_input_number = i + 1
            filename_input_number = pattern_input_number + self.nb_of_data
            pattern_input_port = node.get_input_port(pattern_input_number)
            pattern_input_port_name = pattern_input_port['name']
            if pattern_input_port_name.count('_') != 0:
                data_name = pattern_input_port_name.split('_')[0]
            else:
                data_name = str(pattern_input_number)
            self.method_mapping[pattern_input_number] = {}
            def data_pattern_changed(self, pattern_input_number):
                data_pattern = self.data_mapping[pattern_input_number]['data_pattern_lineedit'].text()
                self.node.set_input(pattern_input_number, str(data_pattern))
            self.method_mapping[pattern_input_number]['data_pattern_changed'] = data_pattern_changed
            def data_filename_changed(self, pattern_input_number):
                filename_input_number = pattern_input_number + self.nb_of_data
                data_filename = self.data_mapping[pattern_input_number]['data_filenames_combobox'].currentText()
                self.node.set_input(filename_input_number, str(data_filename))
            self.method_mapping[pattern_input_number]['data_filename_changed'] = data_filename_changed
            @lock_notify
            def update_data_pattern_input(self, pattern_input_number):
                data_pattern = self.node.get_input(pattern_input_number)
                if self.updating or data_pattern is None: return
                self.updating = True
                self.data_mapping[pattern_input_number]['data_pattern_lineedit'].setText(data_pattern)
                self.updating = False
            self.method_mapping[pattern_input_number]['update_data_pattern_input'] = update_data_pattern_input
            @lock_notify
            def update_data_filenames_combobox(self, pattern_input_number):
                directory = self.node.get_input(0)
                data_pattern = self.node.get_input(pattern_input_number)
                if directory is None or data_pattern is None: return
                directory = path(directory)
                if self.updating: return
                self.updating = True
                self.data_mapping[pattern_input_number]['data_filenames_combobox'].clear()
                self.data_mapping[pattern_input_number]['filename2filepath_mapping'].clear()
                paths = directory.glob(data_pattern)
                filenames = [a_path.name for a_path in paths]        
                self.data_mapping[pattern_input_number]['filename2filepath_mapping']\
                    .update(dict(zip(filenames, paths)))
                self.data_mapping[pattern_input_number]['data_filenames_combobox']\
                    .addItems(self.data_mapping[pattern_input_number]['filename2filepath_mapping'].keys())
                self.updating = False
            self.method_mapping[pattern_input_number]['update_data_filenames_combobox'] \
                = update_data_filenames_combobox    
            @lock_notify
            def update_data_filepath_output(self, pattern_input_number):
                data_output_filepath = None
                directory = self.node.get_input(0)
                filename_input_number = pattern_input_number + self.nb_of_data
                data_input_filename = self.node.get_input(filename_input_number)
                if directory is None or data_input_filename is None: return
                directory = path(directory)
                if self.updating: return
                self.updating = True
                data_input_filename_index = \
                    self.data_mapping[pattern_input_number]['data_filenames_combobox']\
                        .findText(data_input_filename)
                self.data_mapping[pattern_input_number]['data_filenames_combobox']\
                    .setCurrentIndex(data_input_filename_index)
                data_output_filepath = \
                    self.data_mapping[pattern_input_number]['filename2filepath_mapping']\
                        .get(data_input_filename)
                self.node.set_output(i, data_output_filepath)
                self.updating = False
            self.method_mapping[pattern_input_number]['update_data_filepath_output'] \
                = update_data_filepath_output
 
            self.data_mapping[pattern_input_number] = {}
            self.data_mapping[pattern_input_number]['data_label'] = QLabel(data_name, self)
            self.gridlayout.addWidget(self.data_mapping[pattern_input_number]['data_label'], i+2, 0)
            self.data_mapping[pattern_input_number]['data_pattern_lineedit'] = QLineEdit(self)
            self.gridlayout.addWidget(self.data_mapping[pattern_input_number]['data_pattern_lineedit'], i+2, 1, 1, 2)
            self.data_pattern_lineedit_signal_mapper = QSignalMapper(self)
            self.connect(self.data_mapping[pattern_input_number]['data_pattern_lineedit'], 
                         SIGNAL("textChanged(QString)"), 
                         self.data_pattern_lineedit_signal_mapper,
                         SLOT("map()"))
            self.data_pattern_lineedit_signal_mapper.setMapping(self.data_mapping[pattern_input_number]['data_pattern_lineedit'], 
                                          pattern_input_number);
            self.connect(self.data_pattern_lineedit_signal_mapper, 
                         SIGNAL("mapped(int)"), 
                         self.map_to_data_pattern_method)
            self.data_mapping[pattern_input_number]['data_filenames_combobox'] = QComboBox(self)
            self.gridlayout.addWidget(self.data_mapping[pattern_input_number]['data_filenames_combobox'], i+2, 3, 1, 2)
            self.data_filenames_combobox_signal_mapper = QSignalMapper(self)
            self.connect(self.data_mapping[pattern_input_number]['data_filenames_combobox'], 
                         SIGNAL("currentIndexChanged(QString)"), 
                         self.data_filenames_combobox_signal_mapper,
                         SLOT("map()"))
            self.data_filenames_combobox_signal_mapper.setMapping(self.data_mapping[pattern_input_number]['data_filenames_combobox'], 
                                          pattern_input_number);
            self.connect(self.data_filenames_combobox_signal_mapper, 
                         SIGNAL("mapped(int)"), 
                         self.map_to_data_filenames_method)
            self.data_mapping[pattern_input_number]['filename2filepath_mapping'] = {}
        
        self.setWindowTitle("SelectDataFile")
        self.setGeometry(250, 200, 350, 550)
        self.updating = False
        
        self.last_selected_directory = str(path.getcwd())

        self.notify(node, ("input_modified", 'all'))
        self.notify(node, ("caption_modified", node.get_caption()))
    
    def map_to_data_pattern_method(self, pattern_input_number):
        self.method_mapping[pattern_input_number]['data_pattern_changed'](self, pattern_input_number)
        
    def map_to_data_filenames_method(self, pattern_input_number):
        self.method_mapping[pattern_input_number]['data_filename_changed'](self, pattern_input_number)


    def notify(self, sender, event):
        ''' Update the widgets according to the notification sent by the node ''' 
        
        if event[0] == 'caption_modified':
            self.window().setWindowTitle(event[1])

        if(event[0] == "input_modified"):
            if event[1] == 0:
                self.update_directory_input()
                for pattern_number, methods in self.method_mapping.iteritems():
                    methods['update_data_filenames_combobox'](self, pattern_number)
                    methods['update_data_filepath_output'](self, pattern_number)
            elif event[1] <= self.nb_of_data:
                self.method_mapping[event[1]]['update_data_pattern_input'](self, event[1])
                self.method_mapping[event[1]]['update_data_filenames_combobox'](self, event[1])
                self.method_mapping[event[1]]['update_data_filepath_output'](self, event[1])
            elif event[1] < self.nb_of_inputs:
                pattern_number = event[1] - self.nb_of_data
                self.method_mapping[pattern_number]['update_data_filepath_output'](self, pattern_number)
            elif event[1] == 'all':
                self.update_directory_input()
                for pattern_number, methods in self.method_mapping.iteritems():
                    methods['update_data_pattern_input'](self, pattern_number)
                    methods['update_data_filenames_combobox'](self, pattern_number)
                    methods['update_data_filepath_output'](self, pattern_number)

    @lock_notify
    def update_directory_input(self):
        directory = self.node.get_input(0)
        if self.updating or directory is None: return
        self.updating = True
        self.directory_lineedit.setText(directory)
        self.updating = False

    def directory_changed(self, directory):
        self.node.set_input(0, str(directory))
  
    def selectdir(self, checked):
        dirname = QFileDialog.getExistingDirectory(self, 
                                                   'Select an existing data directory',
                                                   self.last_selected_directory,
                                                   QFileDialog.ShowDirsOnly)
        self.last_selected_directory = str(dirname)
        self.node.set_input(0, self.last_selected_directory)      

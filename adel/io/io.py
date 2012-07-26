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


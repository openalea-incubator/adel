import pandas
from openalea.core.path import path

def csv2dataframe(csv_filepath):
    '''Create a pandas.DataFrame from a csv file
    
    :Parameters:
        - `csv_filepath` : The filepath of the csv file to import.
        
    :Types:
        - `csv_filepath` : str
        
    :return: A pandas.DataFrame instance which represents the csv file.
    :rtype: pandas.DataFrame
    
    '''
    return pandas.read_csv(path(csv_filepath), index_col=0, na_values=['NA'], parse_dates=True)
    

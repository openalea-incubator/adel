import pandas
from openalea.core.path import path

def dataframe2csv(dataframe, csv_filepath):
    '''Export a pandas.DataFrame to a csv file.
    
    :Parameters:
        - `dataframe` : The pandas.DataFrame instance to export.
        - `csv_filepath` : The csv filepath where the pandas.DataFrame instance is exported.
        
    :Types:
        - `dataframe` : pandas.DataFrame
        - `csv_filepath` : str
    
    '''
    dataframe.to_csv(path(csv_filepath), na_rep='NA', index_label='date')
    return csv_filepath   

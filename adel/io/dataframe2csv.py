import pandas
from openalea.core.path import path

def dataframe2csv(dataframe, csv_filepath, na_rep='', index=True, index_label=None):
    '''Write a DataFrame to a comma-separated values (csv) file.
    
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

import pandas
from openalea.core.path import path

def csv2dataframe(csv_filepath, index_col=None, na_values=None, parse_dates=False):
    '''Read CSV (comma-separated) file into DataFrame.
    
    :Parameters:
        - `csv_filepath` : The filepath of the csv file to import.
	- `index_col` : Column to use as the row labels of the DataFrame. If a sequence is given, a MultiIndex is used.
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
    return pandas.read_csv(path(csv_filepath), index_col=index_col, na_values=na_values, parse_dates=parse_dates)
    

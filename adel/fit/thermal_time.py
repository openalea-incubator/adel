import pandas
import alinea.echap.ADASWheatReconstruction.thermal_time as recons_thermal_time

def thermal_time(requested_dates, emergence_date, input_data, data_type, thermal_time_increment_method, latitude):
    '''Calculates the thermal time increment. 
  
    :Parameters:
        - `requested_dates` : The dates for which the thermal time increment has to be calculated. These dates must belong to input_data, and be >= to emergence_date.
        - `emergence_date` : The emergence date. The interval [emergence date, emergence date + 1 year[ must be included in input_data. emergence_datetime must be <= than requested_dates
        - `input_data` : The input data. Expect a pandas.DataFrame object. The dataframe index must contain dates.
        - `data_type` : The type of the input data. Can be either 'daily' or 'hourly'.
        - `thermal_time_increment_method` : The method used for thermal time increment calculation. Can be either 'dsT' or 'dsTc' (compensated).
        - `latitude` : The latitude where the data have been measured. Used to convert 'daily' data to 'hourly' data.
         
    :Types:
        - `requested_dates` : pandas.DataFrame
        - `emergence_date` : dict
        - `input_data` : pandas.DataFrame
        - `data_type` : str enumerate
        - `thermal_time_increment_method` : str 
        - `latitude` : float
  
    :return: The thermal time increment calculated.
    :rtype: pandas.DataFrame
    
    '''    
    emergence_datetime = pandas.datetime(emergence_date['year'], emergence_date['month'], emergence_date['day'], emergence_date['hour'], emergence_date['minute'], emergence_date['second'])
    # check that input_data contains 'emergence_datetime' and 'emergence_datetime + 1 year'
    if emergence_datetime not in input_data.index:
        raise Exception(' '.join(['emergence_datetime', emergence_datetime, 'does not belong to input_data']))
    start = emergence_datetime
    end = pandas.datetime(emergence_datetime.year + 1, emergence_datetime.month, emergence_datetime.day, emergence_datetime.hour, emergence_datetime.minute, emergence_datetime.second)
    if end not in input_data.index:
        raise Exception(' '.join(['emergence_datetime + 1 year', end, 'does not belong to input_data']))
    
    def crit_year(index_i):
        # Permits to extract data from start to end
        if index_i >= start and index_i <= end:
            return True
        return False
    
    # extract data from start to end
    extracted_data = input_data.select(crit_year)
    
    if data_type == 'daily':
        # transform to hourly data
        hourly_data = recons_thermal_time.parton_logan(extracted_data, latitude)
    else: # hourly data
        # calculate the mean for each (Tmin,Tmax) tuple
        Tair_mean = (extracted_data['Tmin'] + extracted_data['Tmax']) / 2.0
        hourly_data = pandas.DataFrame(Tair_mean, index=extracted_data.index, columns=['Tair'])
        
    # check that input_data contains requested_dates
    # and that emergence_datetime is <= to requested_dates
    for requested_date in requested_dates.index:
        if requested_date not in input_data.index:
            raise Exception(' '.join(['requested_date', requested_date, 'does not belong to input_data']))
        if requested_date < emergence_datetime:
            raise Exception(' '.join(['requested_date', requested_date, 'is < than emergence_datetime', emergence_datetime]))
        
    # calculate thermal time increment for all dates with the given thermal_time_increment_method
    available_methods = {'dsT': recons_thermal_time.dsT, 
                         'dsTc': recons_thermal_time.dsTc}
    thermal_time_increment_array = available_methods[thermal_time_increment_method](hourly_data['Tair'].values)
    thermal_time_increment_dataframe = pandas.DataFrame(thermal_time_increment_array, index=hourly_data.index, columns=[thermal_time_increment_method])
    
    def crit_dates(index_i):
        # Permits to select requested_dates
        if index_i in requested_dates.index:
            return True
        return False
    
    # select data for requested_dates and return them
    return thermal_time_increment_dataframe.select(crit_dates)
    
    
    
    
def ScatteringOptions(No_Multiple_Scaterring,OptDict):
    '''
    handles multiple scaterring options for Caribu
    '''
    Sleep=No_Multiple_Scaterring

    # return outputs
    return (OptDict['Nz'],OptDict['Zmax'],Sleep,No_Multiple_Scaterring,OptDict['SphereDiameter'],OptDict['keepFF'])

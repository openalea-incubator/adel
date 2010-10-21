def GetLight(SkyObject):
    '''    Get the incident light and the coordinates of sky sectors 
    '''
    
    # write the node code here.
    LightString=""
    for i in range(len(SkyObject.sky)):
        LightString+=((str(SkyObject.sky[i][0])+' '+str(SkyObject.sky[i][1])+' '+str(SkyObject.sky[i][2])+' '+str(SkyObject.sky[i][3])+'\n'))

    # return outputs
    return LightString

class get_meteo_dat(object):
    """  Doc... """ 

    def __init__(self):
        pass


    def __call__(self, Tab_Rg, group, indice):
        DOY, HU = int(Tab_Rg[group][0][indice+1]), int(Tab_Rg[group][1][indice+1]) 
        Rd = float(Tab_Rg[group][2][indice+1])*float(Tab_Rg[group][3][indice+1])*3600/1000000 #integre W.m-2 sur pas de temps (1h)
        Rsun = float(Tab_Rg[group][2][indice+1])*(1-float(Tab_Rg[group][3][indice+1]))*3600/1000000 #integre W.m-2 sur pas de temps (1h)
        return  Rd, Rsun, DOY, HU #Rd, Rsun en MJ.m-2

############ G Louarn - adaptation de spitters.c EGC grignon
from math import *


class spitters_horaire(object):
    """  Doc... """ 

    def __init__(self):
        pass


    def __call__(self, Tab_Rg, latitude):
        """ calcule RdRg et ajoute le resultat dans dans Tab_Rg """
        for group in range(len(Tab_Rg)):
            Tab_Rg[group].append(["RdRg"])
            for i in range(1,len(Tab_Rg[group][0])):
                Rg,DOY,heureTU = float(Tab_Rg[group][2][i]),int(Tab_Rg[group][0][i]),float(Tab_Rg[group][1][i])
                frac = self.RdRsH (Rg,DOY,heureTU,latitude)
                Tab_Rg[group][4].append(str(frac))
        
            Tab_Rg[group][3],Tab_Rg[group][4] = Tab_Rg[group][4],Tab_Rg[group][3]

        return (Tab_Rg,)
    
    def DecliSun (self,DOY):
        """ Declinaison (rad) du soleil en fonction du jour de l'annee """
        alpha=2*3.14*(DOY-1)/365
        return (0.006918-0.399912*cos(alpha)+0.070257*sin(alpha))
    
    def DayLength (self,latitude,decli):
        """ photoperiode en fonction de latitude (degre) et declinaison du soleil (rad) """
        lat=radians(latitude)
        d=acos(-tan(decli)*tan(lat))
        if d<0:
            d=d+3.14
        
        return 2*d
    
    def extra (self,Rg,DOY,heureTU,latitude):
        """ rayonnement extraterrestre horarire """
        hrad=2*3.14/24*(heureTU-12)
        lat=radians(latitude)
        dec=self.DecliSun (DOY)
        costheta=sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(hrad)
        Io=1370*(1+0.033*cos(2*3.14*(DOY-4)/366))#eclairement (w/m2) a la limitte de l'atmosphere dans un plan perpendiculaire aux rayons du soleil, fonction du jour
        So=Io*costheta #eclairement dans un plan parallele a la surface du sol
        return So
    
    def RdRsH (self,Rg,DOY,heureTU,latitude):
        """ fraction diffus/Global en fonction du rapport Global(Sgd)/Extraterrestre(Sod)- pas de temps horaire """
        hrad=2*3.14/24*(heureTU-12)
        lat=radians(latitude)
        dec=self.DecliSun (DOY)
        costheta=sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(hrad)
        Io=1370*(1+0.033*cos(2*3.14*(DOY-4)/366))#eclairement (w/m2) a la limitte de l'atmosphere dans un plan perpendiculaire aux rayons du soleil, fonction du jour
        So=Io*costheta #eclairement dans un plan parallele a la surface du sol
        RsRso=Rg/So
        R=0.847-1.61*costheta+1.04*costheta*costheta
        K=(1.47-R)/1.66
        
        if (RsRso<=0.22) :
            return(1)
        else:
            if (RsRso<=0.35) :
                return(1-6.4*(RsRso-0.22)*(RsRso-0.22))
            else:
                if (RsRso<=K) :
                    return(1.47-1.66*RsRso)
                else:
                    return(R)
            



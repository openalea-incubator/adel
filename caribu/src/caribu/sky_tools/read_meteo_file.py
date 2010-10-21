import IOtable

class read_meteo_file(object):
    """  Doc... """ 

    def __init__(self):
        pass


    def __call__(self, filename):
        f = file(filename, 'r')
        Tab = IOtable.table_csv_str(f) 
        f.close()
        return (self.join_group(Tab),)#lourdeur pour renvoyer tables

    def join_group(self, table):
        """ strucure donnees meteo en colones & selon groupes specifies """
        res=[]
        header=table[0]
        group=-1
        for i in range(1,len(table)):
            if int(table[i][-1])!=group:#cas1 : changement de groupe
                if group != -1: #excepte pour 1ere ligne, ajoute groupe a res
                    res.append(g)
    
                if len(header)==4:#nouveau groupe
                    g = [[],[],[],[]]#fichier meteo sans RdRg
                else:
                    g= [[],[],[],[],[]]#avec RdRg
    
                for j in range(len(g)):
                    g[j].append(header[j])
                
                group = int(table[i][-1])#actualise num groupe
                for j in range(len(g)):
                    g[j].append(table[i][j])    
            else:#cas2 : pas de changement de groupe
                for j in range(len(g)):
                    g[j].append(table[i][j])
    
        res.append(g)#ajoute dernier groupe
        return res



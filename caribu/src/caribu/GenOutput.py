def optical_species(label):
    return int(label[:-11])

def plant_id(label):
    return int(label[-11:-6])

def transparency(label):
    return int(bool(leaf_id(label)))

def leaf_id(label):
    return int(label[-6:-3])

def elt_id(label):
    return int(label[-3:])



def GenOutput(etri, eabs):
    '''    Merge canestra string outputs (Etri,Eabs) into a dictionary of float vectors
    '''
    EtriDic = None;
    
    # write the node code here.

    opt=[]
    opak=[]
    plt=[]
    elt=[]
    area=[]
    eabsm2=[]
    eiinf=[]
    eisup=[]

    eabs=eabs.split('\n')
    #filter empty lines
    eabs=filter(lambda(x): len(x) >0,eabs)
    lg=etri.split('\n')

    
    for l in lg :
        l=l.strip()
        if not l or l.startswith('#'):
            continue
        col=l.split()
        #filter empty lines
        if(not col) : continue
        label=col[0]
        if len(label) < 11:
            label = (12-len(label))*'0'+label    
        opt.append(str(optical_species(label)))
        opak.append(str(leaf_id(label)))
        plt.append(str(plant_id(label)))
        elt.append(str(elt_id(label)))
        area.append(col[1])
        eabsm2.append(col[2])
        eisup.append(col[3])
        eiinf.append(col[4])
        
    EtriDic={'Eabs':eabs,'Area':area,'Eabsm2':eabsm2,'EiInf':eiinf,'EiSup':eisup,'Opt':opt,'Opak':opak,'Plt':plt,'Elt':elt}
    #conversion en float
    for k in EtriDic.iterkeys():
        EtriDic[k] = map(float,EtriDic[k])

    # return outputs
    return EtriDic


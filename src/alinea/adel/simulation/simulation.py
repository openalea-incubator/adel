import os
from openalea.core.path import path

from alinea.adel.AdelR import RunAdel,genString

def simulation(lsystem_filename, nb_plants):
    
    # write AleaPars.h
    fn = path(lsystem_filename)
    dir = fn.dirname()

    param_fn = dir / "AleaPars.h"
    s = '\n#define NB_PLANTS %d\n'%nb_plants
    f = open(param_fn, "w")
    f.write(s)
    f.close()

    view_file = 'view.v'
    leaf_file = 'leaf.a'
    color_file = 'color.map'
    result_file = 'result.str'
    cmd = "cpfg -g -homo -str %s %s %s %s -m %s" % (result_file, fn.basename(), view_file, leaf_file, color_file)
    print(cmd)
    curdir = os.getcwd()
    print(dir)
    os.chdir(dir)
    os.system(cmd)
    os.chdir(curdir)
    result_file = dir / result_file
    if not result_file.exists():
        result_file = dir/'100plants'/'bid050.str'
        
    f = open(result_file)
    s = f.read()
    s = s.replace('\n','')
    f.close()
    return s

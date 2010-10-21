import time
import os, popen2
import string
from SailScene import SailScene

def getNumLayer(file):
    """ return layer number of triangles """
    
    numl=[]
    if os.path.isfile(file):
        fin = open(file)
        lines = fin.readlines()
        fin.close()
        
        for l in lines:
            numl.append(int(l.split(' ')[1]))
    return numl
#implementer la prise en charge de la bande de longueur d'onde
  
def S2v(scene,nlayer,zmax,sleep):
    """  Transform a 3D CaribuScene into 1D Sail Scene  """

    log = "\nS2v :Last call => " + time.asctime() + "\n"
   
    if sleep:
        log += " s2V is sleeping...SailScene has not been generated\n"
        print log
        return None,None,log

    wd = os.getcwd()    
    tmpdir=os.tempnam()
    os.mkdir(tmpdir)
    os.chdir(tmpdir)
    
    print 'tmpdir ', tmpdir
	
    canfile="scene.can"
    patternfile="scene.8"
    optfile="scene.opt"
   
    scene.writeCan(canfile)
    scene.writePattern(patternfile)
    scene.writeOptical(optfile)
        
    cmd = ' '.join(["s2v",canfile, str(nlayer), str(zmax),patternfile,optfile[0:-4]])

    print cmd
    
    fout, fin = popen2.popen2(cmd)
    stderr = ''.join(fout.readlines())
    print stderr + '\n S2v : Success !!'
    log +=  cmd + "\n" + stderr
    fin.close()
    fout.close()
    
    sc = SailScene()
    sc.setScene("cropchar","leafarea","scene.spec")
    sc.setLight(scene.sources)
    numLayer=getNumLayer("s2v.area")
    if 0 in numLayer:
        print 'Warning ! Zmax too low to place all scene triangles !'
        log+='s2V : Warning ! Zmax too low to place all scene triangles !'
    
    for root, dirs, files in os.walk(tmpdir, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    
        
    os.chdir(wd)
    os.rmdir(tmpdir)
        
    return (sc,numLayer,log)

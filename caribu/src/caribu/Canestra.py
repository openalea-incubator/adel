import os
import time
import string
import copy

def Canestra(CaribuScene,
             meanFluxes,
             option_NoScattering,
             sphereDiameter,
             option_keepFF          
             ):
    
    """  Nested radiosity illumination of a 3D CaribuScene """

    log = "\nCanestra :Last call => " + time.asctime() + "\n"
    
           
    wd = os.getcwd()    
    tmpdir=os.tempnam()
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    effargs=['-A'] #generate Etri.vec

    CaribuScene.writeCan("scene.can")
    CaribuScene.writeLight("scene.light")
    CaribuScene.writeOptical("scene.opt")
    effargs+=['-M', "scene.can",'-p','scene.opt','-l', "scene.light"]
    if CaribuScene.hasPattern :
        CaribuScene.writePattern("scene.8")
        effargs+=['-8','scene.8']
    if CaribuScene.hasFF :
        log += "Using FFmatrix of the scene ...\n"
        CaribuScene.writeFF("FF.mat")
        effargs+=["-t", '"%s"'%tmpdir,'-w','FF.mat']
            

    if option_NoScattering:

        log += " No Scattering asked=> enter simple projection mode\n"
        effargs+=["-1"]
	
    else :
        
        ff = open("sail.env","w")
        ff.write(meanFluxes)
        ff.close()
        effargs+=["-e", "sail.env"]       
        effargs+=["-d",str(sphereDiameter)]
        if option_keepFF:
            effargs+=["-t", '"%s"'%tmpdir,'-f','FF.mat']

 
    cmd = "canestrad " + " ".join(effargs) 
    print cmd
    fin, fout = os.popen2(cmd)
    stderr = ''.join(fout.readlines())
    print stderr + '\n Canestra : Success !!'
    log +=  cmd + "\n" + stderr
    fin.close(); fout.close()

    ff=open("Etri.vec")
    etri=ff.read()
    ff.close()
    
    ff=open("Eabs.vec")
    eabs=ff.read()
    ff.close()

    SceneOut = copy.copy(CaribuScene)
    if option_keepFF:
        SceneOut.setFF("FF.mat")
    ##print tmpdir
    ff=open("scene.can")
    SceneOut.setCan(ff.read())
    ff.close()
    
    for root, dirs, files in os.walk(tmpdir, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))

    os.chdir(wd)
    os.rmdir(tmpdir)

    return(SceneOut,etri,eabs,log)

import os
import time
import string


def MCSail(SailScene,sleep):
    """  Compute light fluxes on a layered canopy  """ 
    
    log = "\nMCSail :Last call => " + time.asctime() + "\n"

    if sleep:
        log += " MCSail is sleeping...MeanFluxes have not been generated\n"
        print log
        return None,log
    
    wd = os.getcwd()    
    tmpdir=os.tempnam()
    os.mkdir(tmpdir)
    os.chdir(tmpdir)
    
    SailScene.writeCropChar("cropchar")
    SailScene.writeLeafArea("leafarea")
    SailScene.writeSpectral("spectral")
    SailScene.writeLight("sailIn.light")
    
    effargs = ["sailIn.light"]
    cmd = "mcsail " + ' '.join(effargs)
    print cmd

    fin, fout = os.popen2(cmd)
    stderr = ''.join(fout.readlines())
    log +=  cmd + "\n" + stderr
    fin.close(); fout.close()
    if os.path.isfile("mlsail.env"):
        print stderr + '\n MCSail : Success !!'
        ff=open("mlsail.env")
        fluxes=ff.read()
        ff.close()
    else:
        fluxes=None
        print stderr + '\n MCSail : Failure, sorry !'

    for root, dirs, files in os.walk(tmpdir, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    
    os.chdir(wd)
    os.rmdir(tmpdir)

    return(fluxes,log)

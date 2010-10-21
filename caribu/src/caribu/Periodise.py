import os
import time
import string
import copy

def dQuote(fn):
    return '"%s"'%fn

def Periodise(scene):
    """  Fit a scene within its pattern  """ 
            
    log = "\nPeriodise :Last call => " + time.asctime() + "\n"
    
    if not scene.hasPattern :
        log += " periodise : no pattern => Leaving scene unchanged\n"
        print log
        periodicscene = copy.copy(scene)
        return periodicscene,log

    tmpdir=os.tempnam()
    os.mkdir(tmpdir)
    canfile=os.path.join(tmpdir,"scene.can")
    patternfile=os.path.join(tmpdir,"scene.8")
    outfile=os.path.join(tmpdir,"scene.out")

    scene.writeCan(canfile)
    scene.writePattern(patternfile)
        
    cmd = ' '.join(["periodise","-m", dQuote(canfile), "-8", dQuote(patternfile),"-o",dQuote(outfile)])

    print cmd
    fin, fout = os.popen2(cmd)
    stderr = ''.join(fout.readlines())
    
    try:
        fin.close()
        fout.close()
    except:
        print stderr + '\n Periodise : Failure !!!'
    else:
        print stderr + '\n Periodise : Success !!'

    log +=  cmd + "\n" + stderr
    periodicscene = copy.copy(scene)	

    try:
        ff = open(outfile)

    except IOError:
        print 'erreur execution periodise (empty scene returned) : see log'
        periodicscene.setCan('#periodise error : no scene generated\n')
        for f in [canfile,patternfile] :
            os.remove(f)
    else:
        periodicscene.setCan(ff.read())
        ff.close()    
        for f in [canfile,patternfile,outfile] :
            os.remove(f)

    os.rmdir(tmpdir)
        
    return (periodicscene,log)

""" Chaine les appel de s2v, mcsail et canestra, comme le fait cpfg via caribu
    Syntaxe: caribu.csh Ds file.light file.can nz h file.8 file1.opt ... fileN.opt

  Exemple /home/tld/chelle/QHS/Calgary/dev/Canestra/MC-SAIL/Test 
  caribu.csh 0.2 sky2.light Tout.can 3 18 Tout.8 test nir par

  MC98
  C. Pradal - 2009
  MC09
  Contact: chelle@grinon.inra.fr
  INRA - INRIA - CIRAD
"""
import os
from subprocess import Popen,STDOUT, PIPE
import tempfile
from openalea.core.path import path

def _process(cmd, directory, out):
    """ 
    Run a process in a shell. 
    Return the outputs in a file or string.
    """
    #print ">> caribu.py: process(%s) called..."%(cmd)

    f = open(out,'w')
    p = Popen(cmd, shell=True, cwd=directory,
              stdin=PIPE, stdout=f, stderr=STDOUT)
    status = p.wait()

    f.close()
    
    #print "<<< caribu.py: process finished!"
    return status

def _safe_iter(obj, atomic_types = (basestring, int, float, complex)):
    """Equivalent to iter when obj is iterable and not defined as atomic.
    If obj is defined atomic or found to be not iterable, returns iter((obj,)).
    safe_iter(None) returns an empty iterator"""
    if not isinstance(obj, atomic_types):
        try:
            return iter(obj)
        except TypeError:
            pass
    return iter((obj,) * (obj is not None))

def _abrev(fnc,maxlg=1):
    """
    abreviate a text string containing a path or a file content to the first maxlg lines,
    addind '...' when the number of libnes is greater than maxlg
    """
    if fnc is None or os.path.exists(fnc):
        return(str(fnc))
    lines = fnc.splitlines()
    if maxlg <= 1:
        return(lines[0]+' ... ')
    elif len(lines) <= maxlg:
        return(str(fnc))
    else:
        return('\n'+'\n'.join(lines[0:maxlg]) + "\n...")


class CaribuError(Exception): pass
class CaribuOptionError(CaribuError): pass
class CaribuRunError(CaribuError): pass
    
class Caribu(object):
    
    def __init__(self,
                 canfile = None, 
                 skyfile = None, 
                 optfiles = None,
                 patternfile = None,
                 optnames = None,
                 direct = True,
                 nb_layers = None,
                 can_height = None,
                 sphere_diameter = -1,
                 debug = False,
                 resdir = "./Run",
                 resfile = None
                 ):
        """
        Class fo Nested radiosity illumination on a 3D scene.
        
        canfile: file '.can' (or file content) representing 3d scene
        skyfile: file/file content containing all the light description
        optfiles: list of files/files contents defining optical property
        optnames: list of name to be used as keys for output dict (if None use the name of the opt files or the generic names band0,band1 if optfiles are given as content)
        patternfile: file/file content that defines a domain to till the scene. Consider a toric canopy if pattern is not None
        direct: consider only direct projection
        nb_layers: number of layers to be consider for the scene
        can_height: height of the can scene
        sphere_diameter: used for the radiosity
        debug : print messages and prevent removal of tempdir
        resdir : store caribu results as files in resdir if resdir is not None, store nothing otherwise
        resfile : store caribu output dictionary in file resfile (with pickle) if resfile is not None, store nothing otherwise
        """
        
        print "\n >>>> Caribu.__init__ starts...\n"
        #debug mode
        self.my_dbg = debug
        #print "my_dbg = ",   self.my_dbg
        #tempdir (initialised to allow testing of  existence in del)
        self.tempdir = path('')

        # Input files
        self.scene = canfile
        self.sky = skyfile
        self.opticals = optfiles

        # User options
        self.optnames = optnames
        self.resdir = resdir
        self.resfile = resfile
         # consider a toric canopy
        self.pattern = patternfile
        self.infinity = None
        
        self.form_factor = None
        self.direct = direct           # direct light only
        self.nb_layers = nb_layers     # grid turbid medium
        self.can_height = can_height       # height of the canopy
        self.sphere_diameter = sphere_diameter # parameters for nested radiosity

        self.canestra_name = "canestrad"
        self.sail_name = "mcsail"
        self.periodise_name = "periodise"
        self.s2v_name = "s2v"
        self.ready = True

        print "\n <<<< Caribu.__init__ ends...\n"

    def __del__(self):
        print "Caribu.__del__ called !"
        if not self.my_dbg and self.tempdir.exists():
            print 'Remove tempfile %s'%self.tempdir
            self.tempdir.rmtree()

    def __str__(self):         
        s = """
            scene %s
            sky %s
            optnames %s
            opticals %s
            pattern %s
            infinity %s
            direct %s
            nb_layers %s
            can_height %s
            sphere_diameter %s
            form_factor %s
            ------------
            canestrad: %s
            mcsail: %s
            periodise: %s
            s2v: %s
        """%(_abrev(self.scene), _abrev(self.sky), ' '.join(map(str,_safe_iter(self.optnames))),''.join(map(_abrev,_safe_iter(self.opticals))), self.pattern, self.infinity,self.direct, self.nb_layers, self.can_height, self.sphere_diameter, self.form_factor, self.canestra_name, self.sail_name,self.periodise_name,self.s2v_name)
        if self.my_dbg:
            sopt = """
            -----------
            debug on
            tempdir %s
            resdir %s
            resfile %s
            """%(self.tempdir,self.resdir,self.resfile)
            s += sopt
        return(s)
        
    def show(self, titre="############"):
        print "\n>>> Caribu state in ", titre
        print(self)
        print "<<<<\n\n"

    def init(self):
        if self.scene == None or self.sky==None or self.opticals== None or self.opticals==[] :
            raise CaribuOptionError(">>> The Caribu has not been fully initialized: scene, sky, and opticals have to be defined\n     =>  Caribu can not be run... - MC09")
        
        if self.pattern == None:
            self.infinity = False
            print ">>> pattern not specified => computations done on a NON toric scene"
        else:
            self.infinity = True
            
        # print "infty, pattern", self.infinity, self.pattern    
            
        if not self.direct: 
            if  self.sphere_diameter < 0 :
                # Compute classic radioity without toric scene
                self.pattern = None
                self.infinity = False
                print ">>> sphere_diameter < 0 => classic radiosity called on a NON toric scene"
            else:   ## diameter >=0 
                # consider a toric canopy
                if not self.infinity :   
                    raise CaribuOptionError(">>> incompatible options for nested radiosity: pattern==None (not infinity) &&  sphere_diameter >= 0 ")
                
        self.form_factor=True        
        # self.canestra_1st = True # Boolean that indicates the first or not times, canestra is called thus form factors computed...

        #nrj is a dictionary of dictionary, each containing one simulation outputs. There will be as much dictionaries as optical files given as input 
        self.nrj = {}
        self.show("Caribu::init()")
        
        # Working directory        
        if self.my_dbg:
            print "I'm here ..."
            self.tempdir=path("./Run-tmp")
            if not self.tempdir.exists():
                self.tempdir.mkdir() 
        else:
            # build a temporary directory
            self.tempdir = path(tempfile.mkdtemp())

        # Result directory (if specified)
        if self.resdir is not None:
            self.resdir = path(self.resdir)
            if not self.resdir.exists():
                self.resdir.mkdir()

        # name of band to process (if not given)
        if self.optnames is None:
            #try to derive from filename or use generic name
            optn = []
            for i,opt in enumerate(_safe_iter(self.opticals)):
                if os.path.exists(opt):
                    name = str(path(path(opt).basename()).stripext())
                    optn.append(name)
                else :
                    optn.append('band%d'%(i))
            self.optnames = optn
    
        # Copy the files (or file content) in the tempdir
        self.copyfiles()
           

    def copyfiles(self):
        d = self.tempdir
        
        if os.path.exists(self.scene):
            fn = path(self.scene)
            fn.copy(d/fn.basename())
        else :
            fn = d/'scene.can'
            fn.write_text(self.scene)
        self.scene = path(fn.basename())

        if os.path.exists(self.sky):
            fn = path(self.sky)
            fn.copy(d/fn.basename())
        else :
            fn = d/'sky.light'
            fn.write_text(self.sky)
        self.sky = path(fn.basename())

        if self.infinity:
            if os.path.exists(self.pattern):
                fn = path(self.pattern)
                fn.copy(d/fn.basename())
            else : 
                fn = d/'pattern.8'
                fn.write_text(self.pattern)
            self.pattern = path(fn.basename())

        optn = map(lambda(x): x + '.opt',_safe_iter(self.optnames))
        try:
            for i,opt in enumerate(_safe_iter(self.opticals)):
        #safe_iter allows not to iterate along character composing the optfile name when only one optfile is given
                if os.path.exists(opt):
                #print opt
                    fn = path(opt)
                    fn.copy(d/optn[i])
                else:
                    fn = d/optn[i]
                    fn.write_text(opt)
            self.opticals = map(path,_safe_iter(optn))
        except IndexError:
            raise CaribuOptionError("Optnames list must be None or as long as optfiles list")

    def store_result(self,filename,band_name):
        """
        Add a new entry to the nrj dictionnary, using band_name as key and a dictionary build from filename as value.
        The dictionary build from filename is organised as follow:
            - doc : the first line of filename, that contains informations on the simulation
            - data : a dictionary of vectors, each containing a column of filename
            Columns are:
                - index (float): the polygon index
                - label (str): its can label
                - area (float): its area
                - Eabs,Ei_sup and Ei_inf (float): surfacic density (energy/s/m2) of, respectively, absorbed energy, irradiance on the adaxial side and irradiance on the abaxial side of polygons
        """   
            
        f = open(filename)
        doc = f.readline() # elimine la ligne de commentaire
        f.readline()  # elimine la ligne de commentaire
        idx=[]
        label=[]
        area=[]
        Eabs=[]
        Ei_sup=[]
        Ei_inf=[]
        for line in f:
            elements = line.split() # tu split ta chaine en fonction d'une string de separation: par defaut c'est ' ', '\t', mais tu peux faire split(','), ...
            floats = [float(el) for el in elements] # tu convertis toutes les strings du tableau elements en float. floats est un tableau
            idx.append( floats[0])
            lab = elements[1]
            if len(lab) < 11:
                lab = (12 - len(lab)) * '0' + lab
            label.append(lab)
            area.append(floats[2])
            Eabs.append(floats[3])
            Ei_sup.append(floats[4])
            Ei_inf.append(floats[5])
                 
        f.close()
        data={'index':idx, 'label':label, 'area':area,'Eabs':Eabs,'Ei_sup':Ei_sup, 'Ei_inf':Ei_inf}
        self.nrj[band_name] = {'doc':doc, 'data':data}

    def run(self):
        """
        The main Caribu program.
        1. Periodise: to convert the scene into an infinite one.
        2. s2v: Surface to volume based on the scene, the height of the canopy and optical prop.
        3. mcsail: mean fluxes in the canopy
        4. canestra: compute radiosity
        5. save output on disk if resfile specified
        """
        print "\n >>>> Caribu.run() starts...\n"
        self.init()
        self.periodise()
        self.s2v()
        self.radiative_transfer()
        if self.resfile is not None:
            import pickle 
            file = open(self.resfile, 'w') 
            pickle.dump(self.nrj, file)
            # To restore the value of the object to memory, load the object from the file. Assuming that pickle has not yet been imported for use, start by importing it:
            
            #import pickle 
            #file = open('caribu_run.obj', 'r') 
            #caribu_run = pickle.load(file)
            #x=caribu_run
            #print x['par']['data']['Eabs'][0]
        print "\n <<<< Caribu.run() ends...\n"


    def periodise(self):
        if self.infinity:
            d = self.tempdir
            name, ext = self.scene.splitext()
            outscene = name + '_8' + ext
            cmd = '%s -m %s -8 %s -o %s '%(self.periodise_name,self.scene, self.pattern, outscene)
            print ">>> periodise() : ",cmd
            status = _process(cmd, d, d/"periodise.log")
            if (d/outscene).exists():
                self.scene = outscene
            else:
                raise CaribuRunError("Periodise failed (no scene created).")

    def s2v(self):
        if self.infinity and not self.direct:
            d = self.tempdir
            wavelength = ' '.join([fn.stripext() for fn in self.opticals])
            cmd = "%s %s %d %f %s "%(self.s2v_name,self.scene, self.nb_layers, self.can_height, self.pattern) + wavelength
            print ">>> s2v() : ",cmd
            status = _process(cmd, d, d/"s2v.log")
            # Raise an exception if s2v crashed...
            leafarea= d/'leafarea'
            if not leafarea.exists():
                raise CaribuRunError(">>> s2v(): s2v failed (no leafarea created)")

    def radiative_transfer(self):
        # optics = [fn.stripext() for fn in self.opticals]
        # for opt in optics:
        for opt in self.opticals:
            print ">> radiative_transfer() : %s"%opt
            self.mcsail(opt)
            self.canestra(opt)

    def mcsail(self, opt):
        if self.infinity and not self.direct:
            d = self.tempdir
            optname, ext = path(opt.basename()).splitext()
            d = self.tempdir
            (d/optname+'.spec').copy(d/'spectral')
                
            cmd = "%s %s "%(self.sail_name,self.sky)
       
            print ">>> mcsail(): ",cmd
            logfile = "sail-%s.log"%(optname)
            logfile = d/logfile
            status = _process(cmd, d, logfile)
            
            mcsailenv = d/'mlsail.env'
            if mcsailenv.exists():
                mcsailenv.move(d/optname + '.env')
            else:
                raise CaribuRunError(">>> mcsail(): mcsail failed with %s (no env created)."%opt)

    def canestra(self, opt):
        """Fonction d'appel de l'executable canestrad, code C++ compilee de la radiosite mixte  - MC09"""
        #canestrad -M $Sc -8 $argv[6] -l $argv[2] -p $po.opt -e $po.env -s -r  $argv[1] -1
        d = self.tempdir        
        optname, ext = path(opt.basename()).splitext()
        
        str_pattern = str_direct = str_FF = str_diam = str_env = ""
        
        if self.infinity:
            str_pattern=" -8 %s "%(self.pattern)
            
        if self.direct:
            str_direct=" -1 "
        else:
            str_diam=" -d %s "%(self.sphere_diameter)

            if self.form_factor:
             #compute formfactor
                self.form_factor = False
                self.FF_name= tempfile.mktemp(prefix="", suffix="", dir="")
                str_FF=" -f %s "%(self.FF_name)
            else: 
                str_FF=" -w " + self.FF_name
            if self.sphere_diameter >= 0 :
                str_env=" -e %s.env "%(optname)
                
        cmd = "%s -M %s -l %s -p %s -A %s %s %s %s %s "%(self.canestra_name,self.scene, self.sky,  opt, str_pattern,str_direct,str_diam, str_FF, str_env)  
        print(">>> Canestrad(): %s"%(cmd))
        status = _process(cmd, self.tempdir,d/"nr.log")
        
        if  (d/path("nr.log")).exists():
            # copy log files
            fic = path("nr-" + optname + ".log")
            (d/"nr.log").move(d/fic)

        ficres = d/'Etri.vec0'
        if ficres.exists():
            self.store_result(ficres,str(optname))
            print optname
        
            if self.resdir is not None:
                # copy result files
                fdest = path(optname + ".vec")
                print fdest
                ficres.move(self.resdir/fdest)
        else:
            raise CaribuRunError(">>>  canestra has not finished properly => STOP")
        print ">>> caribu.py: Caribu::canestra(%s) finished !"%(opt)
        

def caribu_test(cas):
    d = path(__file__).dirname()
    scene = d/'data/filterT.can'
    print scene.splitext()
    sim=Caribu(canfile=d/'data/filterT.can', skyfile=d/'data/zenith.light', optfiles=[d/'data/par.opt', d/'data/nir.opt'])

    sim.infinity = False       # consider a toric canopy
    sim.direct = True           # direct light only
    sim.nb_layers = None
    sim.can_height = None
    sim.sphere_diameter = -1
    sim.pattern = None
    run=True
    
    if cas ==1:
        # Case 1: projection on non toric scene
        pass
    elif cas ==2:
        # Case 2: projection on  toric scene
        sim.pattern = d/"data/filter.8"
    elif cas == 3:
        # Case 3: classic radiosity on  non toric scene
        sim.direct = False           # NOT direct light only
        sim.sphere_diameter = -1 # classic radiosity
    elif cas == 4:
        # Case 4: "projected mean fluxes" (SAIL) on  toric scene
        sim.direct = False           # NOT direct light only
        sim.sphere_diameter = 0 # direct + sail (no radiosity)
        sim.pattern = d/"data/filter.8"
        sim.nb_layers = 6
        sim.can_height = 21
    elif cas == 5:
        # Case 5: nested radiosity on   toric scene
        sim.direct = False           # NOT direct light only
        sim.sphere_diameter = 10 
        sim.pattern = d/"data/filter.8"    
        sim.nb_layers = 6
        sim.can_height = 21
    ## Test input incoherencies
    elif cas == -1:
        # Erreur -1: classic radiosity on  toric scene
        sim.direct = False           # NOT direct light only
        sim.sphere_diameter = -1 # classic radiosity
        sim.pattern = d/"data/filter.8"
    elif cas == -2:
        # Erreur -1: classic radiosity on  toric scene
        sim.direct = False           # NOT direct light only
        sim.sphere_diameter = 1 # nested radiosity
        sim.pattern = None

    else:
        print "\n>>>  caribu_test(case): case still not implemented, case in [-1, 1, 2, 3, 4, 5] - MC09"
        run=False
        return

    if run:    
        print "just before Caribu.run() in case ",cas 
        sim.run()  
        print "simulation nb =",len(sim.nrj)
        print sim.nrj.keys()
        return sim.nrj 
    

def main(my_arg):
    """

    MC09
    """
    print ">>> caribu.py :main(%s) starts..."%(my_arg)
    sim=Caribu()

    print ">>> caribu.py :main(): options analysis"
    # parse command line options
    import getopt
    try:
        opts,args = getopt.getopt(my_arg, "hc:s:o:p:XN:Z:D:", ["help","canfile=", "skyfile=", "optfiles=","pattern=","scatter","nb_layers=","can_height=","sphere_diameter="])
    except getopt.GetoptError, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options and  arguments
    print "opts=%s (len=%s)"%(opts, len(opts))
    print "args=%s (len=%s)"%(args, len(args))
    sim.opticals=args
    for opt,arg in opts:
        print "opt=%s, arg=%s"%(opt, arg)
        if opt in ("-h", "--help"):      
            print "caribu.py use...\n 'hc:s:p:1N:Z:D:', ['help','canfile=','skyfile=','pattern=','direct','nb_layers=','can_height=','sphere_diameter=' + optfiles (last args)"
            sys.exit(0)               
        elif opt in ("-c", "--canfile"): 
            sim.scene = arg         
        elif opt in ("-s", "--skyfile"): 
            sim.sky = arg                 
        elif opt in ("-p", "--pattern"): 
            sim.pattern = arg
            sim.infinity = True
        elif opt in ("-X", "--scatter"): 
            sim.direct=False
        elif opt in ("-N", "--nb_layers"): 
            sim.nb_layers = int(arg)
        elif opt in ("-Z", "--can_height"): 
            sim.can_height = float(arg)
        elif opt in ("-D", "--sphere_diameter"): 
            sim.sphere_diameter = arg
        else:
            print "<!> Erreur Option non trouvee: opt=%s, arg=%s"%(opt, arg)

    print ">>> caribu.py :main():  run caribu..."
    sim.show()
    sim.run()
    

if __name__ == "__main__":
    import sys
    print "caribu.py called ..."
    main(sys.argv[1:])
    # USE
    # %run caribu.py  -h
    # %run caribu.py  --help
    # %run caribu.py  -8 # Error !!
    # %run caribu.py  --canfile='data/filterT.can' --skyfile='data/zenith.light' 'data/par.opt'
    # %run caribu.py  --canfile='data/filterT.can' --skyfile='data/zenith.light' 'data/par.opt' 'data/nir.opt'
    # %run caribu.py  --canfile='data/filterT.can' --skyfile='data/zenith.light'  --pattern='data/filter.8'  'data/par.opt' 'data/nir.opt'
    # %run caribu.py  --canfile='data/filterT.can' --skyfile='data/zenith.light'  --pattern='data/filter.8' -X -Z 21 -N 6 -D 10 'data/par.opt' 'data/nir.opt'

'''
# Original caribu.csh in C-shell....
# MC00

mv *.spec temp/
mv spectral temp/
mv cropchar temp/
mv leafarea temp/
mv *.env temp/
mv *.log temp/

#envt
setenv Sc /tmp/virtualis.can

#
if ( $#argv < 7 ) then 
 echo "Syntax error: caribu.csh Ds file.light file.can nz h file.8 file1.opt [... fileN.opt]"
else
 
# Periodise (caribu le fait il ?)
# 
  echo "periodise"  -m $argv[3] -8 $argv[6] -o $Sc
 periodise -m $argv[3] -8 $argv[6] -o $Sc
# S2V
 echo "s2v" $argv[3-]
 #s2v $argv[3-] >& s2v.log
 s2v $Sc $argv[4-] >& s2v.log
 if( $status != 0) then 
   echo "S2V a plante!"
   exit(-1)
 endif # Loop sur le p.o.
 set i=1
    foreach po ($argv[7-])
   echo "==> "$po, $i
   ## MCSAIL
   cp $po".spec" spectral
   echo "  mcsail" $argv[2]
   mcsail $argv[2]>& sail-$po.log
   if( $status != 0) then 
    echo "Sail a plante!"
    exit(-1)
   endif
   mv mlsail.env $po.env
   ## Canestra
   date
   if ( $i == 1 ) then
      ## 1ere fois: calcul des FF et Bfar
      echo "Appel no. "$i ": canestrad" -M $Sc -8 $argv[6] -l $argv[2] -p $po.opt -e $po.env -s -r  $argv[1] -1   
       canestrad -M $Sc -8 $argv[6] -l $argv[2] -p $po.opt -e $po.env -s -r $argv[1]  -t /tmp/ -f caribu -v 2 >& nr-$po.log  
  else
      ## ie fois: pas de calcul des FF et Bfar
      echo "Appel no. "$i ": canestrad" -M $Sc -8 $argv[6] -l $argv[2] -p $po.opt -e $po.env -s -r  $argv[1] -1   
       canestrad -M $Sc -8 $argv[6] -l $argv[2] -p $po.opt -e $po.env -s -r $argv[1]    -t /tmp/ -w caribu  -v 2>& nr-$po.log 
  # nr-$po'_'$i.log
  endif
  date
  cp E0.dat E_$po.dat
  cp B.dat B_$po.dat
  echo " " 
 @ i = $i  + 1
 end
 # \rm leafarea spectral cropchar
endif
## FIN


### Trucs & Astuces
#
# $status = variable de retour 
#
# Exemple de calcul en variable Cshell
#setenv c `cat $1|wc -w`
#@ c = $c  / $r

'''

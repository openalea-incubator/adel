import os

class CaribuScene(object):
    """  Handles CaribuScene """
    # Contient le texte des fichiers necessaires a canestra ou produit par lui (scene, 8, opt, ligtht, FF)
    #Contient les methode pou les brique caribu, l'extraction de donnees et la conversion vers plantGL pour visualisation

    
    def __init__(self):
        self.hasScene = False
        self.scene = ""
        self.hasPattern = False
        self.pattern = "NoPattern"
        self.PO = "#PyCaribu : PO par defaut (~PAR, materiau vert\n#format e : tige,  feuille sup,  feuille inf\n# nbre d'especes\nn 1\n#1 Sol\ns d 0.15\n# espece 1\ne d 0.10   d 0.10 0.05  d 0.10 0.05\n"
        self.wavelength = "defaultPO"
        self.hasSources = False
        self.sources = "NoLightSources"
        self.hasFF = False
        self.FF = "NoFF"

    def setCan(self,canstring):
        """  Set canopy from can file string """
        self.scene = canstring
        self.hasScene = True

    def setPattern(self,pattern_string):
        """  Set pattern """
        self.pattern = pattern_string
        self.hasPattern = True

    def setOptical(self,optstring,wavelength):
        """  Set optical properties """
        self.PO = optstring
        self.wavelength = wavelength

    def setSources(self,sources_string):
        """  Set Light Sources """
        self.sources = sources_string
        self.hasSources = True
       
    def setFF(self,FF_string):
        """  Set Form factor matrix """
        self.FF = FF_string
        self.hasFF = True

    def setFF(self,FFfile):
        """  Set form factors from ff file """
        if os.path.isfile(FFfile):
            fin = open(FFfile)
            self.FF = fin.read()
            self.hasFF = True
            fin.close()


    def writeCan(self,canfile):
        """  write a canfile of the scene """ 
        if not self.hasScene:
            print "!!!Warning!!! CaribuScene has no Scene !"
        fout = open(canfile,"w")
        fout.write(self.scene)
        fout.close()

    def writePattern(self,patternfile):
        """  write a pattern file of  the scene """ 
        if not self.hasPattern:
            print "!!!Warning!!! CaribuScene has no Pattern !"
        fout = open(patternfile,"w")
        fout.write(self.pattern)
        fout.close()

    def writeOptical(self,optfile):
        """  write an optical file of  the scene """
        if self.wavelength == "defaultPO":
            print "!!!Warning!!! No PO specified, using default PO (see CaribuScene) !"
        fout = open(optfile,"w")
        fout.write(self.PO)
        fout.close()

    def getIncidentEnergy(self):
        """ return the total ammount of energy emitted by the lights of the scene """
        Ei=0
        if self.hasSources:
            sources=self.sources.splitlines()
            for s in sources :
                l=s.strip()
                if not l or l.startswith('#'):
                    continue
                col=l.split()
                #filter empty lines
                if(not col) : continue
                Ei += float(col[0])

        return(Ei)

    def writeLight(self,lightfile):
        """  write a lightfile of the sources """ 
        if not self.hasSources:
            print "!!!Warning!!! CaribuScene has no ligth sources !"
        fout = open(lightfile,"w")
        fout.write(self.sources)
        fout.close()

    def writeFF(self,FFfile):
        """  write a FFfile of the scene """ 
        if not self.hasFF:
            print "!!!Warning!!! CaribuScene has no FormFactor yet !"
        fout = open(FFfile,"w")
        fout.write(self.FF)
        fout.close()
    
    def __str__(self):
        s = """
Pattern: 
%s

Current Wavelength : %s
PO:
%s

has FF: %s
Light Sources:
%s
Scene:
%s
"""%(self.pattern,
    self.wavelength,
    self.PO,
     str(self.hasFF),
    '\n'.join(self.sources.splitlines()[0:5])+'...',
     '\n'.join(self.scene.splitlines()[0:7])+'...')
        return s

    
    
class FileCaribuScene(CaribuScene):
    """Adaptor to contruct CaribuScenes from files"""

    def __init__(self, canfile,lightfile,patternfile=None,optfile=None):
        CaribuScene.__init__(self)

        if os.path.isfile(canfile):
            fin = open(canfile)
            self.setCan(fin.read())
            fin.close()

        if os.path.isfile(lightfile):
            fin = open(lightfile)
            self.setSources(fin.read())
            fin.close()
       
        if (patternfile is not None) and os.path.isfile(patternfile):
            fin = open(patternfile)
            self.setPattern(fin.read())
            fin.close()

        if (optfile is not None) and os.path.isfile(optfile):
            waveLength=os.path.basename(optfile).split('.')[0]
            fin = open(optfile)
            self.setOptical(fin.read(),waveLength)
            fin.close()

def newFileCaribuScene(canfile,lightfile,patternfile=None,optfile=None):
    return FileCaribuScene(canfile,lightfile,patternfile,optfile)


class ObjCaribuScene(CaribuScene):
    """Adaptor to construct caribuScene from objects"""

    def __init__(self, scene_obj,light_string,pattern_tuple=None,opt_string=None,waveLength=None):
        CaribuScene.__init__(self)
        if scene_obj is not None:
            try:
                self.scene = scene_obj.to_canestra()
            except AttributeError:
                print("Scene object input to ObjCaribuScene should have a to_canestra method")
                raise
            else:
                self.hasScene = True
        if light_string is not None:
            self.setSources(light_string)
        if pattern_tuple is not None:
            pat = '\n'.join([' '.join(map(str,pattern_tuple[0])),' '.join(map(str,pattern_tuple[1])),' '])
            self.setPattern(pat)
        if opt_string is not None and waveLength is not None:
            self.setOptical(opt_string,waveLength)

def newObjCaribuScene(scene_obj=None,ligth_string=None,pattern_tuple=None,opt_string=None,waveLength=None):
    return ObjCaribuScene(scene_obj,ligth_string,pattern_tuple,opt_string,waveLength)

def getIncidentEnergy(caribu_scene):
    return caribu_scene.getIncidentEnergy(),


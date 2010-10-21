import os

class SailScene(object):
    """  Constructs and Handles MultiLayered Sail Scene """
    
    def __init__(self):
        self.layers = "No Layers"
        self.layerdesc = "No layer desc" 
        self.optical = "No Optical Undefined"
        self.light = "No lights"
       

    def setScene(self,cropcharfile,leafareafile,specfile):
        """  Set scene from cropchar and leafarea files """
        if os.path.isfile(leafareafile):
            fin = open(leafareafile)
            self.layers = fin.read()
            fin.close()
        if os.path.isfile(cropcharfile):
            fin = open(cropcharfile)
            self.layerdesc = fin.read()
            fin.close()
        if os.path.isfile(specfile):
            fin = open(specfile)
            self.optical = fin.read()
            fin.close()

    def setLight(self,lightstring):
        """ set Scen ligths from string """
        self.light = lightstring

    def writeLight(self,file):
        """ write ligth file """
        fout = open(file,"w")
        fout.write(self.lights)
        fout.close()

    def writeCropChar(self,file):
        """  write cropchar file for Sail input """ 
        fout = open(file,"w")
        fout.write(self.layerdesc)
        fout.close()

    def writeLeafArea(self,file):
        """  write leafarea file for Sail input """ 
        fout = open(file,"w")
        fout.write(self.layers)
        fout.close()

    def writeSpectral(self,file):
        """  write spectral file for Sail input """ 
        fout = open(file,"w")
        fout.write(self.optical)
        fout.close()

    

    def __str__(self):
        ldesc = self.layerdesc.splitlines()
        return '\nLayer description :\n' + self.layerdesc + '\nOptical properties of layers:\n' + self.optical + '\nLayer content (LAI,LIDF):\n' + self.layers

    
    

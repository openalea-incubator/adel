from alinea.caribu.CaribuScene import newStringCaribuScene

def StrCaribuScene(canstring, lightstring, patternstring, optstring, wavelength):
    '''    Create CaribuScene from string content of Caribu files
    '''
    return newStringCaribuScene(canstring, lightstring, patternstring, optstring, wavelength)

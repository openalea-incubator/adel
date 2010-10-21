import rpy2.robjects as robj

r = robj.r

from pkg_resources import resource_filename as rfn

# Set Numeric Locale value 
r('Sys.setlocale(category="LC_NUMERIC",locale="C")')

r.source(rfn('alinea.sensitivity','effectsPB.R'))

inputmatrixpy = robj.globalEnv['input.matrix']
#renamepy = robj.globalEnv['rename']
#saNgrainpy = robj.globalEnv['sensitivity.Ngrain']
#saNlaminapy = robj.globalEnv['sensitivity.Nlamina']
effectspy = robj.globalEnv['effects.results']


def SensitivityAnalysisPB(inputfile, outdir) :
	inputmatrixpy(inputfile)
#	renamepy()
#	saNgrainpy()
#	saNlaminapy()
	effectspy(outdir)
	
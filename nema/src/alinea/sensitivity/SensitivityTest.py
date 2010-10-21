import rpy2.robjects as robj

r = robj.r

from pkg_resources import resource_filename as rfn

# Set Numeric Locale value 
r('Sys.setlocale(category="LC_NUMERIC",locale="C")')

#rcode = resource_string('alinea.sensitivity', 'plotOutSim.R')
#r(rcode)

r.source(rfn('alinea.sensitivity','plotSA2.R'))

readOut = robj.globalEnv['ReadOut']
#sapb = robj.globalEnv['PBNEMA']
#effects = robj.globalEnv['EFFECTSNgrain']
##te = robj.globalEnv['tes']
##re = robj.globalEnv['res']

co = robj.globalEnv['cod']
resu = robj.globalEnv['resul']
re = robj.globalEnv['res']


def SensitivityTest(outdir) :
	readOut(outdir)
#	r('dev.new()')
#	r('dev.new()')
##	te()
##	re()
	co()
	resu()
	re(outdir)
#	effects()
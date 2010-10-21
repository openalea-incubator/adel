import rpy2.robjects as robj

r = robj.r

from pkg_resources import resource_filename as rfn

# Set Numeric Locale value 
r('Sys.setlocale(category="LC_NUMERIC",locale="C")')

#rcode = resource_string('alinea.nema', 'plotOutSim.R')
#r(rcode)

r.source(rfn('alinea.nema','plotOutSim.R'))

readOut = robj.globalEnv['ReadOut']
plotNdyn = robj.globalEnv['NAgreenDynamics']
plotAbs = robj.globalEnv['Absorb']
plotDM = robj.globalEnv['DMdynamics']


def plotOut(outdir) :
	readOut(outdir)
	r('dev.new()')
	plotNdyn()
	r('dev.new()')
	plotAbs()
	r('dev.new()')
	plotDM()



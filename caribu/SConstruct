# -*-python-*-
# SConstruct file replacing Makefile

from openalea.sconsx import config, environ
import os


pj = os.path.join

name = 'Caribu'

SConsignFile()
options = Variables( ['options.py'], ARGUMENTS )

conf = config.ALEAConfig(name, ['install'])
conf.UpdateOptions(options)

# Cpradal january 2009 for OpenGL compiling
# aconfigurer pour Macosx - MC09
# tools = ['opengl']
# env = config.ALEASolution(options, tools)
# oldies
env = Environment(options=options)

conf.Update(env)

# Generate Help available with the cmd scons -h
Help(options.GenerateHelpText(env))

# Set build directory
prefix = env['build_prefix']
BuildDir( prefix, '.' )

env.Prepend(CPPPATH='#/src/cpp/include')



# Build Stage

bibliotek = env.SConscript( pj(prefix, "src/cpp/bibliotek/SConscript"), exports='env')
meschach = env.SConscript( pj(prefix, "src/cpp/meschach/mesch12a/src/SConscript"), exports='env')

#src/cpp/GLprojection
#src/cpp/glCanestra/src
dirs = """
src/cpp/mc-sail/src
src/cpp/Periodise/src
src/cpp/s2v/src
src/cpp/Canestra/src
"""
# src/cpp/glCanestra/src - MC09

dirs = map(lambda x:pj(prefix,x, 'SConscript'), Split(dirs))
env.SConscript(dirs, exports='env bibliotek meschach')

Default("build")

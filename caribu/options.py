import sys

if('win32' in sys.platform):
    compiler='mingw'

#My option - MCjul08
if ('darwin' in sys.platform):
    #Modified MC july 2008
    EXTRA_CXXFLAGS=' -g '
    EXTRA_LINKFLAGS=' -g '
    debug = True


# Redirect path
import os

cdir = os.path.dirname(__file__)
pdir = os.path.join(cdir, "../../caribu")
pdir = os.path.abspath(pdir)

__path__ = [pdir] + __path__[:]

from alinea.caribu.__init__ import *

#!/usr/bin/env python

import per2py

print "test deux : "+str(per2py.deux())

for i in range(5):
    print "\nLoop #"+str(i+1)+" : ******************\n"
    per2py.periodise(['-m','Test/bac1.can','-8','Test/bac1.8','-o','plop.can'])

import sys
print sys.copyright

print "----- END of SCRIPT -----"



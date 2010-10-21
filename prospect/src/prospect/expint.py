#This script is designed to compute MATLAB expint function
#Prospect_5.py is the script where scalar_expint will be vectorize.

from scipy import exp, Inf, vectorize
from scipy import integrate as integ
from numpy import ndarray, array

def integrand(t,n,x):  
        return exp(-x*t) / t**n

def scalar_expint(n,x):  
        return integ.quad(integrand, 1, Inf, args=(n, x))[0]


def expint(x):
    if not isinstance(x, ndarray):
        x = array(x) 
    # Real MATLAB expint, first input must be 1, second the vector to compute
    a = vectorize(scalar_expint) 
    return a(1,x)


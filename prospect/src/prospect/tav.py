#February 2010 Convert the matlab tav function in Python for integration into OpenAlea.
#Converted By Alexis Comar: alexis.comar@avignon.inra.fr

# Stern F. (1964), Transmission of isotropic radiation across an
# interface between two dielectrics, Appl. Opt., 3(1):111-113.
# Allen W.A. (1973), Transmission of isotropic light across a
# dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
# 63(6):664-666.
# ***********************************************************************
#Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
# Properties Model Separating Photosynthetic Pigments, Remote Sensing of
# Environment
# ***********************************************************************


from numpy import ndarray, array, sqrt,log, zeros
import math as m

def tav(tetha, ref):
    if not isinstance(ref, ndarray):
        ref = array(ref) 
    
    ref=ref*1. #Transform potential ref in float
    
    s=ref.shape[0]
    tetha=tetha*m.pi/180.
    r2=ref**2
    rp=r2+1
    rm=r2-1
    a=(ref+1)**2/2.
    k=-(r2-1)**2/4.
    ds=m.sin(tetha)

    k2=k**(2)
    rm2=rm**2

    if tetha==0:
        f=4*ref/(ref+1.)**2
    elif tetha==m.pi/2.:
        b1=zeros((1,s))
    else:
        b1=sqrt(((ds**2-rp/2)**2+k))

    b2=ds**2-rp/2
    b=b1-b2
    ts=(k2/(6*b**3)+k/b-b/2)-(k2/(6*a**3)+k/a-a/2);
    tp1=-2*r2*(b-a)/(rp**2);
    tp2=-2*r2*rp*log(b/a)/rm2;
    tp3=r2*(b**(-1)-a**(-1))/2;
    tp4=16*r2**(2)*(r2**2+1)*log((2*rp*b-rm2)/(2*rp*a-rm2))/(rp**(3)*rm2);
    tp5=16*r2**(3)*((2*rp*b-rm2)**(-1)-(2*rp*a-rm2)**(-1))/rp**3;
    tp=tp1+tp2+tp3+tp4+tp5;
    
    return((ts+tp)/(2*ds**2))

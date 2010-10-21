#define _Multicou
#include <iostream> //.h>
using namespace std ;

#include <fstream> //.h>
#include <assert.h>

//#include <stdio> //.h>
//#include <stdlib> //.h>
#include <math.h>

#ifndef __GNUG__
//#ifndef BCC32
#include "bool.h"
//#endif
#endif  

// using namespace std ; // dans le .cpp ?

#include <T_utilitaires.h>
#define REEL double

// Type declarations (ex fortran common)

struct Msailin{
  REEL *l, *ttl, *roo, *tau, *bmu, *bnu;/* tableau de n couches */
  REEL tts, tto, psi; 
  int nbang;
  Tabdyn<REEL,2> f;	/* was [n][45] */

  Msailin(){nbang=-1;tts=tto=psi=0;}
  void alloue(int n){
    n++;
    register int i;
    if(nbang>0){
      l=new REEL[n];   ttl=new REEL[n];
      roo=new REEL[n]; tau=new REEL[n];
      bmu=new REEL[n]; bnu=new REEL[n];
      f.alloue(n,nbang+1); f.maj(0);
      for(i=0;i<n;i++)
	l[i]=ttl[i]=roo[i]=tau[i]=bmu[i]=bnu[i]=0.;
    }
    else{
      fprintf(stderr,"<!> Msailin.alloue() : allocation de f impossible car nbang pas initialise\n");
      exit(-1);
    }
  }//alloue()
  ~Msailin(){
    if (nbang >0) {
      delete [] l;delete [] ttl;delete [] roo;
      delete [] tau;delete [] bmu;delete [] bnu ;
    } // sinon plante à vide
  }//destructor
} ;
struct Limit{
    REEL ed, RSdd, RSsd, RSdo, RSso, es;
};

struct Profout{// a allouer pour n couches
  REEL transdir, transdif, trans, refdif, refdir, absc;
} ;
struct Mlayout{// a allouer pour n couches
    REEL tss, too, rdd, tdd, rsd, tsd, rdo, tdo, rso;
};

// Prototypes
#ifdef _Msail
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN void msailad(Msailin &);
EXTERN void mlayer(Msailin &, Mlayout*);

#undef EXTERN
#ifdef _Mprof
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN void mprofil(Limit& ,Profout*,Mlayout*);

#undef EXTERN
#ifdef _Multicou
#define EXTERN
#else
#define EXTERN extern
#endif
//Global variables
/* Ne fonctionne plus avec gcc-3.2
 * EXTERN double pi, rd;
 * EXTERN int N; 
 * donc declaration de pi, rd et N dans multicou.cpp */
extern double pi, rd;
extern int N;

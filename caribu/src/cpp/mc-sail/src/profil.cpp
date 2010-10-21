#define _Mprof

#include <iostream> // définit des namespaces
using namespace std ;

#include "multicou.h"

/*******************************************************************************
*             Albert OLIOSO                                                    *
*             Mars 1994                                                        *
*             Calcul des matrices de reflectance et de transmittance           *
*             et des parametres radiatifs d'un couvert multicouche             *
*                                                                              *
*******************************************************************************/
struct Coeff{
  REEL  RMdd, TMdd, RMdo, RMsd, TMsd, TMds, RMso, TMss,
    TMTdd, TMTsd, TMTds, TMTss,den;
};

void mprofil(Limit& limit,Profout* Tprof,Mlayout* Tlay){
  register int ic;
  Coeff *T,*o;
  Mlayout l;
  
  printf("profil.C::mprofil() : ~debut\n");   fflush(stdout);fflush(stderr);
  T=NULL;
  T=new Coeff[N+1];
  //sources
  REEL *ETu,*ETd;
/*  Mlayout l; redéclaration AH 21/02/2001 */
  ETu=new REEL[N+1];
  ETd=new REEL[N+1];
  ETu[0]=ETu[0]=0;

  printf("profil.C::mprofil() : ~debut\n");   fflush(stdout);fflush(stderr);
  if(T==NULL){
    fprintf(stderr,"mprofil() : T=new Coeff[N+1] no success\n");
  }
  // limite inferieure = sol nu  (reflectance)
  T[0].RMdd = limit.RSdd;
  T[0].RMsd = limit.RSsd;
  T[0].RMdo = limit.RSdo;
  T[0].RMso = limit.RSso;
  // limite superieure = rayonnement incident
  T[N].TMTss = 1.;
  T[N].TMTds = 0.;
  T[N].TMTsd = 0.;
  T[N].TMTdd = 1.;
  // Transdir(N+2) = es/(es+ed)
  Tprof[N].transdir = limit.es;
  // Transdif(N+2) = ed/(es+ed)
  Tprof[N].transdif= limit.ed;
  //  Matrices de reflectance et de transmittance  de chaque couche
  // Rajout de sources ETu, ETd - MC1201
  o=&(T[1]);
  for (ic =1; ic <N+1; ic++,o++) {
    l=Tlay[ic];
    o->den = 1.-T[ic-1].RMdd*l.rdd;
    
    o->RMdd = l.rdd+l.tdd*T[ic-1].RMdd*l.tdd/o->den;
    o->RMsd = l.rsd+(l.tss*T[ic-1].RMsd+l.tsd*T[ic-1].RMdd)*l.tdd/o->den;
    o->RMdo = l.rdo+l.tdd*(T[ic-1].RMdd*l.tdo+T[ic-1].RMdo*l.too)/o->den;
    o->RMso = (l.tss*T[ic-1].RMsd+l.tsd*T[ic-1].RMdd)*l.tdo;
    o->RMso = o->RMso + (l.tsd + l.tss*T[ic-1].RMsd*l.rdd)*T[ic-1].RMdo*l.too;
    o->RMso = l.rso + l.tss*T[ic-1].RMso*l.too+o->RMso/o->den;
    o->TMss = l.tss;
    o->TMsd = (l.tss*l.rdd*T[ic-1].RMsd+l.tsd)/o->den;
    o->TMds = 0.;
    o->TMdd = l.tdd/o->den;
    // source
    
  }//for ic
  //  Matrice de transmittance a partir du toit du couvert
  for (ic =N-1; ic >0; ic--){
    T[ic].TMTss = T[ic].TMss*T[ic+1].TMTss;
    T[ic].TMTsd = T[ic].TMsd*T[ic+1].TMTss + T[ic].TMdd*T[ic+1].TMTsd;
    T[ic].TMTds = T[ic].TMss*T[ic+1].TMTds;
    T[ic].TMTdd = T[ic].TMdd*T[ic+1].TMTdd + T[ic].TMsd*T[ic+1].TMTds;
  }//for ic
  
  /*  Transmittance :
      elle se compose de deux termes : rayonnement direct non intercepte
      rayonnement diffusant vers le sol     
      reflectance diffuse au dessus de la couche
      */
  for (ic =1; ic <N+1; ic++) {
    Tprof[ic].transdir = T[ic].TMTss*limit.es;
    //  Transdir(ic) = TMTss(ic)*es/(es+ed)
    Tprof[ic].transdif = (T[ic].TMTsd*limit.es+T[ic].TMTdd*limit.ed);
    //  Transdif(ic) = (TMTsd(ic)*es+TMTdd(ic)*ed)/(es+ed)
    Tprof[ic].trans = Tprof[ic].transdir + Tprof[ic].transdif;
  }//for ic
  for (ic =0; ic <N; ic++) {
    Tprof[ic].refdif = T[ic].RMsd * Tprof[ic+1].transdir
                     + T[ic].RMdd * Tprof[ic+1].transdif;
    Tprof[ic].refdir = T[ic].RMso * Tprof[ic+1].transdir
                     + T[ic].RMdo * Tprof[ic+1].transdif;
  }//for ic
  //for (ic =0; ic <N; ic++) { <= bizarre ic-1 en 0
  for (ic =1; ic <N; ic++) {
    Tprof[ic].absc  =   Tprof[ic+1].transdir*(1.-T[ic].RMsd);
    Tprof[ic].absc +=   Tprof[ic+1].transdif*(1.-T[ic].RMdd);
    Tprof[ic].absc += - Tprof[ic].transdir*(1.-T[ic-1].RMsd);
    Tprof[ic].absc += - Tprof[ic].transdif*(1.-T[ic-1].RMdd);
  }//for ic
   printf("profil.C::mprofil() : avant  delete [] T \n");   fflush(stdout);fflush(stderr);
  delete [] T;

  /**********************************************************************/
  // Version matricielle primaire (12/01/2000)
  //octave: M=read("Mcoef.dat");x=M\E0;y=reshape(x,4,4);y=y';y(4:-1:1,:)

  printf("profil.C::mprofil() : creation et ecriture de la matrice  \n"); 
  // limite inferieure = sol nu  (reflectance) T[0]
  // limite superieure = rayonnement incident T[N]
  int i,j;
  Tabdyn<REEL,2> M;
  M.alloue(4*(N),4*(N));
  M.maj(0);
  // Limite atm
  i=0;
  for (ic =N; ic>0 ; ic--) {
    j=N-ic;
    // comme E=ME+E0 => E0=(I-M)E
    M(i,i)= M(i+1,i+1)= M(i+2,i+2)=M(i+3,i+3)=1;		
    if(ic<N){// hors couche atm
      //Es
      M(i++,j*4-4)=-l.tss;
      //E-
      M(i  ,j*4-4)=-l.tsd;
      M(i  ,j*4-3)=-l.tdd;
      M(i++,j*4+2)=-l.rdd;
    }else// couche atm
      i+=2;

    if(ic>1){// Couche non au-dessus du sol
      l=Tlay[ic-1];
      //E+
      M(i  ,j*4  )=-l.rsd;
      M(i  ,j*4+1)=-l.rdd;
      M(i++,j*4+6)=-l.tdd;
      //Eo
      M(i  ,j*4  )=-l.rso;
      M(i  ,j*4+1)=-l.rdo;
      M(i++,j*4+6)=-l.tdo;
    }
    else{//Sol Lambertien
      printf("Sol:i=%d,j=%d, rs=%2f\n",i,j,limit.RSsd);
      //E+
      M(i  ,j*4  )=-limit.RSsd;
      M(i++,j*4+1)=-limit.RSdd;
      //Eo
      M(i  ,j*4  )=-limit.RSso;
      M(i++,j*4+1)=-limit.RSdo;
    }
  }//for ic

  //Ecriture des fichers
  FILE *mat;
  mat=fopen("Mcoef.dat","w"); 
  for (i=0; i<4*(N); i++) {
    for (j=0; j<4*(N); j++) {
       fprintf(mat,"%.4f  ",M(i,j));
    }
    fprintf(mat,"\n");
  }
  fclose(mat);
  mat=fopen("Mvec.dat","w"); 
  fprintf(mat,"1.\n");
  for (j=1; j<4*(N);j++) {
    fprintf(mat,"0.\n");
  }
  fclose(mat);
  printf("profil.C::mprofil() : matrice ecrite \n"); 
  
}//mprofil()

#undef _Mprof

/*******************************************************************
*                 canopy.H - mike 95                               *
*            classe scene de la radioxity                          *
*             Precalcul pour le solver                             *
********************************************************************/

#ifndef CANOPY
#define CANOPY

#include <iostream>
using namespace std ;

#include <cstdlib> // pour exit

#include<fstream> //.h>
#include<iomanip> //.h>

using namespace std ;

#define  _shihan
extern "C" {
#include "sparse.h"
}
#include "verbose.h"
#include "diffuseur.h"
#include "voxel.h"

// Canopy : contient les caracteristiques de la scene
// Elle contiendra les resultats du lance de la simulation
class Canopy{
  double Etot,Einit;
  reel bmin[3],vmin[3];
  reel bmax[3],vmax[3];
  bool infty;
  reel delta[2]; //sert a l'infini
 public:
  //temporary public variable
  Voxel mesh;
  ListeD<Diffuseur *> Ldiff;
  ListeD<double> Ldiff0; //liste des labels des diffuseurs du .can (bon et pas bons) - MC10
  int Timg; //Resolution de l'image projplan (Avant en #define) - 0699 (default 1536)
  //member function
  unsigned int radim; // nombre de faces visibles de la scene
  // necessaire au capteur virtuel
  unsigned int nb_face; 
  unsigned int nbcell; 
  unsigned int nbprim; 
  
  Canopy() {Etot=Einit=0.0;}
  // cree la liste des diffuseurs de la scene
  long int  parse_can(char *,char *,char *,reel *,reel*,int,char *,Diffuseur **&);
  long int  read_shm(int,char *,char *,reel *,reel*,int,char *,Diffuseur **&);
  void cstruit_grille(double Renv) {mesh.construction(bmin,bmax,Renv,Ldiff);}
  void sail_pur(VEC **Cfar,double *Esource,char* envname);

#ifdef _HD
  void calc_FF_Bfar(
		    VEC **Cfar,
		    double *Esource,
		    char * envname,
		    bool bias,
		    double denv=0.3,
		    int nbsim=1) ;
#else
  void calc_FF_Bfar(SPMAT *FF,
		    VEC **Cfar,
		    double *Esource,
		    char * envname,
		    bool bias,
		    double denv=0.3,
		    int nbsim=1) ; 
  void calc_FF(SPMAT *FF);
#endif
  
  void projplan(Vecteur &,bool,double *);
  void data3d(int tx,int ty,Vecteur &visee,bool infty,long int ** &Zno) ;
  // bool converge(double seuil);
  //void xabs(char*,double *,double*,bool normme=false);
  //void xrad(char*,double *,double*,bool normme=false);
};// Canopy


//+************ pt2pkt()
typedef  double Punkt[3];
inline void pt2pkt(Point &P, Punkt p) {
  //attention au tableau et pointeur
  p[0]=P[0];
  p[1]=P[1];
  p[2]=P[2];
}//pt2pkt()


#endif

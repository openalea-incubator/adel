//          Capteur.H , me , today

#ifndef _Capteur
#define _Capteur

#include <math.h>

#include "primitive.h"
#include "T_utilitaires.h"
#include "image.h"
#include  "rayon.h"
#include "outils.h"
#include <diffuseur.h>

class Camera {
 private:
//   Image &img; pas clean l'affectation d'un * new ... a un Type &
   Image *img;
   Image *moy;
   Image *sigma;
//   Primitive &prim;
   Primitive *prim;
   Point oeil;
   double taille;
   Rayon &ray;
   double foc;//al 
//repere image
   Point E; // reperes {prim[1] ou E , (u,v,prim.normal()) }
   Vecteur u;
   Vecteur v;   
 public:
   Camera(char * ficname="std.cam",char * imgname="out.ppm")
     { init(ficname,imgname); }
   ~Camera()
     { delete img; delete moy;delete sigma; delete prim;}
   void init(char * ficname, char * imgname);
   void calc_visi(Liste <Diffuseur*>&);
   void capte_visi(Rayon &strahl,Diffuseur *inter,int nbray);  
   void capte(Rayon &,Diffuseur *,int);
   void developpe()
    { img->sauve();}
   void maj_stat();
   void calc_stat(unsigned int );
   void variance(char *nomvar="sigma.ppm")
     { sigma->sauve(nomvar);}
 private: //fonctions-membres utilitaires
  void colorie_triangle(Diffuseur * pdif,Point &a,Point &b,Point &c, Tabdyn<double, 2>& Zbuf, Tabdyn<Diffuseur *, 2>& Zprim,Image* pimg=NULL);
};// Camera
#endif



#ifndef _VOXEL
#define _VOXEL

#include "bsp.h"

/* Voxel : creation d'un arbre BSP greffe' sur une grille cubique dont la
   taille des voxels est fonction de l'environnement proche d'un diffuseur*/ 

class Voxel{
protected:
  Tabdyn<BSP *,3> SI;
  double taille_vox;
  int nb_vox[3];
  reel *O;//origine
  
public:
  ~Voxel();
  BSP* operator() (int, int, int);
  Voxel& operator = (Voxel&);
  // Fonctions d'acces au membre prive de la classe
  unsigned int nb_voxel(int i) { return SI.maxi()[i]; }
  double taille() {return taille_vox;}
  reel size() {return taille_vox;}
  reel * origine() {return O;}
  int coord(char  axe, double P) { return (int) ((P-O[axe])/taille_vox);}
  double milieu_vox(int axe,int indi) {return (indi+0.5)*taille_vox+O[axe];}
  double sommet(int axe,int indi) {return indi*taille_vox+O[axe];}
  // effectue la creation et la numezrotation de la Grille
  void construction(reel*, reel*,double, ListeD<Diffuseur*>&); 
  void visu();
//************************************************************
private:
  // Fonction privee utile pour construction()
  //effectue la division de la grille
  void division(char axe, unsigned char *, BSP *);
  //effectue l'initialisationet la construction de la grille
  void creation(reel *, reel *,int*,ListeD<Diffuseur *>&);
  // initialise les indices de SI
  void numerotation(BSP*);
};

#endif



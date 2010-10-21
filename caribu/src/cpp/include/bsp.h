#ifndef _BSP
#define _BSP

//#include "utilitaires.cc"
#include "diffuseur.h"

// BSP

class BSP{
protected:
  int nb_diffuseurs;
  reel position_max[3];
  reel position_min[3];
  int label;
public:
  Liste<Diffuseur *> Ldiff;
public:
  ~BSP();
  // Fonctions d'acces aux membres prives de la classe
  int nb_diffuseur() { return nb_diffuseurs; }
  int nom() { return label; }
  void baptise(int nom) {label=nom;}
  reel maxi(int i) { return position_max[i]; }
  reel mini(int i) { return position_min[i]; }
  void what_in();
  
  void destruction_boite(); // detruit une boite
  // initialise la 1ere boite au volume englobant de la scene
  void volume_englobant_scene(reel *, reel *,ListeD<Diffuseur*>&);
  // determine dans quelle boite se trouve un diffuseur
  int partage_diffuseur(Diffuseur*, BSP*, int);
  // ajoute un diffuseur a la liste des diffuseur de la boite
  void ajouter(Diffuseur*);
  // copie les champs de la boite mere dans les champs de la boite fille
  void copie_boite(BSP*);
  // remplit la liste des diffuseurs de la boite
  void remplir_diffuseur(BSP*, BSP*, int);
  // decoupe une boite en deux boites filles
  void decoupe_boite(int, BSP*, BSP*, double);
  bool contient(Point &P);
};

#endif

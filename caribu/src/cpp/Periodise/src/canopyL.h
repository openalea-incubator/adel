#ifndef CANOPY
#define CANOPY

#include <iostream>
using namespace std;

#include <cstdlib> // pour exit
#include <fstream>
#include <iomanip>

#include "diffuseur.h"
// Canopy : contient les caracteristiques de la scene
// Elle contiendra les resultats du lance de rayons
class Canopy{ 
public:
  Liste<Diffuseur *> liste_diff_scene;
public:
  Canopy(bool, bool, bool, int) {
  }
  // cree la liste des diffuseurs de la scene
  void parse_can(char *,char *,reel *,reel*,bool,char *);
  void gencan(char *);
};// Canopy

#endif

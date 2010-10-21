#ifndef LUMIERE
#define LUMIERE

#include "outils.h"

// SOURCE LUMINEUSE

class Source_Lumiere
{
   public:
	int nb_photons;
        virtual Param_Inter lache_photon()=0;
        Source_Lumiere(int alancer=10)
	 { nb_photons=alancer;} 
};

// SOLEILSS

class SoleilSS : public Source_Lumiere
{
   private:
	Point origine_soleil;
	Vecteur direction_soleil;
        double bornemin[3];
        double bornemax[3];
	int n[2]; 		// Variable de travail
	reel Delta[2];		// Variable de travail
	Point pt;		// Variable de travail
	Alea a;			// Variable de travail
	int xi;			// Variable de travail

   public:
	SoleilSS(int, Vecteur&, double*, double*);
	SoleilSS () {}
	void init (int&, Vecteur&, double*, double*); // meme utilisation que le 1er constructeur
	virtual Param_Inter lache_photon(); // donne un photon a lancer
   private :
    // fonctions internes
	Point perturbe(Point&); // perturbation du point de depart
	void init_grille2d(); // calcul une repartition des points de depart  suivant une grille 2D 
	void parcours_grille2d(); // parcours de la grille 2D
};


// SOLEIL

class Soleil : public Source_Lumiere
{
   private:
	Point origine_soleil;
	Vecteur direction_soleil;
        double bornemin[3];
        double bornemax[3];
	Alea a;			// Variable de travail

   public:
	Soleil(int a, Vecteur b, double*, double*);
	Soleil () {}
	void init (int&, Vecteur&, double*, double*); // meme utilisation que le 1er constructeur
	 Param_Inter lache_photon(); // donne un photon a lancer
};
#endif





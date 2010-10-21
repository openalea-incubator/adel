#include "lumiere.h"

#define EPSILON 1e-6
//$$$$$$$$$$$$$$    SOLEIL   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
Soleil :: Soleil(int a, Vecteur b, double* bmin, double* bmax)
	: Source_Lumiere(a)
{ //  cout<<"Soleil [cstructor] bmax[2] = "<<bmax[2]<<" - bmaxZ = "<<bornemax[2]<<endl;
  init(a,  b, bmin, bmax);
}

void Soleil :: init (int& n,Vecteur& dir, double* bmin, double* bmax){
  Vecteur V;
  dir.normalise();
  // cas d'un direction oblique
  for(register int i=0;i<3;i++){
    bornemin[i]=bmin[i];
    bornemax[i]=bmax[i];
  }
  V=dir*(bornemax[2]-bornemin[2]);
  if(dir[0]>0.0)
    bornemin[0]-= V[0];
  else bornemax[0]-= V[0];
  if(dir[1]>0.0)
    bornemin[1]-= V[1];
  else bornemax[1]-= V[1];
  Point P( bornemin[0],
	   bornemin[1],
	   bornemax[2]- (bornemax[2]-bornemin[2])/500.0);
  
  nb_photons=n;
  direction_soleil=dir;
  origine_soleil=P;
}//Soleil :: init ()

Param_Inter Soleil :: lache_photon()
 {  Point P;
    Vecteur V(a.tirage()*(bornemax[0]-bornemin[0]),a.tirage()*(bornemax[1]-bornemin[1]),0.0);
    reel poids=1.0;
    Param_Inter photon;
    P=origine_soleil+V;
//    cout<<"Soleil[lache_photon] alea?="<<a.tirage()<<"\n\t  origine_photon ";P.show();V.show();
    photon.change_origine(P);
    photon.change_direction(direction_soleil);
    photon.change_poids(poids);

    return photon;
 }//Soleil::lache_photon()

//$$$$$$$$$$$$$$   SOLEILSS  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

Point SoleilSS :: perturbe(Point& pt)
{
	Point res;

	res[0]=pt[0] + a.tirage()*Delta[0];
	res[1]=pt[1] + a.tirage()*Delta[1];
	res[2]=pt[2];

	return res;
}
	
SoleilSS :: SoleilSS(int a, Vecteur& b, double* bmin, double* bmax)
	: Source_Lumiere(a)
{ 
   init(a,b,bmin, bmax);
}

void SoleilSS :: init (int& n, Vecteur& dir, double* bmin, double* bmax){
  Vecteur V;
  dir.normalise();
  // cas d'un direction oblique
  for(register int i=0;i<3;i++){
    bornemin[i]=bmin[i];
    bornemax[i]=bmax[i];
  }
  if(dir[0]>0.0)
    bornemin[0]-= V[0];
  else bornemax[0]-= V[0];
  if(dir[1]>0.0)
    bornemin[1]-= V[1];
  else bornemax[1]-= V[1];
  
  nb_photons=n;
  direction_soleil=dir;
  init_grille2d();
}

void SoleilSS :: init_grille2d()
{
	double interm;
	double taille_scene[2];

	/* calcul de la taille de la scene en x et y */
	taille_scene[0]=bornemax[0]-bornemin[0];
	taille_scene[1]=bornemax[1]-bornemin[1];

	interm = sqrt(taille_scene[0] * taille_scene[1]);
	/* calcul du nombre de carreaux de la grille 2D en x et en y */
	n[0] = (int) ((taille_scene[0] * sqrt(nb_photons)) / interm);
	n[1] = (int) ((taille_scene[1] * sqrt(nb_photons)) / interm);
	/* calcul de la taille d'un carreau en x et en y */
	Delta[0] = taille_scene[0] / n[0];
	Delta[1] = taille_scene[1] / n[1];
	nb_photons = n[0]*n[1];

	/* initialisation du point utilitaire */
	pt[0]=bornemin[0]-Delta[0];
	pt[1]=bornemin[1];
//	pt[2]=bornemax[2]-0.001;//-EPSILON;
	pt[2]=bornemax[2]-(bornemax[2]-bornemin[2])/500.0;//-EPSILON;   
	xi=-1;
        cout<<"SoleilSS[init_2d] pt[2] ="<<pt[2]<<" - bmaxZ = "<<bornemax[2]<<endl;
}

void SoleilSS :: parcours_grille2d()
{
	pt[0]+=Delta[0];
	xi++;
	if (xi >= n[0])
	{
		xi=0;
		pt[0]=bornemin[0];
		pt[1]+=Delta[1];
	}
	origine_soleil=perturbe(pt);
}

Param_Inter SoleilSS :: lache_photon()
{
	reel poids=1.0;
	Param_Inter photon;

	parcours_grille2d();
	photon.change_origine(origine_soleil);
	photon.change_direction(direction_soleil);
	photon.change_poids(poids);

	return photon;
}







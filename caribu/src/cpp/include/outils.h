#ifndef OUTILS
#define OUTILS

#include <iostream>
using namespace std ;

#include "ferrlog.h"

#ifndef __GNUG__
#ifndef WIN32
#include "bool.h"
#include "verbose.h"
#endif
#endif

//#include "matlib.h"
#include <stdio.h>
#include <time.h>
#ifndef _OSF_SOURCE
/*
 *      Useful mathmatical constants:
 *
 * M_2PI        - 2*pi
 *
 */
#define M_2PI      6.2831853071795864769E0  /*Hex  2^ 2 * 1.921FB
54442D18 */
#endif /* _OSF_SOURCE */

#include "T_geometrie.h"
#define ER_MAEZ 999999999.9


//enum Boolean {falsch, richtig};
//itoa primitif
// alloc de la chaine fait par la fonction
char * itoa(int i);

//raus() : si cond vraie alors affiche msg et ciao
void raus(bool cond, const char *msg="Cause mysterieuse ...");
inline void STOP(int i)
{ Ferr <<"$#***>  Break "<<i<<"\t\t Srike any key to continue...";
  getchar();
  Ferr<< "\n";
}//STOP()

// Si la condition est remplie, on sort en faute en signalant
// le fichier et la ligne ou la faute s'est produite. 
// Le nom du fichier de log permet de retrouver le programme fautif.
void Dehors(bool cond, const char *msg, const int ligne);

// Si  condition est remplie, on sort en faute en signalant
// un manque de mémoire, le fichier et la ligne ou la faute 
// s'est produite, et on conseille de rebooter.
// Le nom du fichier de log permet de retrouver le programme fautif.
void TestMemoire(bool cond, const char *msg, const int ligne) ;

//-******** T_MAX
template <class Type>
inline   Type T_max(Type v1, Type v2)
 {return((v1>=v2) ? v1 : v2);}

//-******** T_MIN
template <class Type>
inline   Type T_min(Type v1, Type v2)
 {return((v1<v2) ? v1 : v2);}

//-******** T_imax : gaffe au depassement
template <class Type>
inline  int T_imax(Type tab, int i, int j)
 {return((tab[i]<tab[j]) ? j : i);}

//-******** T_imin : gaffe au depassement
template <class Type>
inline  int T_imin(Type tab, int i, int j)
 {return((tab[i]<tab[j]) ? i : j);}

// PARAMETRES d'INTERSECTION
class Param_Inter{
private:
  Point origine;
  Vecteur direction;
  double poids;
public:
  int ordre;
  inline Param_Inter(Point, Vecteur, double);
  inline Param_Inter() {}
  inline Param_Inter&   operator = (Param_Inter&);
  inline Point   origin();
  inline Vecteur direct();
  inline double  poid();
  inline void    change_origine( Point&);
  inline void    change_direction( Vecteur&);
  inline void    change_poids(const double&);
  void show(){
    cout <<"paraminter : Origine";origine.show();cout<<"\ndirection ";direction.show();cout<<"\npoids = "<<poids<<" - ordre = "<<ordre<<endl;}
};//Class Param_Intetr 
 
// Alea : permet de tirer des nombres aleatoires
/*class Alea
{   public:
	Alea() { srand48(0); }
	double  tirage() { return (drand48()); }
};
*/



/***********************************************************
***********     PARAMETRES d'INTERSECTION    ***************
************************************************************/

inline Param_Inter::Param_Inter(Point a, Vecteur b, double c=1)
  : origine(a), direction(b), poids(c)
{}

inline Param_Inter& Param_Inter::operator = (Param_Inter& a){
  origine=a.origine;
  direction=a.direction;
  poids=a.poids;
  ordre=a.ordre;
  return *this;
}

inline Point Param_Inter::origin(){
  return origine;
}

Vecteur Param_Inter :: direct(){
  return direction;
}

inline double Param_Inter::poid(){
  return poids;
}

inline void Param_Inter::change_origine( Point& ori){
  origine=ori;
}

inline void Param_Inter::change_direction( Vecteur& dir){
  //	direction=dir;
  direction[0]=dir[0];
  direction[1]=dir[1];
  direction[2]=dir[2];       
}//change_direction()

inline void Param_Inter::change_poids(const double& p){
  poids=p;
}


#endif



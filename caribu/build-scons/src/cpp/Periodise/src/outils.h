#ifndef OUTILS
#define OUTILS
//#include "matlib.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>
using namespace std;
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
#define OUT 999999999.9

#define srand48 srand
#define drand48() ((double)rand()/(double)RAND_MAX)

//enum Boolean {falsch, richtig};
//itoa primitif
// alloc de la chaine fait par la fonction
char * itoa(int i);

//raus() : si cond vraie alors affiche msg et ciao
void raus(bool cond, const char *msg="Cause mysterieuse ...");
inline void STOP(int i)
{ cerr <<"$#***>  Break "<<i<<"\t\t Srike any key to continue...";
  getchar();
  cerr<<endl;
}//STOP()
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
  bool prems;
  inline Param_Inter(Point, Vecteur, double);
  inline Param_Inter() {prems=true;}
  inline Param_Inter&   operator = (Param_Inter&);
  inline Point  & origin();
  inline Vecteur &direct();
  inline double  poid();
  inline void    change_origine( Point&);
  inline void    change_direction( Vecteur&);
  inline void    change_poids(const double&);
  void show(){
    cout <<"paraminter : Origine";origine.show();cout<<"\ndirection ";direction.show();cout<<"\npoids = "<<poids<<" - ordre = "<<ordre<<endl;}
};//Class Param_Intetr 
 
// Alea : permet de tirer des nombres aleatoires
class Alea{
private:
  double unif(){return (drand48());}
public:
  Alea() { srand48(0); }
  double  tirage() { return unif(); }
  double  operator()() { return unif(); }
  ~Alea(){ }; 
};

/***********************************************************
***********     PARAMETRES d'INTERSECTION    ***************
************************************************************/

inline Param_Inter::Param_Inter(Point a, Vecteur b, double c=1)
  : origine(a), direction(b), poids(c)
{prems=true;}

inline Param_Inter& Param_Inter::operator = (Param_Inter& a){
  origine=a.origine;
  direction=a.direction;
  poids=a.poids;
  ordre=a.ordre;
  prems=true;
  return *this;
}

inline Point &Param_Inter::origin(){
  return origine;
}

Vecteur &Param_Inter :: direct(){
  return direction;
}

inline double Param_Inter::poid(){
  return poids;
}

inline void Param_Inter::change_origine( Point& ori){
  origine=ori;
  prems=true;
}

inline void Param_Inter::change_direction( Vecteur& dir){
  //	direction=dir;
  direction[0]=dir[0];
  direction[1]=dir[1];
  direction[2]=dir[2];
   prems=true;      
}//change_direction()

inline void Param_Inter::change_poids(const double& p){
  poids=p;
}


#endif



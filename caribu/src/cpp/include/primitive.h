#ifndef PRIMITIVE
#define PRIMITIVE

#define LONG_LIGNE_CAN 1024

#include <iostream>
using namespace std;

#include <cassert>

#include "T_geometrie.h"
#include "T_utilitaires.h"
#include "outils.h"

// PRIMITIVE

class Primitive{
protected:
  Point pt_intersection;
  double nom;
public:
  double name(){return nom;}
  void called(double no){nom=no;}
  void qui(){cout <<" Primitive : "<<nom<<endl;}
  Primitive(Point);
  Primitive() {}
  virtual ~Primitive() {}
  virtual Point& operator [] (int) = 0;
  virtual int nb_sommet() = 0;
  virtual Point &sommets(int) = 0;
  virtual Vecteur normal() = 0;
  virtual reel cst_equ() = 0;
  virtual void cst_equ(reel) = 0;
  virtual void calc_cst_eq()=0;
  virtual int tout_point_inf(reel&, const int&) = 0;
  virtual int tout_point_sup(reel&, const int&) = 0;
  virtual int nb_in(reel&,reel&, const int&) = 0;
  virtual void show(const char* msg="",ostream & out=cout)=0;
  virtual double intersect(Param_Inter& parag,Point *I) = 0;
  virtual double surface ()=0;
  virtual Point centre()=0;
  virtual double distance2_point(Point &)=0;// carre de la distance du point P a la primitive
  virtual bool appart_sphere(Point & O, double R2)=0;
  virtual Vecteur azi()=0; //azimuth zero
};

// POLYGONE

class Polygone : public Primitive{
protected:
  int nb_sommets;
  Point *sommet,isob;
  Vecteur normale;
  reel cst_equ_plan; /*constante de l'equation du plan du polygone : P*N-D=0 */
  float  dGS;
public:
  Polygone(){}
  Polygone( Polygone&);
  Polygone(Polygone*);
  Polygone(Liste<Point>&Lpt,double libel=0) {init(Lpt,libel);}
  void init(Liste<Point>& ,double);
  Polygone (char*, double,
	    reel* mini=NULL,reel* maxi=NULL,
	    bool valid=true);
  void init(char*, double,
	    reel* mini=NULL,reel* maxi=NULL,
	    bool valid=true);
  Polygone (string, double,
	    reel* mini=NULL,reel* maxi=NULL,
	    bool valid=true); // HA 2003

  Polygone (float(*)[3],double,reel* mini=NULL,reel* maxi=NULL);
 ~Polygone();
  void free();
  //	Polygone& operator = (Polygone&);
  Point& operator [] (int);
  // Fonctions d'acces aux membres prives de la classe
  int nb_sommet()
    { return nb_sommets; }
  Point& sommets(int i)
    { return sommet[i]; }
  Vecteur normal()
    { return normale; }
  reel cst_equ()
    { return cst_equ_plan; }
  void cst_equ(reel cst)
    { cst_equ_plan=cst; }
  // Fonctions qui permet de savoir si un polygone est contenu dans une boite ou non
  int tout_point_inf(reel&, const int&); 
  int tout_point_sup(reel&, const int&);
  int nb_in(reel&,reel&, const int&) ;
  void calc_cst_eq() {
    calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
  }
  void calcul_normale_cst_equ(Point&, Point&, Point&);
  void show(const char* msg="",ostream & out =cout);
  // { out <<"Polygone[show] "<<msg<<endl;}
  double intersect(Param_Inter& parag,Point* I);
  double surface ();
  Point centre(){return isob;}
  Vecteur azi();
  double distance2_point(Point &); // carre de la distance du point P au polygone (convexe)
 bool appart_sphere(Point & O, double R2); 
protected:
  bool dedans(const Point &I);          
        
};//Class Polygone

// TRIANGLE

class Triangle : public Polygone
{
   public:
  Triangle(){sommet = new Point[3];nb_sommets=3;nom=-1;}
	Triangle(Liste<Point>&,double);
        Triangle(Point &, Point &,Point &,double);
  // Triangle (char*,double,double*,double*){Ferr<<"a faire!\n");}
//    double intersect(Param_Inter& parag,Point *I); 
        void show(const char* msg="",ostream &out=cout)
         { out << msg<<"-Triangle :"; qui();
           out << "          " <<sommet[0][0]<<" "<<sommet[0][1]<<" "<<sommet[0][2]<<endl;
           out << "          " <<sommet[1][0]<<" "<<sommet[1][1]<<" "<<sommet[1][2]<<endl;
           out << "          " <<sommet[2][0]<<" "<<sommet[2][1]<<" "<<sommet[2][2]<<endl;
	 }
  double surface ();
};

#endif





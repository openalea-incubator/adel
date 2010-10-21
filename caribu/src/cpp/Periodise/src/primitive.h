#ifndef PRIMITIVE
#define PRIMITIVE

//#include <strstream.h>
#include <sstream>
using namespace std;
#include <assert.h>

#include "T_geometrie.h"
#include "T_utilitaires.h"
#include "outils.h"

// Active l'accleration de l'intersection Segura_Feito(1998)
#define _SEGURA


// PRIMITIVE

class Primitive{
protected:
  double nom;
  reel psurf;
  
public:
  Point pt_intersection;
  double name(){return nom;}
  void called(double no){nom=no;}
  void qui(){printf(" Primitive : %.0f\n",nom);}
  Primitive(Point);
  Primitive( Primitive &);
  Primitive() {}
  virtual ~Primitive() {}
  virtual Point& operator [] (int) = 0;
  virtual int nb_sommet() = 0;
  virtual Point &sommets(int) = 0;
  virtual Vecteur normal() = 0;
  virtual reel cst_equ() = 0;
  virtual void cst_equ(reel) = 0;
  virtual void calc_cst_eq()=0;
  virtual int tout_point_inf(reel&, char&) = 0;
  virtual int tout_point_sup(reel&, char&) = 0;
  virtual int nb_in(reel&,reel&, char&) = 0;
  virtual void show(const char* msg="",ostream & out=cout)=0;
  virtual double intersect(Param_Inter& parag,Point *I) = 0;
  double surface (){return psurf;}
  virtual Point centre()=0;
  virtual Vecteur azi()=0; //azimuth zero
};

// POLYGONE

class Polygone : public Primitive{
protected:
  int nb_sommets;
  Point *sommet;
  Vecteur normale;
  reel cst_equ_plan; /*constante de l'equation du plan du polygone : P*N-D=0 */

public:
  bool correct;
  Polygone(){}
  Polygone(Polygone&);
  Polygone(Liste<Point>&Lpt,double libel=0) {init(Lpt,libel);}
  void init(Liste<Point>& ,double);
  Polygone (char*,double,reel* mini=NULL,reel* maxi=NULL);
  Polygone (float(*)[3],double,reel* mini=NULL,reel* maxi=NULL);
  ~Polygone();
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
  int tout_point_inf(reel&, char&); 
  int tout_point_sup(reel&, char&);
  int nb_in(reel&,reel&, char&) ;
  void calc_cst_eq() {
    calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
  }
  void calcul_normale_cst_equ(Point&, Point&, Point&);
  void show(const char* msg="",ostream & out =cout);
  // { out <<"Polygone[show] "<<msg<<endl;}
  double intersect(Param_Inter& parag,Point* I);
  void  calc_surface();
  Point centre();
  Vecteur azi();
        
protected:
  bool dedans(Point I);          
        
};

// TRIANGLE

class Triangle : public Polygone{
   public:
        Triangle(){}
	Triangle(Liste<Point>&,double);
        Triangle(Point &, Point &,Point &,double);
Triangle (char*str,double nome,reel* mini=NULL,reel* maxi=NULL);
  // Triangle (char*,double,double*,double*){cerr<<"a faire!\n");}

#ifdef _SEGURA
  double intersect(Param_Inter& parag,Point* I);
#endif
//    double intersect(Param_Inter& parag,Point *I); 
        void show(const char* msg="",ostream &out=cout)
         { out << msg<<"-Triangle :"; qui();
           out << "          " <<sommet[0][0]<<" "<<sommet[0][1]<<" "<<sommet[0][2]<<endl;
           out << "          " <<sommet[1][0]<<" "<<sommet[1][1]<<" "<<sommet[1][2]<<endl;
           out << "          " <<sommet[2][0]<<" "<<sommet[2][1]<<" "<<sommet[2][2]<<endl;
	 }
        void calc_surface ();      
};

#endif





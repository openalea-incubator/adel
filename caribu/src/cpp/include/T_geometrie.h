#ifndef GEOMETRIE
#define GEOMETRIE

// OC 14 4 98 : "constification", ie ajout des const manquants et
//		creation de deux [] distincts, cf ligne 26 Vs 27
//    29 4 98 : mise en comm de copycons et op= superflus

#include<iostream> //.h>
using namespace std ;

// flux d'erreur
#include "ferrlog.h"

#ifndef __GNUG__
//#  ifndef BCC32
#    include "bool.h"
#    include "verbose.h"
//#  endif 
#endif

extern "C" {
#include <stdlib.h>
#include <assert.h>
}
#include <math.h>



// MC98 typedef double reel;
#ifdef _SGI
typedef REELLE reel;
#else
//Pour linux
typedef float  reel; 
#endif
// Point ou vecteur defini de facon homogene
class Vecteur;
class Homogene{
protected:
  reel homo[4];

public:
  Homogene (const reel a=0, const reel b=0, const reel c=0, const reel d=0);
//OC  Homogene (const Homogene&);  // constructeur par copie // vire par OC le 29 4 98
  const reel& operator [] (const int i) const { return homo[i]; }
	reel& operator [] (const int i) {return homo[i]; }
  bool operator== (const Homogene X) const 
  { if(homo[0]!=X.homo[0]) return false;
  if(homo[1]!=X.homo[1]) return false;
  if(homo[2]!=X.homo[2]) return false;
  if(homo[3]!=X.homo[3]) return false;
  return true;
  }
//OC  Homogene& operator = ( const Homogene&); // vire par OC le 29 4 98
  Homogene operator * (const double&)const ;  // multiplication par un scalaire
  Homogene& operator *= (const double&);
  Homogene operator / (const double&)const ;  // division par un scalaire
  Homogene& operator /= (const double&); 
  //  Homogene operator * ( class Matrice4&); //multiplication par une matrice
  Homogene chgt_base(const Vecteur &u, const Vecteur &v, const Vecteur &w) const;
  void show() const {
    Ferr <<" Homogene : "<<homo[0]<<"  "<<homo[1]<<"  "<<homo[2]<<"  "<<homo[3]<<'\n';}
};

// POINT

class Point : public Homogene{
public:
  Point (const reel a=0, const reel b=0, const reel c=0);
  Point (const  Vecteur &V);
  Point (const Homogene &);          
  Point operator + (const Homogene&)const; // addition d'un point et d'un vecteur
  Point operator - (const Homogene&)const; // soustraction d'un point et d'un vecteur
  Point& operator += (const Vecteur&);  // Translation de V
  Point& operator -= (const Vecteur&);  // Translation de -V   
  reel distance_pt_plan (const Point&, const Point&, const Point&) const;
  double dist(const Point& A) const
  { return sqrt(((homo[0]-A[0])*(homo[0]-A[0]))+((homo[1]-A[1])*(homo[1]-A[1]))+((homo[2]-A[2])*(homo[2]-A[2])));}
  double dist2(const Point& A) const
  { return ((homo[0]-A[0])*(homo[0]-A[0]))+((homo[1]-A[1])*(homo[1]-A[1]))+((homo[2]-A[2])*(homo[2]-A[2]));}
};
// VECTEUR

class Vecteur : public Homogene{
public:
  Vecteur (const reel a=0, const reel b=0, const reel c=0);
  Vecteur (const Point &P, const Point &P2);
  Vecteur (const Point &P);
  Vecteur(const Homogene &);
  Vecteur operator + (const Vecteur&)const ; // addition de deux vecteurs
  Vecteur& operator += (const Vecteur&);
  Vecteur operator - (const Vecteur&)const ; // soustraction de deux vecteurs
  Vecteur& operator -= (const Vecteur&);
  bool colineaire(const Vecteur& V)const 
  { return (fabs(prod_scalaire(V))<1.0)? false:true;  }
  Vecteur operator - () const; // operateur unaire
  double prod_scalaire(const Vecteur&)const ; // produit scalaire
  Vecteur prod_vectoriel(const Vecteur&)const ; // produit vectoriel
  double prod_mixte(const Vecteur&, const Vecteur&)const ; // produit mixte
  double norme()const ; // renvoie la norme d'un vecteur
  Vecteur normalisation()const ; // normalise un vecteur
  void normalise(); // normalise un vecteur
  void formation_vecteur(const Point&, const Point&);
  Vecteur rotate(const Vecteur &axe, const double teta)const ;
};



/*************************************************
          Definitions des fonctions-membres
**************************************************/
// la defintition des fonctions-membres non Template
// a etre deportee dans _fgeometrie.C 

#define _HOMLINE
#ifdef _HOMLINE
//Homogene

inline Homogene::Homogene (const reel a, const reel b, const reel c, const reel d){
  homo[0]=a;
  homo[1]=b;
  homo[2]=c;
  homo[3]=d;
}

//OCinline Homogene::Homogene (const Homogene &A){
//OC  homo[0]=A[0];
//OC  homo[1]=A[1];
//OC  homo[2]=A[2];
//OC  homo[3]=A[3];
//OC}

/*inline reel& Homogene::operator[] (const int i){ 
  //assert(i >= 0 && i < 4);
	return homo[i]; 
}
*/
//OC inline Homogene&  Homogene::operator = (const Homogene& A){
//OC   homo[0]=A[0];
//OC   homo[1]=A[1];
//OC   homo[2]=A[2];
//OC   homo[3]=A[3];
//OC 
//OC   return *this;
//OC }

inline Homogene Homogene::operator * (const double& a)const {
  Homogene res;

  res[0]=homo[0]*a;
  res[1]=homo[1]*a;
  res[2]=homo[2]*a;
  res[3]=homo[3];
  return res;
}

inline Homogene& Homogene::operator *= (const double& a){
  homo[0]*=a;
  homo[1]*=a;
  homo[2]*=a;
  return *this;
}

inline Homogene Homogene::operator / (const double& a)const {
  Homogene res;

  assert(a!=0);
  res[0]=homo[0]/a;
  res[1]=homo[1]/a;
  res[2]=homo[2]/a;
  res[3]=homo[3];
  return res;
}

inline Homogene&  Homogene::operator /= (const double& a){
  assert(a!=0);
  homo[0]/=a;
  homo[1]/=a;
  homo[2]/=a;

  return *this;
}
template <class Type>
int min3 (Type a, Type b, Type c)
{
        if (a<=b && a<=c) return 0;
        else
                if (b<=c) return 1;
                else
			return 2;
}

#ifdef _FUC
inline Homogene Homogene::operator * (Matrice4& M){
  Homogene res;
  register int i, j;

  for (i=0; i<3; i++)
    for (j=0; j<4; j++)
      res[i]+=M[i][j]*homo[j];
  res[3]=homo[3];
  return res;
}
#endif
inline Homogene Homogene::chgt_base(const Vecteur &u, const Vecteur &v, const Vecteur &w) const{
  register int i;
  Homogene tmp;

  tmp[0]=u[0]*homo[0]+u[1]*homo[1]+u[2]*homo[2];
  tmp[1]=v[0]*homo[0]+v[1]*homo[1]+v[2]*homo[2];
  tmp[2]=w[0]*homo[0]+w[1]*homo[1]+w[2]*homo[2];
  tmp[3]=homo[3];
  return tmp;  
}// Homogene.chgt_base()

#endif //_HOMLINE

#ifdef _FUC
// VECTEUR UNITE

class Vecteur_Unite : public Vecteur
{
   public:
	Vecteur_Unite(const reel a=0, const reel b=0, const reel c=0);
	Vecteur_Unite(const Vecteur&); // constructeur par copie
};


// BASE ORTHONORMALE

class Base_orthonormale
{
   private:
	Vecteur_Unite u, v, w;

   public:
	Base_orthonormale(Vecteur_Unite a, Vecteur_Unite b, Vecteur_Unite c)
	: u(a), v(b), w(c)
	{}
	Base_orthonormale() {}
};


// MATRICE 4*4

class Matrice4
{
   private:
	Homogene M[4];

   public:
	Matrice4 (Homogene a, Homogene b, Homogene c, Homogene d);
	Matrice4 () {}
	Homogene& operator [] (int i);
	Matrice4& operator = (Matrice4& A);
};


// MATRICE de taille indefini

template <class Type>
class Matrice
{
   public:
	int n, m;
	Type *M;

   public:
	Matrice(int a=4, int b=4) : n(a), m(b)
	{
		M = new Type[a*b];
	}
	~Matrice();
	Type& operator () (int, int);
	Matrice& operator = (Matrice&);
	Matrice operator * (Matrice&); // multiplication de deux matrices
};

template <class Type>
int min3 (Type a, Type b, Type c)
{
        if (a<=b && a<=c) return 0;
        else
                if (b<=c) return 1;
                else
			return 2;
}

// MATRICE de taille indefini

template <class Type>
Matrice<Type> :: ~Matrice()
{
	delete M;
}

template <class Type>
Type& Matrice<Type> :: operator () (int i, int j)
{
	if (i<n && j<m)
		return (M[m*i + j]);
	else
	{
		Ferr << "ERREUR d'indice dans la matrice"<<'\n';
		exit (31);
	}
}

template <class Type>
Matrice<Type>& Matrice<Type> :: operator =
	(Matrice<Type>& A)
{
	int i, j;

	if(n != A.n || m != A.m)
	{
		Ferr << "ERREUR - Impossible d'egaliser 2 matrices de taille differente"<<'\n';
		cout.flush();
		exit (32);
	}
	else
	{
		for (i=0; i<n; i++)
			for (j=0; j<m; j++)
				M[m*i+j]=A(i,j);
	}

	return *this;
}

template <class Type>
Matrice<Type>  Matrice<Type> :: operator * (Matrice<Type>& A)
{
	int i, j, k;
	Matrice<Type> res(n, A.m);


	if (m != A.n)
	{
		Ferr<<"ERREUR-multiplication de matrices imcompatibles"<<'\n';
		exit(33);
	}
	else
	{
		for (i=0; i<n; i++)
			for (j=0; j<A.m; j++)
			{
				Type somme=0;
				for (k=0; k<m; k++)
					somme+=M[m*i+k]*A(k,j);
				res(i,j)=somme;
			}
	}
	
	return res; 
}
#endif //_FUC
#endif




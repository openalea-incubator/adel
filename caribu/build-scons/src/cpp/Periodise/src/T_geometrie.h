#ifndef GEOMETRIE
#define GEOMETRIE

extern "C" {
#include <stdlib.h>
#include <assert.h>
}
#include <cmath>
#include<iostream>
using namespace std;
#include "bool.h"
typedef float reel;

// Point ou vecteur defini de facon homogene
class Vecteur;
class Homogene{
  /*protected:
  reel homo[4];
  */
public:
 reel homo[4];
  Homogene (reel a=0, reel b=0, reel c=0, reel d=0);
  // Homogene (Homogene&);  // constructeur par copie
  Homogene (const Homogene&);  // constructeur par copie
  reel& operator [] (int i) {return homo[i]; }
  bool operator== (Homogene X)
  { if(homo[0]!=X.homo[0]) return false;
  if(homo[1]!=X.homo[1]) return false;
  if(homo[2]!=X.homo[2]) return false;
  if(homo[3]!=X.homo[3]) return false;
  return true;
  }
  Homogene& operator = ( Homogene A){
    homo[0]=A.homo[0];
    homo[1]=A.homo[1];
    homo[2]=A.homo[2];
    homo[3]=A.homo[3];
    return *this;
  };
  Homogene operator * ( double);  // multiplication par un scalaire
  Homogene& operator *= ( double);
  Homogene operator / ( double);  // division par un scalaire
  Homogene& operator /= (double); 
  //  Homogene operator * ( class Matrice4&); //multiplication par une matrice
  Homogene chgt_base(Vecteur &u,Vecteur &v,Vecteur &w);
  void show()
  {cout <<" Homogene : "<<homo[0]<<"  "<<homo[1]<<"  "<<homo[2]<<"  "<<homo[3]<<"  "<<endl;}
};

// POINT

class Point : public Homogene{
public:
  Point (reel a=0, reel b=0, reel c=0);
  Point (class Vecteur &V);
  // Point (const Homogene &);
  Point (const Homogene );
  
  Point operator + (const Homogene&); // addition d'un point et d'un vecteur
  Point operator - (const Homogene&); // soustraction d'un point et d'un vecteur
  Point& operator += (const Vecteur&);  // Translation de V
  Point& operator -= (const Vecteur&);  // Translation de -V   
  reel distance_pt_plan (Point&, Point&, Point&);
  double dist(Point& A)
  { return sqrt(((homo[0]-A[0])*(homo[0]-A[0]))+((homo[1]-A[1])*(homo[1]-A[1]))+((homo[2]-A[2])*(homo[2]-A[2])));}
  double dist2(Point& A)
  { return ((homo[0]-A[0])*(homo[0]-A[0]))+((homo[1]-A[1])*(homo[1]-A[1]))+((homo[2]-A[2])*(homo[2]-A[2]));}
};
// VECTEUR

class Vecteur : public Homogene{
public:
  Vecteur (reel a=0, reel b=0, reel c=0);
  Vecteur (Point &P,Point &P2);
  Vecteur (const Point P);
  Vecteur(const Homogene);
  //  Vecteur& operator = (Vecteur&u){return Homogene::operator=((Homogene&)u);}

 Vecteur operator + (const Vecteur&); // addition de deux vecteurs
  Vecteur& operator += (const Vecteur&);
  Vecteur operator - (const Vecteur&); // soustraction de deux vecteurs
  Vecteur& operator -= (const Vecteur&);
  bool colineaire(Vecteur& V)
  { return (fabs(prod_scalaire(V))<1.0)? false:true;  }
  Vecteur operator - (); // operateur unaire
  double prod_scalaire(const Vecteur&); // produit scalaire
  Vecteur prod_vectoriel(const Vecteur&)const; // produit vectoriel
  double prod_mixte(const Vecteur&, const  Vecteur&)const; // produit mixte
  double norme(); // renvoie la norme d'un vecteur
  const Vecteur normalisation(); // normalise un vecteur
  void normalise(); // normalise un vecteur
  void formation_vecteur(Point&,Point&);
  Vecteur rotate(Vecteur &axe,double teta);
};



/*************************************************
          Definitions des fonctions-membres
**************************************************/
// la defintition des fonctions-membres non Template
// a etre deportee dans _fgeometrie.C 

//#define _HOMLINE
#ifdef _HOMLINE
//Homogene

inline Homogene::Homogene (reel a, reel b, reel c, reel d){
  homo[0]=a;
  homo[1]=b;
  homo[2]=c;
  homo[3]=d;
}

inline Homogene::Homogene (const Homogene &A){
  homo[0]=A.homo[0];
  homo[1]=A.homo[1];
  homo[2]=A.homo[2];
  homo[3]=A.homo[3];
}

/*inline reel& Homogene::operator[] (int i){ 
  //assert(i >= 0 && i < 4);
	return homo[i]; 
}
*/
/*inline Homogene&  Homogene::operator = (  Homogene& A){
  homo[0]=A.homo[0];
  homo[1]=A.homo[1];
  homo[2]=A.homo[2];
  homo[3]=A.homo[3];

  return *this;
}
*/
inline Homogene Homogene::operator * (double a){
  Homogene res;

  res[0]=homo[0]*a;
  res[1]=homo[1]*a;
  res[2]=homo[2]*a;
  res[3]=homo[3];
  return res;
}

inline Homogene& Homogene::operator *= (double a){
  homo[0]*=a;
  homo[1]*=a;
  homo[2]*=a;
  return *this;
}

inline Homogene Homogene::operator / (double a){
  Homogene res;

  assert(a!=0);
  res[0]=homo[0]/a;
  res[1]=homo[1]/a;
  res[2]=homo[2]/a;
  res[3]=homo[3];
  return res;
}

inline Homogene&  Homogene::operator /= (double a){
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
inline Homogene Homogene::chgt_base(Vecteur &u,Vecteur &v,Vecteur &w){
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
	Vecteur_Unite(reel a=0, reel b=0, reel c=0);
	Vecteur_Unite(Vecteur&); // constructeur par copie
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
		cout << "ERREUR d'indice dans la matrice\n"; cout.flush();
		exit (1);
	}
}

template <class Type>
Matrice<Type>& Matrice<Type> :: operator =
	(Matrice<Type>& A)
{
	int i, j;

	if(n != A.n || m != A.m)
	{
		cout << "ERREUR - Impossible d'egaliser 2 matrices de taille differente\n";
		cout.flush();
		exit (1);
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
		cout << "ERREUR-multiplication de matrices imcompatibles\n"; cout.flush();
		exit(1);
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






































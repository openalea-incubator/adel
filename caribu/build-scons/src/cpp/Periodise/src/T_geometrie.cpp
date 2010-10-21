#include <cstdio>
using namespace std;

#include "T_geometrie.h"

#ifndef _HOMLINE
//Homogene

Homogene::Homogene (reel a, reel b, reel c, reel d){
  homo[0]=a;
  homo[1]=b;
  homo[2]=c;
  homo[3]=d;
}

/*Homogene::Homogene (Homogene &A){
  homo[0]=A.homo[0];
  homo[1]=A.homo[1];
  homo[2]=A.homo[2];
  homo[3]=A.homo[3];
  //printf("Homogene:Cstructeur par copie (mis PAS EN  const) z=%g\n",A[2]);
}
*/
Homogene::Homogene (const Homogene &A){
  homo[0]=A.homo[0];
  homo[1]=A.homo[1];
  homo[2]=A.homo[2];
  homo[3]=A.homo[3];
  //printf("Homogene:Cstructeur par copie (mis en const) z=%g\n",A[2]);
}

/*reel& Homogene::operator[] (int i){ 
  //assert(i >= 0 && i < 4);
	return homo[i]; 
}

Homogene& Homogene::operator = ( Homogene A){
  homo[0]=A.homo[0];
  homo[1]=A.homo[1];
  homo[2]=A.homo[2];
  homo[3]=A.homo[3];

  return *this;
}
*/
Homogene Homogene::operator * (double a){
  Homogene res;

  res[0]=homo[0]*a;
  res[1]=homo[1]*a;
  res[2]=homo[2]*a;
  res[3]=homo[3];
  return res;
}

Homogene& Homogene::operator *= (double a){
  homo[0]*=a;
  homo[1]*=a;
  homo[2]*=a;
  return *this;
}

Homogene Homogene::operator / (double a){
  Homogene res;

  assert(a!=0);
  res[0]=homo[0]/a;
  res[1]=homo[1]/a;
  res[2]=homo[2]/a;
  res[3]=homo[3];
  return res;
}

Homogene& Homogene::operator /= (double a){
  assert(a!=0);
  homo[0]/=a;
  homo[1]/=a;
  homo[2]/=a;

  return *this;
}
#ifdef _FUC
Homogene Homogene::operator * (Matrice4& M){
  Homogene res;
  register int i, j;

  for (i=0; i<3; i++)
    for (j=0; j<4; j++)
      res[i]+=M[i][j]*homo[j];
  res[3]=homo[3];
  return res;
}
#endif
Homogene  Homogene::chgt_base(Vecteur &u,Vecteur &v,Vecteur &w){
  Homogene tmp;

  tmp[0]=u[0]*homo[0]+u[1]*homo[1]+u[2]*homo[2];
  tmp[1]=v[0]*homo[0]+v[1]*homo[1]+v[2]*homo[2];
  tmp[2]=w[0]*homo[0]+w[1]*homo[1]+w[2]*homo[2];
  tmp[3]=homo[3];
  return tmp;  
}// Homogene.chgt_base()

#endif
// VECTEUR

Vecteur::Vecteur (reel a, reel b, reel c){
  homo[0]=a;
  homo[1]=b;
  homo[2]=c;
  homo[3]=0;
}
Vecteur::Vecteur (const Point P){
  homo[0]=P.homo[0];
  homo[1]=P.homo[1];
  homo[2]=P.homo[2];
  homo[3]=0;
}

Vecteur::Vecteur (Point &P, Point &P2){
  homo[0]=P2.homo[0]-P.homo[0];
  homo[1]=P2.homo[1]-P.homo[1];
  homo[2]=P2.homo[2]-P.homo[2];
  homo[3]=0;
}
Vecteur::Vecteur (const Homogene P){
  homo[0]=P.homo[0];
  homo[1]=P.homo[1];
  homo[2]=P.homo[2];
  homo[3]=0;
}
Vecteur Vecteur::operator + (const Vecteur &v1){
  Vecteur res;
  res[0] = v1.homo[0]+homo[0];
  res[1] = v1.homo[1]+homo[1];
  res[2] = v1.homo[2]+homo[2];
  res[3] = 0;
  return res;
}

Vecteur& Vecteur::operator += (const Vecteur &vect){
  homo[0]+=vect.homo[0];
  homo[1]+=vect.homo[1];
  homo[2]+=vect.homo[2];
  homo[3]=0;

  return *this;
}

Vecteur Vecteur::operator - (const Vecteur &v1){
  Vecteur res;

  res[0] = homo[0]-v1.homo[0];
  res[1] = homo[1]-v1.homo[1];
  res[2] = homo[2]-v1.homo[2];
  res[3] = 0;

  return res;
}

Vecteur& Vecteur::operator -= (const Vecteur &vect){
  homo[0]-=vect.homo[0];
  homo[1]-=vect.homo[1];
  homo[2]-=vect.homo[2];
  homo[3]=0;

  return *this;
}

Vecteur Vecteur::operator - (){
  Vecteur res;

  res[0]=-homo[0];
  res[1]=-homo[1];
  res[2]=-homo[2];
  res[3]=0;

  return res;
}

double Vecteur::prod_scalaire(const Vecteur &v1){
  return (v1.homo[0]*homo[0] + v1.homo[1]*homo[1] + v1.homo[2]*homo[2]);
}

Vecteur Vecteur::prod_vectoriel(const Vecteur &v1)const{
  Vecteur res;

  res[0] = v1.homo[2]*homo[1] - v1.homo[1]*homo[2];
  res[1] = v1.homo[0]*homo[2] - v1.homo[2]*homo[0];
  res[2] = v1.homo[1]*homo[0] - v1.homo[0]*homo[1];
  res[3] = 0;

  return res;
}

/*double Vecteur::prod_mixte(const Vecteur &v1,  const Vecteur &v2){
  double res;
  Vecteur u=v2;
  Vecteur inter=u.prod_vectoriel(v1);
  const Vecteur u2(homo[0],homo[1],homo[2]);
  
  res   = inter.prod_scalaire(u2);

  return res;
}
*/
double Vecteur::prod_mixte(const Vecteur &v1,  const Vecteur &v2)const{
  double res;
  Vecteur u=v2;
  Vecteur inter=u.prod_vectoriel(v1);
  
  res   = inter.prod_scalaire(*this);

  return res;
}


double Vecteur::norme(){
	return (sqrt(homo[0]*homo[0] + homo[1]*homo[1] + homo[2]*homo[2]));
}

const Vecteur Vecteur::normalisation(){
  double norm;
  Vecteur res;

  norm = norme();
  res[0] = homo[0]/norm;
  res[1] = homo[1]/norm;
  res[2] = homo[2]/norm;
  res[3] = 0;

  return res;
};
void Vecteur::normalise(){
  double norm=norme();
  if(norm<=0.0)
  { cout<<"Vecteur [normalise] norme <=0\n";
  exit(-1);
  }
  homo[0]/=norm;
  homo[1]/=norm;
  homo[2]/=norm;  
}
void Vecteur::formation_vecteur(Point& p1, Point& p2){
  //fprintf(stderr,"Vecteur::formation_vecteur : DEBUT\n");
  homo[0]=p2[0]-p1[0];
  //fprintf(stderr,"Vecteur:: BUG 1 just after\n"); fflush(stderr);
  homo[1]=p2[1]-p1[1];
  //fprintf(stderr,"Vecteur:: BUG 2 just after=> %lf, %lf\n",p2[2],p1[2] ); fflush(stderr);
  
   homo[2]=p2[2]-p1[2];
  //fprintf(stderr,"Vecteur:: BUG 3 just after\n"); fflush(stderr);
  homo[3]=0;
  //fprintf(stderr,"Vecteur::formation_vecteur : FIN\n");
}

Vecteur Vecteur::rotate( Vecteur &axe,double teta){
  Vecteur res;
  reel  cosa, cosf, sina;
  reel  ax11, ax12, ax13, ax22, ax23, ax33, ax1, ax2, ax3;

  cosa = cos(teta);
  cosf = 1.0 - cosa;
  sina = sin(teta);
  ax11 = axe[0] * axe[0] * cosf;
  ax12 = axe[0] * axe[1] * cosf;
  ax13 = axe[0] * axe[2] * cosf;
  ax22 = axe[1] * axe[1] * cosf;
  ax23 = axe[1] * axe[2] * cosf;
  ax33 = axe[2] * axe[2] * cosf;
  ax1 = axe[0] * sina;
  ax2 = axe[1] * sina;
  ax3 = axe[2] * sina;
  res[0] = (ax11+cosa)*homo[0] + (ax12-ax3)*homo[1]  + (ax13+ax2)*homo[2] ;
  res[1] = (ax12+ax3)*homo[0]  + (ax22+cosa)*homo[1] + (ax23-ax1)*homo[2] ;
  res[2] = (ax13-ax2)*homo[0]  + (ax23+ax1)*homo[1]  + (ax33+cosa)*homo[2];
  return res;
}//rotate()

// POINT

Point::Point (reel a, reel b, reel c){
  homo[0]=a;
  homo[1]=b;
  homo[2]=c;
  homo[3]=1;
}
Point::Point(Vecteur &V){
  homo[0]=V.homo[0];
  homo[1]=V.homo[1];
  homo[2]=V.homo[2];
  homo[3]=1;
}
/*Point::Point(const Homogene &V){
  homo[0]=V.homo[0];
  homo[1]=V.homo[1];
  homo[2]=V.homo[2];
  homo[3]=1;
}
*/
Point::Point(const Homogene V){
  homo[0]=V.homo[0];
  homo[1]=V.homo[1];
  homo[2]=V.homo[2];
  homo[3]=1;
}

Point Point::operator + (const Homogene &vect){
  Point res;

  if (vect.homo[3]==1){
    cout << "ERREUR - Impossible d'additionner deux points"; cout.flush();
    exit (1);
  }
  else{
    res[0]=vect.homo[0]+homo[0];
    res[1]=vect.homo[1]+homo[1];
    res[2]=vect.homo[2]+homo[2];
    res[3]=1;

    return res;
  }
}

Point Point::operator - (const Homogene &vect){
  Point res;

  if (vect.homo[3]==1)	{
    cout << "ERREUR - Impossible de soustraire deux points"; cout.flush();
    exit (1);
  }
  else	{
    res.homo[0]=homo[0]-vect.homo[0];
    res.homo[1]=homo[1]-vect.homo[1];
    res.homo[2]=homo[2]-vect.homo[2];
    res.homo[3]=1;

    return res;
  }
}

Point& Point::operator += (const Vecteur& vect){  // Translation de vect :MC 
  homo[0]+=vect.homo[0];
  homo[1]+=vect.homo[1];
  homo[2]+=vect.homo[2]; 

  return *this;
}//+=
Point& Point::operator -= (const Vecteur& vect) { // Translation de -vect   
  homo[0]-=vect.homo[0];
  homo[1]-=vect.homo[1];
  homo[2]-=vect.homo[2];

  return *this;
}//-=
reel Point::distance_pt_plan (Point& a, Point &b, Point &c){
  Vecteur AB, AC, AP;
  Vecteur ABN, ACN;
  reel res;

  /* Formation des vecteurs */
  AB.formation_vecteur(a,b);
  AC.formation_vecteur(a,c);
  AP.formation_vecteur(a,*this);

  /* Formation de la normale */
  Vecteur normale= ABN.prod_vectoriel(ACN);
  normale.normalise();

  res=AP.prod_scalaire(normale);
	
  return (res);
}


#ifdef _FUC
// VECTEUR UNITE

Vecteur_Unite :: Vecteur_Unite(reel a, reel b, reel c){
  homo[0]=a;
  homo[1]=b;
  homo[2]=c;
  homo[3]=0;

  normalisation();
}

Vecteur_Unite :: Vecteur_Unite(Vecteur &a){
  homo[0]=a[0];
  homo[1]=a[1];
  homo[2]=a[2];
  homo[3]=0;

  normalisation();
}


// MATRICE 4*4

Matrice4 :: Matrice4 (Homogene a, Homogene b, Homogene c, Homogene d){
  M[0]=a;
  M[1]=b;
  M[2]=c;
  M[3]=d;
}

Homogene& Matrice4 :: operator [] (int i){
  return (M[i]);
}

Matrice4& Matrice4 :: operator = (Matrice4& A){
  register int i, j;

  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      M[i][j]=A[i][j];

  return *this;
}

#endif

// OC 14 4 98 constification generalise du module, cf le .H
// OC 29 4 98 mise en comm de copycons et op= superflus

#include <cstdio>
using namespace std;

#include <T_geometrie.h>

#ifndef _HOMLINE
//Homogene

Homogene::Homogene (const reel a, const reel b, const reel c, const reel d){
  homo[0]=a;
  homo[1]=b;
  homo[2]=c;
  homo[3]=d;
}

//OC Homogene::Homogene (const Homogene &A){
//OC  homo[0]=A[0];
//OC  homo[1]=A[1];
//OC  homo[2]=A[2];
//OC  homo[3]=A[3];
//OC }

reel& Homogene::operator[] (const int i){ 
  //assert(i >= 0 && i < 4);
	return homo[i]; 
}

//OC Homogene& Homogene::operator = ( const Homogene& A){
//OC   homo[0]=A[0];
//OC   homo[1]=A[1];
//OC   homo[2]=A[2];
//OC   homo[3]=A[3];
//OC 
//OC   return *this;
//OC }

Homogene Homogene::operator * (const double& a)const {
  Homogene res;

  res[0]=homo[0]*a;
  res[1]=homo[1]*a;
  res[2]=homo[2]*a;
  res[3]=homo[3];
  return res;
}

Homogene& Homogene::operator *= (const double& a){
  homo[0]*=a;
  homo[1]*=a;
  homo[2]*=a;
  return *this;
}

Homogene Homogene::operator / (const double& a)const {
  Homogene res;

  assert(a!=0);
  res[0]=homo[0]/a;
  res[1]=homo[1]/a;
  res[2]=homo[2]/a;
  res[3]=homo[3];
  return res;
}

Homogene& Homogene::operator /= (const double& a){
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
Homogene  Homogene::chgt_base(const Vecteur &u, const Vecteur &v, const Vecteur &w) const{
  register int i;
  Homogene tmp;

  tmp[0]=u[0]*homo[0]+u[1]*homo[1]+u[2]*homo[2];
  tmp[1]=v[0]*homo[0]+v[1]*homo[1]+v[2]*homo[2];
  tmp[2]=w[0]*homo[0]+w[1]*homo[1]+w[2]*homo[2];
  tmp[3]=homo[3];
  return tmp;  
}// Homogene.chgt_base()

#endif
// VECTEUR

Vecteur::Vecteur (const reel a, const reel b, const reel c){
  homo[0]=a;
  homo[1]=b;
  homo[2]=c;
  homo[3]=0;
}
Vecteur::Vecteur (const Point &P){
  homo[0]=P[0];
  homo[1]=P[1];
  homo[2]=P[2];
  homo[3]=0;
}

Vecteur::Vecteur (const Point &P, const Point &P2){
  homo[0]=P2[0]-P[0];
  homo[1]=P2[1]-P[1];
  homo[2]=P2[2]-P[2];
  homo[3]=0;
}
Vecteur::Vecteur (const Homogene &P){
  homo[0]=P[0];
  homo[1]=P[1];
  homo[2]=P[2];
  homo[3]=0;
}
Vecteur Vecteur::operator + (const Vecteur &v1) const{
  Vecteur res;
  res[0] = v1[0]+homo[0];
  res[1] = v1[1]+homo[1];
  res[2] = v1[2]+homo[2];
  res[3] = 0;
  return res;
}

Vecteur& Vecteur::operator += (const Vecteur &vect){
  homo[0]+=vect[0];
  homo[1]+=vect[1];
  homo[2]+=vect[2];
  homo[3]=0;

  return *this;
}

Vecteur Vecteur::operator - (const Vecteur &v1) const{
  Vecteur res;

  res[0] = homo[0]-v1[0];
  res[1] = homo[1]-v1[1];
  res[2] = homo[2]-v1[2];
  res[3] = 0;

  return res;
}

Vecteur& Vecteur::operator -= (const Vecteur &vect){
  homo[0]-=vect[0];
  homo[1]-=vect[1];
  homo[2]-=vect[2];
  homo[3]=0;

  return *this;
}

Vecteur Vecteur::operator - () const{
  Vecteur res;

  res[0]=-homo[0];
  res[1]=-homo[1];
  res[2]=-homo[2];
  res[3]=0;

  return res;
}

double Vecteur::prod_scalaire(const Vecteur &v1) const{
  return (v1[0]*homo[0] + v1[1]*homo[1] + v1[2]*homo[2]);
}

Vecteur Vecteur::prod_vectoriel(const Vecteur &v1)const{
  Vecteur res;

  res[0] = v1[2]*homo[1] - v1[1]*homo[2];
  res[1] = v1[0]*homo[2] - v1[2]*homo[0];
  res[2] = v1[1]*homo[0] - v1[0]*homo[1];
  res[3] = 0;

  return res;
}

double Vecteur::prod_mixte(const Vecteur &v1, const Vecteur &v2) const{
  double res;
  Vecteur inter;

  inter = v2.prod_vectoriel(v1);
  res   = inter.prod_scalaire(*this);

  return res;
}

double Vecteur::norme() const{
	return (sqrt(homo[0]*homo[0] + homo[1]*homo[1] + homo[2]*homo[2]));
}

Vecteur Vecteur::normalisation()const{
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
    { Ferr<<"Vecteur [normalise] norme <=0"<<'\n';
  exit(27);
  }
  homo[0]/=norm;
  homo[1]/=norm;
  homo[2]/=norm;  
}
void Vecteur::formation_vecteur(const Point& p1, const Point& p2){
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

Vecteur Vecteur::rotate(const  Vecteur &axe,const double teta) const{
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

Point::Point (const reel a, const reel b, const reel c){
  homo[0]=a;
  homo[1]=b;
  homo[2]=c;
  homo[3]=1;
}
Point::Point(const Vecteur &V){
  homo[0]=V[0];
  homo[1]=V[1];
  homo[2]=V[2];
  homo[3]=1;
}
Point::Point(const Homogene &V){
  homo[0]=V[0];
  homo[1]=V[1];
  homo[2]=V[2];
  homo[3]=1;
}
Point Point::operator + (const Homogene &vect) const {
  Point res;

  if (vect[3]==1){
    Ferr << "ERREUR - Impossible d'additionner deux points"<<'\n';
    exit(28) ;//exit (1);
  }
  else{
    res[0]=vect[0]+homo[0];
    res[1]=vect[1]+homo[1];
    res[2]=vect[2]+homo[2];
    res[3]=1;

    return res;
  }
  
}

Point Point::operator - (const Homogene &vect) const{
  Point res;

  if (vect[3]==1)	{
    Ferr << "ERREUR - Impossible de soustraire deux points"<<'\n';
    exit (29);
  }
  else	{
    res[0]=homo[0]-vect[0];
    res[1]=homo[1]-vect[1];
    res[2]=homo[2]-vect[2];
    res[3]=1;

    return res;
  }
}

Point& Point::operator += (const Vecteur& vect){  // Translation de vect :MC 
  homo[0]+=vect[0];
  homo[1]+=vect[1];
  homo[2]+=vect[2]; 

  return *this;
}//+=
Point& Point::operator -= (const Vecteur& vect) { // Translation de -vect   
  homo[0]-=vect[0];
  homo[1]-=vect[1];
  homo[2]-=vect[2];

  return *this;
}//-=
reel Point::distance_pt_plan (const Point& a, const Point &b, const Point &c)const{
  Vecteur AB, AC, AP, normale;
  Vecteur ABN, ACN;
  reel res;

  /* Formation des vecteurs */
  AB.formation_vecteur(a,b);
  AC.formation_vecteur(a,c);
  AP.formation_vecteur(a,*this);

  /* Formation de la normale */
  normale = ABN.prod_vectoriel(ACN);
  normale = normale.normalisation();

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

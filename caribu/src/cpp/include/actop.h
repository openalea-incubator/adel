//                       actop.H

#ifndef _ACTOP
#define _ACTOP

#include <math.h>
#include <stdio.h>

//#include "outils.h"

// Actop = Acteur Optique
// Uniquement objets lambertiens

//enum Face{sup,inf};
#define Face char 
#define sup 0
#define inf 1
//cause que CC connait pas enum !
struct Transf{ 
  union {
    double teta;
    double abs;
  };
  union {
    double phi;
    double refl;
  };
  double trans;
};

class Actop{
protected:
  //  Alea rand;
public:
  Actop(){}
  virtual double rho()=0;
  virtual double tau()=0; 
  void show()
  {printf("Actop[show]\n");}    
};// Actop

// Diffus (Lambertian)
class Lambert:
 virtual public Actop{
protected:
  double reflectance;
  double transmit;
public:
  Lambert()
//    { reflectance[sup]=0.05;reflectance[inf]=0.05; transmit=0.05;} //Vis
//    { reflectance[sup]=0.49;reflectance[inf]=0.49; transmit=0.49;} //IR
//    { reflectance[sup]=0.0;reflectance[inf]=0.99; transmit=0.01;} //piege
  { reflectance=1.0; transmit=0.0;} 
  
  Lambert(double refl,double transmi)
  { reflectance=refl;  transmit=transmi;}
  double rho() {return  reflectance;}
  double tau() { return transmit;}
  void show()
  {printf("Lambert[show]\n");}  
protected:
  //  inline void lambert(double &teta, double &phi);
};// Lambert

#ifdef _NON_LAMBERT
// speculaire (Fresnel) + attenuation par "poils noirs"
class Fresnel:
  virtual public Actop{
protected:
  double Fr0;      // Fresnel en 0 [Schlick] (indice de refraction n)
  double hair;        // "hair" index
  double  fr_0(double n);
public:
  Fresnel()
  { Fr0=fr0(1.4); hair=0.1;} 
  
  Fresnel(double n, double h=0.1)
  { Fr0=fr_0(n); hair=h;} 
  void interact (Transf &param);    
  double  lux (Transf &param);//pipo (a modif qd temps : pb de directionalite -> catpeur
  void show()
  {cout<<"Fresnel[show]\n";}
protected:
inline double fresnel(double &teta, double &phi);  
};// Fresnel

// Speculaire + diffus 
class Specdifu: 
  public Lambert, 
  public Fresnel{
protected:
public:
  Specdifu(double n, double h, double refl,double transmi=0):
    Lambert(refl, transmi),Fresnel(n,h)
  {} 
  void interact (Transf &param);    
  double  lux (Transf &param); // que diffus , cf. pb Fresnel 
  void show()
  {cout<<"Specdifu[show]\n";}  
};// Specdifu

// Neural Network Bidirectional Reflectance
class NNBR : public Actop{
protected:

public:
  NNBR()
  {}
  void show()
  {cout<<"NNBR[show]\n";}  
};// NNBR
#endif

/*  ***************** Archeinformatik *******************
class Transp : public Actop{
private:
  double reflectance[2];
  double transmit[2];
  
public:
  Transp()
//    { reflectance[sup]=0.05;reflectance[inf]=0.05; transmit=0.05;} //Vis
//    { reflectance[sup]=0.49;reflectance[inf]=0.49; transmit=0.49;} //IR
//    { reflectance[sup]=0.0;reflectance[inf]=0.99; transmit=0.01;} //piege
  { reflectance[sup]=0.12;reflectance[inf]=0.12;
  transmit[sup]=0.06; transmit[inf]=0.06;} //std.opt
  
  Transp(double reflsup,double reflinf, double transup, double transinf)
  { reflectance[sup]=reflsup;reflectance[inf]=reflinf;
  transmit[sup]=transup; transmit[inf]=transinf;} 
  void interact (Transf &param);    
  double  lux (Transf &param);
  void show()
  {cout<<"Transp[show]\n";}  
  
};// Transp

class Opaque : public Actop{
private:
  double reflectance;
  double absorbtance;
  
public:
  Opaque()
  { reflectance=0.0; absorbtance=1.0;}
  Opaque(double refl)
  { reflectance=refl; absorbtance=1.0-refl;}
  void interact (Transf &param);	
  double  lux (Transf &param);   
  void show()
  {cout<<"Opaque[show]\n";}  
  
};// Opaque
******************Fin  Archeinformatik ******************
*/
#endif








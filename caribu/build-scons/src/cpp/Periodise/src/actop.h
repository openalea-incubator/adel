//                       actop.H
extern "C"{
  #include "soil.h"
}

#ifndef _ACTOP
#define _ACTOP

#include <math.h>

#include "outils.h"
#include "Mmath.h"

// Actop = Acteur Optique
// Uniquement objets lambertiens

enum Face{sup,inf};

struct Transf{
  double teta;
  double phi;
  union {
    double abs;
      double tetaeye;
  };
  union {
    double rediffuse;
    double phieye;
  };
};

class Actop{
protected:
  Alea rand;
public:
  Actop(){}
  virtual  void interact (Transf &)=0;
  virtual  double lux (Transf &)=0;
  //virtual void fini()=0;
  virtual void show()=0; // {cout<<"Actop[show]\n";}
  virtual void koi()=0;
  virtual  ~Actop(){};
};// Actop

// Diffus (Lambertian)
class Lambert:
 virtual public Actop{
protected:
  double reflectance;
  double transmit;
  //double stat[36];
public:
  Lambert()
//    { reflectance[sup]=0.05;reflectance[inf]=0.05; transmit=0.05;} //Vis
//    { reflectance[sup]=0.49;reflectance[inf]=0.49; transmit=0.49;} //IR
//    { reflectance[sup]=0.0;reflectance[inf]=0.99; transmit=0.01;} //piege
  { reflectance=0.5; transmit=0.0;} 
  
  Lambert(double refl,double transmi){
    reflectance=refl;  transmit=transmi;
    // cout<<"Lambert cstruct :ref - transm = "<<reflectance<<", "<<transmit<<endl;
    //for(int i=0;i<18;i++) stat[i]=0.0;
  } 
  void interact (Transf &param);    
  double  lux (Transf &param);
  void show()
  {cout<<"Lambert[show]\n";}
   virtual void koi() {cout<<"\n==> Fils d'Actop, I'm  Lambert\n";}
  //void fini();
    virtual ~Lambert(){};
protected:
  inline void lambert(double &teta, double &phi);
};// Lambert

// speculaire (Fresnel) + attenuation par "poils noirs"
class Fresnel:
  virtual public Actop{
protected:
  double Fr0;      // Fresnel en 0 [Schlick] (indice de refraction n)
  double hair;        // "hair" index
  double  fr_0(double n);
public:
  Fresnel()
  { Fr0=fr_0(1.4); hair=0.1;printf("Fresnel <!> Cstructeur par defaut!\n");} 
  
  Fresnel(double n, double h=0.1)
  { Fr0=fr_0(n); hair=h;printf("Fresnel cstr : Fr0 =%.4g\n",Fr0);} 
  void interact (Transf &param);    
  double  lux (Transf &param);//pipo (a modif qd temps : pb de directionalite -> catpeur
  void show()
  {cout<<"Fresnel[show]\n";}
 void koi() {cout<<"\n==> Fils d'Actop, I'm  Fresnel\n";}
    //void fini() {}
    virtual ~Fresnel(){};
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
     void koi() {cout<<"\n==> Fils d'Actop, I'm  Specdifu\n";}
  //void fini() {}
virtual ~Specdifu(){};
};// Specdifu

// lobe Gaussien:  (Fresnel) + attenuation par "poils noirs"
class Gauss:
  virtual public Fresnel{
protected:
  double zeta;
  double zeta2; //spread parameter du lobe gaussien (zeta^2)
    inline double lobeG(double &teta, double &phi);
public:
  Gauss()
  { Fr0=fr_0(1.4); hair=0.1;zeta2=0.04;printf("Gauss <!> cstructeur par defaut\n");} 
  Gauss(double n, double h, double xi): Fresnel(n,h)
  { zeta=xi; zeta2=xi*xi;
  //printf("Gauss ;-) cstructeur a param\n");
  } 
  void interact (Transf &param);    
  double  lux (Transf &param);//pipo (a modif qd temps : pb de directionalite -> capteur
  void show()
  {cout<<"Gauss[show]\n";}
  void koi() {cout<<"\n==> Fils d'Actop, I'm  Gauss\n";}
    //void fini() {}
    virtual ~Gauss(){};
};// Gauss

// Speculaire Gaussien+ diffus 
class Ross: 
  public Lambert, 
  public Gauss{
protected:
public:
  Ross(double n, double h, double xi, double refl,double transmi=0):
    Fresnel(n,h),Lambert(refl, transmi),Gauss(n,h,xi)
  {} //<!> le cstr de Fresnel n'est pas appele par Gauss !
  void interact (Transf &param);    
  double  lux (Transf &param); // que diffus , cf. pb Fresnel 
  void show()
  {cout<<"Ross[show]\n";}  
    void koi() {cout<<"\n==> Fils d'Actop, I'm  Ross\n";}  //void fini() {}
virtual ~Ross(){};
};
/****  RossL *****/
// Distribution lambertienne mais avec meme valeur de Rdh que Ross (Rspec+Rdiff)
class RossL: public Ross{ 
protected:
public:
  RossL(double n, double h, double xi, double refl,double transmi=0):
    Ross(n,h,xi,refl,transmi){} 
  void interact (Transf &param);    
  double  lux (Transf &param); // que diffus , cf. pb Fresnel 
  void show()
  {cout<<"RossL[show]\n";}  
    void koi() {cout<<"\n==> Fils d'Actop, I'm  RossL\n";} 
 //void fini() {}
virtual ~RossL(){};
};// RossL

// Hot Spot Gaussien (fresnel pour la valeur)+ diffus 
class Sol: 
  public Lambert{
protected:
    double w,h,b,c,bb,cc;// parametres de SOILSPEC
    double xi;// spread de l'approx gaussienne
    //fonction interne
    double refbd(double ts, double to, double psi);
    double refdh(double ts);
public:
    Sol(double w_, double h_, double b_, double c_, double bb_, double cc_,double xi_, double refl,double transmi=0):
      Lambert(refl, transmi){
	w=w_;
	h=h_;
	b=b_;
	c=c_;
	bb=bb_;
	cc=cc_;
	xi=xi_*xi_;
      }
    void interact (Transf &param);    
    double  lux (Transf &param); 
    void show()
  {cout<<"Sol[show]\n";}  
    //void fini() {}
    void koi() {cout<<"\n==> Fils d'Actop, I'm  Sol\n";}
    virtual ~Sol(){};
};// Sol

// Neural Network Bidirectional Reflectance
class NNBR : public Actop{
protected:
  
public:
  NNBR()
  {}
  void show()
  {cout<<"NNBR[show]\n";}
  //void fini() {}
  void koi() {cout<<"\n==> Fils d'Actop, I'm  NNBR\n";}
  virtual ~NNBR(){};
};// NNBR


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

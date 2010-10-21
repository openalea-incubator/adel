//                       actop.C

#include "actop.h"    

/************************************************************
*                      Prototypes                           *
*************************************************************/



/************************************************************
*       Fonctions utiles (locales a actop.C)                *
*************************************************************/

//-*************** lambert() **********************
#ifdef RAD
inline void  Lambert::lambert(double &teta, double &phi){
  phi=rand.tirage()*2*M_PI;
  teta=0.5*acos(1-2*rand.tirage());
}//Lambert::lambert()

#endif
/************************************************************
*                       Lambert                             *
*************************************************************/



#ifdef _NON_LAMBERT
//-*************** Fr0() ***********************
 double  Fresnel::fr_0(double n){
  double x;
    
  x=(n-1)/(n+1);
  x*=x;
  return x;
}//Fresnel::fr_0

//-*************** fresnel() ***********************
inline double Fresnel::fresnel(double &teta, double &phi){
  double Fr;
  
  phi=phi+M_PI;
  // Formule de Fresnel, par approx. Schlick (fraction rationnelle) error <1%
  Fr=1-cos(teta);
  Fr*=Fr;
  Fr*=Fr*Fr;
  Fr*=(1-Fr0);
  Fr+=Fr0;
  //Attenuation Polis Noirs (Nilson 91)
  if(hair>0.0)
    Fr*=exp(-2*hair*tan(teta)/M_PI);
  return Fr;
}//Fresnel::fresnel()

     
/************************************************************
*                       Fresnel                             *
*************************************************************/

//-*************** Fresnel::interact() ***********************
void Fresnel::interact (Transf &param){
  double Rspec=fresnel(param.teta,param.phi);
  param.abs*=1.0-Rspec;
  param.redif*=Rspec;
}//Fresnel::interact()

//-*************** Fresnel::lux() ****************************
double Fresnel::lux (Transf &param){
  // faux : fait pour pas avoir de plantage si l'image est nulle
 return  fresnel(param.teta, param.phi);

  // rajouter l'angle solide de l'objectif pour la proba d'interception du spec
/*  if ((fabs(param.tetaeye-param.teta)<1e-6)&& (fabs(fabs(param.phieye-param.phi)-M_PI)<1e-6))
    return fresnel(param.teta, param.phi);
  else
    return 0; //pas dans la dir refl. speculaire
*/
}//Fresnel::lux()
   
/************************************************************
*                       Specdifu                            *
*************************************************************/

//-*************** Specdifu::interact() ***********************
void Specdifu::interact (Transf &param){
  double Rspec=fresnel(param.teta,param.phi);
  double Rdiff=reflectance+transmit;

  if(Rspec+Rdiff>1.0){
    Ferr <<"\n\n\t#*> Violation de la loi de conservation de l'E : Rspec +Rdiff="<<Rspec+Rdiff<<endl;
    exit(1);
  }
  //cout <<"Specdifu[interact] verif conservation . E = "<<Rspec+Rdiff<<endl; 
  if(rand.tirage()<Rdiff/(Rspec+Rdiff)){
    //cas diffus
    param.phi=param.phi+M_PI; // pour raz des chgt de fresnel()
    lambert(param.teta,param.phi);
    if(transmit>0.0){
      double tir=rand.tirage();
      //cout <<" Specdifu[interact] tir = "<<tir<<" - ref = "<<reflectance<<endl;
      if(tir>=reflectance/(transmit+reflectance)){
	param.teta*=-1.0;
	// cout<<"\t\t **** TRANSMISSION ****\n";
      }//if transmis
    }//if transparent
  }//if diffus
  param.abs*=1.0-Rspec-Rdiff;
  param.redif*=Rspec+Rdiff; 
}//Specdifu::interact()

//-*************** Specdifu::lux() ****************************
double Specdifu::lux (Transf &param){
  Lambert::lux(param);
}//Specdifu::lux()


/************************************************************
*                       NNBR                               *
*************************************************************/

// a implementer

//-*************** NNBR::interact() ***********************
void NNBR::interact (Transf &param){
  NNBR(param.teta,param.phi);  }
  param.abs*=1.0-reflectance-transmit;
  param.redif*=reflectance+transmit;
}//NNBR::interact()

//-*************** NNBR::lux() ****************************
double NNBR::lux (Transf &param){
  if (param.tetaeye>=0.0)
    return reflectance/M_PI;
  else
    // cas opaque avec transmit=0 traite au niveau diffuseur (pas de pb de /0.0)
    return transmit/M_PI;
}//NNBR     ::lux()
   
#endif

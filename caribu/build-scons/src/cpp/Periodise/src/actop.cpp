//                       actop.C

#include "actop.h"

/************************************************************
*                      Prototypes                           *
*************************************************************/

/************************************************************
*                      Variables                            *
*************************************************************/
//Alea Actop::rand;
static double rad=M_PI/180.;
/************************************************************
*       Fonctions utiles (locales a actop.C)                *
*************************************************************/

/* void Lambert::fini() {
   int tetai=-90;
   FILE *fstat;
   if(stat[16]!=0.0) {
   fstat=fopen("lambert.dat","w");
   for(int i=0;i<36;i++) {
   if(i==18) tetai=5;
   fprintf(fstat,"%d    %g\n ",tetai,stat[i]);
   tetai+=5;
   }
   fclose(fstat);
   printf(" ****** ~Lambert\n");
   }
   }
   */
//-*************** lambert() **********************
inline void  Lambert::lambert(double &teta, double &phi){
  phi=rand.tirage()*M_2PI;
  //phi=0;((int)(rand.tirage()*0.99))*M_PI;
  //teta=rand.tirage()*M_PI_2;
  teta=0.5*Macos(1-2*rand.tirage());
  //teta=Masin(rand.tirage());
    /* while(fabs(teta-M_PI_2)<1e-8)
    teta=0.5*Macos(1-2*rand.tirage());
    */
  //cout<<" Lambert::lambert : teta = "<<teta<<endl;
  
}//Lambert::lambert()

//-*************** Fr0() ***********************stat[tetai+18]+=param.abs;
 double  Fresnel::fr_0(double n){
  double x;
  printf("DEBUG Fresnel::fr0() n=%.4g\n",n);    
  x=(n-1)/(n+1);
  x*=x;
  return x;
}//Fresnel::fr_0

//-*************** fresnel() ***********************
inline double Fresnel::fresnel(double &teta, double &phi){
  double Fr;
  static bool fst=true;
  phi=phi;
  
  // Formule de Fresnel, par approx. Schlick (fraction rationnelle) error <1%
  Fr=1-cos(teta);//rigueur?
  Fr*=Fr;
  Fr*=Fr*Fr;
  Fr*=(1-Fr0);
  Fr+=Fr0;
  //if(Fr>1) printf(" Fr >1 - Fr= %.4g - teta = %.4g\n",Fr,teta/rad);
  //Attenuation Poils Noirs (Nilson 91)
  if(hair>0.0) {
    Fr*=exp(-2*hair*tan(teta)/M_PI);
    //if(Fr>1) printf(" Fr_a >1 - Fr_a= %.4g - teta = %.4g\n",Fr,teta/rad);
  }
  if(fst){
    printf("=====> Fresnel(%g) = %g - Fr0=%g\n",teta/rad,Fr,Fr0);
    fst=false;
  }
  return Fr;
}//Fresnel::fresnel()

//-*************** lobeG() ***********************
inline double Gauss::lobeG(double &teta, double &phi){
  double Fr,um,g,Z, spread;
  static bool fst=true;
  double t=teta,p=phi,cost,sint,sing,sc;//necessaire a cause de la pertub angulaire de Fresnel  
  int iter=0;

  cost=cos(teta);
  sint=sin(teta);
  // Formule de Fresnel, par approx. Schlick (fraction rationnelle) error <1%
  Fr=fresnel(t,p);
  //printf(" <!> Fr=%g\n",Fr);
  // hard coder pur debuguer le cone...
  //Fr= 0.0332473; //0deg
  //Fr=0.0332529;//30deg
  //Fr*=1-cost*0.25;
  //these de Y. Govaerts (Annexe C., p. 156)
  // spread=zeta2*(0.001+cost*sqrt(cost));
  // modif 030699
  spread=zeta2*(0.001+pow(cost,(1.+2*zeta)));

  //phi_r
  phi=2.*M_PI*rand();
  //gamma
  if(fst){
    fst=false;
    printf("teta=%g, zeat2=%g ==> spread=%g\n",teta/rad, zeta2, spread);
  }

  if(0){// MC0499: marche pas terr
    if(cost==1.){
      printf("cost=1 ==> um=0\n");
      um=0;
    }
    else{
      //um=exp(-1/((1-cost*cost)/spread));
      //MODIF MIKE
      um=exp(-1/(1-cost*cost)*spread);
      //Test a la c....
      //um=0;
    }
  }
  else 
    um=0;
  do{  
    g=rand();// (1.-um)+um
    g=atan(sqrt(-spread*log(g)));
    //if(rand()>0.5) g=-g; // pour avoir les 2 brins de la gaussienne
    //teta_r
    // printf("lobeG() : g=%g,p=%g, ",g/rad,p/rad);
    sc=sin(g)*cos(phi);
    Z=cos(g)*cost - sc*sint;   
    iter++;
    if(iter==10)
      fprintf(stderr,"<10> Gauss:lobeG() boucle de calcul de g: iter= 10 => teta=%.1g, zeat2=%.6g ==> spread=%.6g\n",teta/rad, zeta2, spread);  
    if(iter==50)
      fprintf(stderr,"<50> Gauss:lobeG() boucle de calcul de g: iter= 50 => teta=%.1g, zeat2=%.6g ==> spread=%.6g\n",teta/rad, zeta2, spread);  
    if(iter==100)
      fprintf(stderr,"<100> Gauss:lobeG() boucle de calcul de g: iter= 100 => teta=%.1g, zeat2=%.6g ==> spread=%.6g\n",teta/rad, zeta2, spread);  
    if(iter==500)
      fprintf(stderr,"<500> Gauss:lobeG() boucle de calcul de g: iter= 500 => teta=%.1g, zeat2=%.6g ==> spread=%.6g\n",teta/rad, zeta2, spread);  
    
  }while(Z<=0);
  /*  whil(t>M_2PI);
      NON tester que le z est negatif sinon avec le phi ca peut passer 
	      */
      /*   printf("teta=%g, g=%g, Z=%g\n",teta/rad, g/rad,Z); */

      //Calcul du tetha
      t=Macos(Z);
      //Calcul du phi
      if(fabs(teta)<0.0001){
	//printf("phi=%g\n",phi/rad);
	p=phi;
      }
      else{
	sing=sin(g);
	p=Masin(sin(phi)*sing/sin(t));
	if((cost*sc+sint*cos(g))<0){
	  p=M_PI-p;
	  //printf("**** cas x<0 ***\n");
	}
      }
      //printf(" : teta=%g, phi=%g - t = %g\n",teta/rad,phi/rad,t/rad);
  

      teta=t;
      phi=p+M_PI;

  
      return Fr;
}//Gauss::lobeG()



/************************************************************
*                       Lambert                             *
*************************************************************/

//-*************** Lambert::interact() ***********************
void Lambert::interact (Transf &param){
  //printf(" Lambert[interact] DEBUT \n"); fflush(stdout);
  lambert(param.teta,param.phi);
    if(transmit>0.0){
      double tir=rand.tirage();
      //cout <<" Lambert[interact] tir = "<<tir<<" - ref = "<<reflectance<<endl;
      if(tir>=reflectance/(transmit+reflectance)){
	param.teta=M_PI-param.teta;
	// cout<<"\t\t **** TRANSMISSION ****\n";
      }
    }
    /* double fi;
       fi=param.phi;
       tetai=((int) (param.teta*180.0/M_PI))/5;
       fi*=180.0/M_PI;
       if((fi<=2.5) || (fi>=357.5))
       stat[17-tetai]+=param.abs;
       if(fabs(fi-180.0)<=2.5)
       stat[tetai+18]+=param.abs;
       */
  param.abs*=1.0-reflectance-transmit;
  param.rediffuse*=reflectance;//a completer dans Difusueur
}//Lambert::interact()

//-*************** Lambert::lux() ****************************
double Lambert::lux (Transf &param){
  if (param.tetaeye>=0.0)
    return reflectance/M_PI;
  else
    // cas opaque avec transmit=0 traite au niveau diffuseur (pas de pb de /0.0)
    return transmit/M_PI;
}//Lambert::lux()
   

/************************************************************
*                       Fresnel                             *
*************************************************************/

//-*************** Fresnel::interact() ***********************
void Fresnel::interact (Transf &param){
  double Rspec=fresnel(param.teta,param.phi);
  param.phi+=M_PI;
  param.abs*=1.0-Rspec;
  param.rediffuse*=Rspec;
}//Fresnel::interact()

//-*************** Fresnel::lux() ****************************
double Fresnel::lux (Transf &param){
  // faux : fait pour pas avoir de plantage si l'image est nulle
    return fresnel(param.teta, param.phi);
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
  //printf("Specdifu[interact] DEBUT rho+tau=%g\n",reflectance+transmit);
  double Rspec=fresnel(param.teta,param.phi);
  double Rdiff=(1.0-Rspec)*(reflectance+transmit);

  /*  if(Rspec+Rdiff>1.0){    
      cerr<<"\n\n\t#*> Violation de la loi de conservation de l'E : Rspec +Rdiff="<<Rspec+Rdiff<<endl;
      fprintf(stderr,"Rf= %.4g - Tf= %.4g - Rspec= %.4g \n",reflectance,transmit, Rspec);
      exit(-1);
      }
  */
  //printf("Specdifu[interact] Rspec=%g, verif conservation . E = %g\n",Rspec,Rspec+Rdiff); 
  if(rand.tirage()<Rdiff/(Rspec+Rdiff)){
    //cas diffus
    //printf("Cas diffus\n");
    lambert(param.teta,param.phi);
    if(transmit>0.0){
      double tir=rand.tirage();
      //cout <<" Specdifu[interact] tir = "<<tir<<" - ref = "<<reflectance<<endl;
      if(tir>=reflectance/(transmit+reflectance)){
	param.teta=M_PI-param.teta;
	// cout<<"\t\t **** TRANSMISSION ****\n";
      }//if transmis
    }//if transparent
  }//if diffus
  else{
    param.phi+=M_PI*(1+ 5.*(0.5-rand())/180.0);
    param.teta+=5.*(0.5-rand())*rad;
    if(param.teta>M_PI_2)
      printf(" Gargione : t=%g!\n", param.teta/M_PI*180.);
    param.teta=(param.teta<M_PI_2)?param.teta:M_PI_2;
  }
  param.abs*=1.0-Rspec-Rdiff;
  //param.rediffuse*=Rspec+Rdiff;
  param.rediffuse*=Rspec+(1.0-Rspec)*reflectance; //a boucler pae eq. bilan d'E
  //Einc=Eabs+Erefl+Etrans=Eabs+(Espec+ErefD)+Etrans
}//Specdifu::interact()

//-*************** Specdifu::lux() ****************************
double Specdifu::lux (Transf &param){
 return  Lambert::lux(param);
}//Specdifu::lux()

/************************************************************
*                       Gauss                             *
*************************************************************/

//-*************** Gauss::interact() ***********************
void Gauss::interact (Transf &param){
  double Rspec=lobeG(param.teta,param.phi);
  param.abs*=1.0-Rspec;
  param.rediffuse*=Rspec;
}//Gauss::interact()

//-*************** Gauss::lux() ****************************
double Gauss::lux (Transf &param){
  // faux : fait pour pas avoir de plantage si l'image est nulle
 return  lobeG(param.teta, param.phi);

  // rajouter l'angle solide de l'objectif pour la proba d'interception du spec
/*  if ((fabs(param.tetaeye-param.teta)<1e-6)&& (fabs(fabs(param.phieye-param.phi)-M_PI)<1e-6))
    return Gauss(param.teta, param.phi);
  else
    return 0; //pas dans la dir refl. speculaire
*/
}//Gauss::lux()
   
/************************************************************
*                       Ross                            *
*************************************************************/

//-*************** Ross::interact() ***********************
void Ross::interact (Transf &param){
  double Rspec=lobeG(param.teta,param.phi);
  double Sdiff=reflectance+transmit;
  double Rdiff=(1.0-Rspec)*Sdiff;
  
  static bool fst=true;
  if(fst){
    printf("Ross: Rspec = %.4g - Rdiff = %.4g (rd=%.4g/td=%.4g)\n",Rspec, Rdiff,(1.0-Rspec)*reflectance, (1.0-Rspec)*transmit);
    fst=false;
  }

  /*  if(Rspec+Rdiff>1.0){    
      cerr<<"\n\n\t#*> Violation de la loi de conservation de l'E : Rspec +Rdiff="<<Rspec+Rdiff<<endl;
      fprintf(stderr,"Rf= %.4g - Tf= %.4g - Rspec= %.4g \n",reflectance,transmit, Rspec);
      exit(-1);
      }
  */
  //cout <<"Ross[interact] verif conservation . E = "<<Rspec+Rdiff<<endl; 
  if(rand.tirage()<Rdiff/(Rspec+Rdiff)){
    //cas diffus
    lambert(param.teta,param.phi);
    if(transmit>0.0){
      double tir=rand.tirage();
      //cout <<" Ross[interact] tir = "<<tir<<" - ref = "<<reflectance<<endl;
      if(tir>=reflectance/Sdiff){
	param.teta=M_PI-param.teta;
	//cout<<"\t\t **** TRANSMISSION ****\n";
      }//if transmis

    }//if transparent
  }//if diffus

  param.abs*=1.0-Rspec-Rdiff;
  param.rediffuse*=Rspec+(1.0-Rspec)*reflectance; //a boucler pae eq. bilan d'E
  //Einc=Eabs+Erefl+Etrans=Eabs+(Espec+ErefD)+Etrans
  //Necessaire pour calculer  dans Difusueur le  transmis
}//Ross::interact()

//-*************** Ross::lux() ****************************
double Ross::lux (Transf &param){
 return  Lambert::lux(param);
}//Ross::lux()

/************************************************************
*                       RossL                            *
*************************************************************/

//-*************** RossL::interact() ***********************
void RossL::interact (Transf &param){
  double Rspec= fresnel(param.teta,param.phi);
  double R = Rspec+(1-Rspec)*reflectance;
  double T= (1.0-Rspec)*transmit;
  double Sdiff=reflectance+transmit;
  double Rdiff=(1.0-Rspec)*Sdiff;
  
  static bool fst=true;
  if(fst){
    printf("RossL: Rspec = %.4g - Rdiff = %.4g (rd=%.4g/td=%.4g)\n",Rspec, Rdiff,(1.0-Rspec)*reflectance, (1.0-Rspec)*transmit);
    fst=false;
  }
  //distribue angulairement la reflectance et la transmission de facon lambertienne
  lambert(param.teta,param.phi);
  if(transmit>0.0){
    double tir=rand.tirage();
    //cout <<" RossL[interact] tir = "<<tir<<" - ref = "<<reflectance<<endl;
    if(tir>=R/(R+T)){
      param.teta=M_PI-param.teta;
      //cout<<"\t\t **** TRANSMISSION ****\n";
    }//if transmis    
  }//if transparent
  
  param.abs*=1.0-Rspec-Rdiff;
  param.rediffuse*=Rspec+(1.0-Rspec)*reflectance; //a boucler pae eq. bilan d'E
  //Einc=Eabs+Erefl+Etrans=Eabs+(Espec+ErefD)+Etrans
  //Necessaire pour calculer  dans Difusueur le  transmis
}//RossL::interact()

//-*************** RossL::lux() ****************************
double RossL::lux (Transf &param){
 return  Lambert::lux(param);
}//RossL::lux()

/************************************************************
*                       Sol                            *
*************************************************************/
// Fonctions internes : calcul de SOILSPEC
/*******************************************/
/*     reflectance bidirectionnelle       */
/******************************************/
double Sol::refbd(double ts, double to, double psi){
  double g, g1, g2, ce, bg, ci, he, gg,  hi, pg, se, si;

  ci = cos(ts);
  ce = cos(to);
  si = sin(ts);
  se = sin(to);
  g = sqrt(1. - w);
  gg = g * 2.;
  g1 = ci * ce + si * se * cos(psi);
  g2 = ci * ce - si * se * cos(psi);
  pg = 1+ b*g1 + c*((3*g1*g1-1)/2.) + bb*g2 + cc*((3*g2*g2-1)/2.);
  bg = 1./(1.+tan(Macos(g1)/2.)/h);
  hi = (1+ci*2) / (1+gg*ci);
  he = (1+ce*2) / (1+gg*ce);

  return w * ((1+bg)*pg + hi*he-1) / (4.*M_PI*(ci + ce));
}
/*********************************************/
/*     reflectance directe/hemispherique     */
/*********************************************/
double Sol::refdh(double ts){
  double tmp,g, r1, r2, r3, r4, r5, r6, ci, gg,  rb, rc, dg1, dg2,Pci[3];
  
  ci = cos(ts);
  g = sqrt(1. - w);
  gg = g * 2.;
  
  Pci[1]=ci*ci;
  Pci[2]=ci*Pci[1];
  rb = b + bb;
  rc = c + cc;
  dg1 = log(gg + 1.);
  dg2 = log(gg - 1.);

  r1 = 2*(g+1)*Pci[2]*log((1+ci)/ci)/(gg*ci-1);
  r2 = (ci * 2. + 1.) * dg1 / (gg*gg * (gg* ci - 1.));
  r3 = (ci * 2. * (g + 1.) + 1.) / (gg);
  tmp=3*Pci[1]-1;
  r4 = ci * (4*rb*Pci[1] - rc*tmp*tmp - 4) * log((ci + 1.) / ci) / 4.; /* /4 ou /4. ?*/
  r5 = rb * ci * (ci * 2. - 1.) / 2.;/* /2 ou /2. ?*/
  r6 = 0.375*rc*ci* (6.*Pci[2]-3*Pci[1]-2+1);
  
  return (r1-r2-r3)*w*(g-1.)/(gg*ci+1.)+(r4-r5+r6+1.)*w/2.;
}

//Fonctions exteriorisables
//-*************** Sol::interact() ***********************
void Sol::interact (Transf &param){
  //1st sol : Comme Ross sauf que je reaxe dans le hot spot par un phi-=M_PI
  // Gauusien + lambert, mais a partir de valeur de SOILSPEC et un calage a la mano
  double ts=param.teta;
  // <!>double Rhs=refbd(ts,ts,0.);
  double Rdh=refdh(ts);
  double Rbd=refbd(ts,ts,0.)*M_PI;
  double Rdiff=reflectance;
  static bool fst=true;

  if(fst){
    printf("Sol: Rdh(%.4g) = %.4g - Rbd = %.4g - Rdh/Rbd=  %.4g\n",ts/rad,Rdh,Rbd,Rdh/Rbd);
    fst=false;
  }

  if(Rdh>1.0){    
    cerr<<"\n\n\t#*> Violation de la loi de conservation de l'E : Rdh="<<Rdh<<endl;
    exit(-1);
  }
  //cout <<"Sol[interact] verif conservation . E = "<<Rdh<<endl; 
  if(rand.tirage()<Rdiff/Rbd){
    //cas diffus
    lambert(param.teta,param.phi);
  }//if diffus
  else{ // Lobe Gaussien
    double g,t,p,cost,sint,sinp,sc;
    do{
      p=0;//M_PI*rand()-M_PI/2.;
      cost=cos(ts);
      sint=sin(ts);
      //tirage gaussien du teta /direction HS
      g=atan(sqrt(-xi*log(rand())));
      //if(rand()>0.5) g=-g;
      // chgt repere HS -> feuille
      sc=sin(g)*cos(p);
      t=Macos(cos(g)*cost - sc*sint);
      if(fabs(param.teta)>0.0001){
	p=Macos(( -cos(g)*sint - sc*cost)/sin(t));
	sinp=sin(p);
	if(fabs(sinp)<0.0001)
	  p+=M_PI;
	else
	  if(sinp>0) //car Macos -> [0:pi] et sin(g)>0
	    p=M_PI*2.-p;
      }
    }while(t>M_2PI);
    param.teta=t;
    param.phi=p-M_PI;
  }
  param.abs*=1.0-Rdh;
  param.rediffuse*=Rdh; //a boucler pae eq. bilan d'E
  //Einc=Eabs+Erefl+Etrans=Eabs+(Espec+ErefD)+Etrans
}//Sol::interact()

//-*************** Sol::lux() ****************************
double Sol::lux (Transf &param){
    if (param.tetaeye>=0.0)
    return refbd(param.teta,param.tetaeye,param.phieye);
  else
    // cas opaque avec transmit=0 traite au niveau diffuseur (pas de pb de /0.0)
    return transmit/M_PI;
}//Sol::lux()

/************************************************************
*                       NNBR                               *
*************************************************************/

// a implementer
/*
//-*************** NNBR::interact() ***********************
void NNBR::interact (Transf &param){
  NNBR(param.teta,param.phi);  }
  param.abs*=1.0-reflectance-transmit;
  param.rediffuse*=reflectance;
}//NNBR::interact()

//-*************** NNBR::lux() ****************************
double NNBR::lux (Transf &param){
  if (param.tetaeye>=0.0)
    return reflectance/M_PI;
  else
    // cas opaque avec transmit=0 traite au niveau diffuseur (pas de pb de /0.0)
    return transmit/M_PI;
}//NNBR     ::lux()
   
*/
   

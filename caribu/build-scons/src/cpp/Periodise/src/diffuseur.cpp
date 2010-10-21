#include <cmath> 
using namespace std;
#include "diffuseur.h"

//#define CO(A) A
#define CO(A)  

//-*************** face() ***********************************
inline Face face(double scal)
 { return ((scal>0.0)?inf:sup);
 }
 
//-*************** DiffO::lux() *****************************

double DiffO::lux(Param_Inter& parag){
  Transf parao;
  Vecteur Vi(parag.origin()),Ve,N;
   
   Ve=parag.direct();
   N=prim->normal();

   Ve.normalise();
   Vi.normalise();   
   N.normalise();

   //cout<<"Diff0[lux()]poid = "<<parag.poid();prim->qui();
   parao.teta=Macos(-Vi.prod_scalaire(N));
   parao.phi=0.0;
   parao.tetaeye=M_PI_2-Macos(N.prod_scalaire(Ve));
      // a revoir si mat anisotrope
   if( Ve.colineaire(N) || Vi.colineaire(N) )
      parao.phieye=0.0;
    else
     { if (Vi.colineaire(Ve))
           parao.phieye= parao.phi;
        else
	 { Vecteur Ui,Ue;
	   Ui=Vi.prod_vectoriel(N);
	   Ue=Ve.prod_vectoriel(N);
           parao.phieye= Ui.prod_scalaire(Ue);
           if(N.prod_scalaire(Ui.prod_vectoriel(Ue))>0.0)
             parao.phieye=2.0*M_PI-parao.phieye;
	 }
     }
 
   return opti->lux(parao)*parag.poid()*fabs(cos(M_PI_2-parao.tetaeye));

 }//DiffO::lux()

//-*************** DiffT::lux() *****************************
//$$$$$$$$  Vi : incident ; Ve : eye ; N : normale 
double DiffT::lux(Param_Inter& parag){
  Transf parao;
  double scal;
  Vecteur Vi(parag.origin()),Ve,N;
  Face side;

  Ve=parag.direct();
  N=prim->normal();
  
  Ve.normalise();
  Vi.normalise();   
  N.normalise();
  
  //   cout<<"DiffT[lux()] poids = "<<parag.poid();prim->qui();
  scal=Vi.prod_scalaire(N);
  side=face(scal);
  if(side==inf){
    scal=-scal;
    N=-N;
  } 
  //parao.teta=M_PI_2-Macos(-scal);
  parao.teta=Macos(-scal);
  parao.phi=0.0;
  parao.tetaeye=M_PI_2-Macos(N.prod_scalaire(Ve));
  // a revoir si mat anisotrope
  if( Ve.colineaire(N) || Vi.colineaire(N) )
    parao.phieye=0.0;
  else{
    if (Vi.colineaire(Ve))
      parao.phieye= parao.phi;
    else{
      Vecteur Ui,Ue;
      Ui=Vi.prod_vectoriel(N);
      Ue=Ve.prod_vectoriel(N);
      parao.phieye= Ui.prod_scalaire(Ue);
      if(N.prod_scalaire(Ui.prod_vectoriel(Ue))>0.0)
	parao.phieye=2.0*M_PI-parao.phieye;
    }
  }
  return popt(side)->lux(parao)*parag.poid()*fabs(cos(M_PI_2-parao.tetaeye));
}//DiffT::lux()

//-*************** DiffO::interact() ************************

Param_Inter DiffO::interact(Param_Inter &parag,bool inside,int order){
  Transf parao;
  Vecteur rediff;
  double scal=parag.direct().prod_scalaire(prim->normal());

  // cout<<"Diff0[interact()]poid = "<<parag.poid();prim->qui();
  //  cout<<" Normale  : ";prim->normal().show();
  //  cout<<" Incident : ";parag.direct().show();
   if(scal>0.0){
     //  printf(" Si pas infini, bug... Opaque attaque' par de'rrie`re");
     //prim->qui();
     //cout<<" Normale  : ";prim->normal().show();
     //cout<<" Incident : ";parag.direct().show();
     //cout<<" Origine  : ";parag.origin().show();
  inside=true;
  //sortie inchangee ( eg transparence)
  return parag;
  }
  //raus(scal>0.0,"DiffO[interact] Opaque attaque' par de'rrie`re...");
  Einc+=parao.abs=parao.rediffuse=parag.poid();
  parao.teta=Macos(-scal);
  parao.phi=0.0;

  opti->interact(parao);
  //   cout <<"DiffO[interact] teta redif = "<<parao.teta<<" - phi = "<<parao.phi<<endl; 
  rediff=parag.direct(); 
  if( fabs( rediff.prod_scalaire(prim->normal() ) )>1.0-1e-6 ) {
    rediff[0]=primi().sommets(0)[0]-primi().sommets(1)[0];
    rediff[1]=primi().sommets(0)[1]-primi().sommets(1)[1];
    rediff[2]=primi().sommets(0)[2]-primi().sommets(1)[2];
  }
  else
    rediff=rediff.prod_vectoriel(prim->normal());
  rediff.normalise();
    //parao.phi=rambo.tirage()*M_2PI;
  /*
  rediff=rediff.rotate(prim->normal(), parao.phi);
  rediff.normalise();
  rediff=prim->normal().rotate(rediff,(parao.teta<0.0)?M_PI+parao.teta:parao.teta);
  rediff.normalise();
  */
  Vecteur u,v,w,n;
  
  n[0]=n[1]=sin(parao.teta);
  n[0]*=cos(parao.phi);
  n[1]*=sin(parao.phi);
  n[2]=cos(parao.teta);     
  v=rediff;
  w=prim->normal();
  u=v.prod_vectoriel(w);
  for(register int i=0;i<3;i++){
    rediff[i]= u[i]*n[0] + v[i]*n[1] + w[i]*n[2] ;
  }
  
  //rediff=prim->normal().rotate(rediff,M_PI_2);
  //cout<<"DiffO::interact() : teta = "<<parao.teta<<endl;

  /*
    rediff[0]=cos(parao.phi)*sin(parao.teta);
    rediff[1]=sin(parao.phi)*sin(parao.teta);
    rediff[2]=cos(parao.teta);     
    */
  
  parag.change_direction(rediff);
  // if(parao.rediffuse>0.5) cout<<"interact poid rediffuse = "<<parao.rediffuse<<endl;

  //maj de l'energie du diffuseur
  if(order==0){
     B1st+=parao.rediffuse;
      Eabs[0]+=parao.abs;
  }else
     Eabs[1]+=parao.abs;
  radio+=parao.rediffuse;
  nk++;//cout<<"nk++\n";
  parag.change_poids(parao.rediffuse);
  // cout<<"DiffO[interact()] abs="<<parao.abs<<"=> Eabs ="<<Eabs<<" - Poids = "<<parag.poid()<<endl;
  // cout <<" DiffO[interact()] new direct ";parag.direct().show(); 
  return parag;
}

//-*************** DiffT::interact() ************************

Param_Inter DiffT::interact(Param_Inter& parag,bool inside,int order){
  Transf parao;
  Vecteur rediff,N=prim->normal();
  double Ei,Etrans;
  double scal=parag.direct().prod_scalaire(prim->normal());
  Face side=face(scal);

  inside=inside;
  //cout<<"DiffT[interact()]poid = "<<parag.poid()<<" Face :"<<side;prim->qui();
  Einc[side]+=parao.abs=parao.rediffuse=Ei=parag.poid();
  //parao.teta=M_PI_2-Macos(-scal);
  if(side==inf){
    scal=-scal;
    N=-N;
  }
  parao.teta= Macos(-scal);
  //cout <<"DiffT[interact] teta inc = "<<parao.teta<<endl;
  //cout <<"DiffT[interact] scal = "<<scal<<" normale = ";N.show();fflush(stdout);
  parao.phi=0.0;
  // appel a Actop[side]

  /*Bug virtual MC09/98 Test2 => du a un delete opti; dans ~DiffT() !!!
   printf("DiffT[interact] Appel de opti->koi \n");fflush(stdout);
  opti->koi();
  printf(" Apres Actop(opti)::koi() \n");fflush(stdout);

  printf("DiffT[interact] Appel de popt(%d)->koi() \n",(int)side);fflush(stdout);
  popt(side)->koi();
  printf(" Apres Actop::koi() \n");fflush(stdout);
  

  */ 

  
  
  popt(side)->interact(parao);
  // printf("DiffT[interact] Retour de l'Appel de popt->interact() \n");fflush(stdout);
  //parao.teta*=1.0-2.0*side;
  //cout <<"DiffT[interact] teta redif = "<<parao.teta<<endl;fflush(stdout);

  rediff=parag.direct();
  if( fabs( rediff.prod_scalaire(prim->normal() ) )>1.0-1e-6) {
    rediff[0]=primi().sommets(0)[0]-primi().sommets(1)[0];
    rediff[1]=primi().sommets(0)[1]-primi().sommets(1)[1];
    rediff[2]=primi().sommets(0)[2]-primi().sommets(1)[2];
    rediff.normalise();
  }
  else
    rediff=rediff.prod_vectoriel(N);
  rediff.normalise();
  /* rediff=rediff.rotate(N,parao.phi);
     rediff.normalise();
     //rediff=N.rotate(rediff,-M_PI_2+parao.teta);
     rediff=N.rotate(rediff,(parao.teta<0.0)?M_PI+parao.teta:parao.teta);
     */
  Vecteur u,v,w,n;
  
  n[0]=n[1]=sin(parao.teta);
  n[0]*=cos(parao.phi);
  n[1]*=sin(parao.phi);
  n[2]=cos(parao.teta);     
  v=rediff;
  w=N;
  u=v.prod_vectoriel(w);
  for(char i=0;i<3;i++){
    rediff[i]= u[i]*n[0] + v[i]*n[1] + w[i]*n[2] ;
  }

  parag.change_direction(rediff);

  //maj E du diffuseur
  Ei-=parao.abs;// E rediffusee
  Etrans=Ei-parao.rediffuse;//E transmise
  if(order==0)
    B1st[side]+=parao.rediffuse;
  radio[side]+=parao.rediffuse;  //Erroneous arithmetic operation !
  //if(radio[side]<10) {cout<<"DiffT[interact()]  radio ="<<radio[side]<<endl;show();}
  
  // side=(int)((side+1)%2);
  if(side==sup) side=inf; else side=sup;
  if(order==0){
    B1st[side]+=Etrans;
    Eabs[0]+=parao.abs; 
  }else
    Eabs[1]+=parao.abs;  
  //if(radio[side]<10) cout<<"                   radio ="<<radio[side]<<endl;
  radio[side]+=Etrans;
  nk[side]+=1;//cout<<"nk+=1\n";
  parag.change_poids(Ei);
 
  // cout<<"DiffT[interact()] abs="<<parao.abs<<"=> Eabs ="<<Eabs<<" - Poids = "<<parag.poid()<<endl;
  //parag.direct().show();   
  return parag;
}//DiffT::interact()


//-*************** DiffT::maj_stat() ************************

void DiffT::maj_stat(int nbr){
  Eabs[0]/=(double) nbr;
  Eabs[1]/=(double) nbr;
  radio[0]/=(double) nbr;
  B1st[0]/=(double) nbr;
  Einc[0]/=(double) nbr;
  radio[1]/=(double) nbr;
  B1st[1]/=(double) nbr;
  Einc[1]/=(double) nbr;
  nk[0]/=(double) nbr;
  nk[1]/=(double) nbr;
  E1a[0]+=Eabs[0];E1a[1]+=Eabs[1];
  E1i[0]+=Einc[0];
  E1i[1]+=Einc[1];
  radio1[0]+=radio[0];
  radio1[1]+=radio[1];
  B1st1[0]+=B1st[0];
  B1st1[1]+=B1st[1];
  E2a[0]+=Eabs[0]*Eabs[0]; E2a[1]+=Eabs[1]*Eabs[1];
  E2i[0]+=Einc[0]*Einc[0];
  E2i[1]+=Einc[1]*Einc[1];
  radio2[0]+=radio[0]*radio[0];
  radio2[1]+=radio[1]*radio[1];
  B1st2[0]+=B1st[0]*B1st[0];
  B1st2[1]+=B1st[1]*B1st[1];
  Einc[0]=Einc[1]=radio[0]=radio[1]=B1st[0]=B1st[1]=0.0;
  Eabs[0]=Eabs[1]=0.;
}// DiffT::maj_stat() 


//-*************** Diff0::maj_stat() ************************

void DiffO::maj_stat(int nbr){
  Eabs[0]/=(double) nbr;
  Eabs[1]/=(double) nbr;
  radio/=(double) nbr;
  nk/=(double) nbr;
  B1st/=(double) nbr;
  Einc/=(double) nbr;
  E1a[0]+=Eabs[0]; E1a[1]+=Eabs[1];
  E1i+=Einc;
  radio1+=radio;
  B1st1+=B1st;
  E2a[0]+=Eabs[0]*Eabs[0]; E2a[1]+=Eabs[1]*Eabs[1];
  E2i+=Einc*Einc;
  radio2+=radio*radio;
  B1st2+=B1st*B1st;
  Einc=radio=B1st=0.0;
  Eabs[0]=Eabs[1]=0.;
}// DiffO::maj_stat() 


//-*************** DiffT::calc_stat() ************************

void DiffT::calc_stat(int niter,double Stoit){
  raus((niter<1),"Diffuseur[calc_stat] appel sans specifier le nb d'iterations!");
  // calcul de la moyenne et de la variance
  double n=(double)niter,surf ;
  
  surf=prim->surface()/Stoit;
  // Absorbee
  Eabs[0]=E1a[0]/= n/Stoit;
  Eabs[1]=E1a[1]/= n/Stoit;
  E2a[0]=sqrt((E2a[0]/(n/Stoit/Stoit)-( E1a[0]* E1a[0])) /(n-1.5));
  E2a[1]=sqrt((E2a[1]/(n/Stoit/Stoit)-( E1a[1]* E1a[1])) /(n-1.5));
  //Incident
  Einc[0]=E1i[0]/= n*surf;

  E2i[0]=sqrt((E2i[0]/(n*surf*surf)-( E1i[0]* E1i[0])) /(n-1.5));
  Einc[1]=E1i[1]/= n*surf;
  E2i[1]=sqrt((E2i[1]/(n*surf*surf)-( E1i[1]* E1i[1])) /(n-1.5));
  //radiosite
  radio[0]=radio1[0]/= n*surf;
  radio2[0]=sqrt((radio2[0]/(n*surf*surf)-( radio1[0]* radio1[0])) /(n-1.5));
  radio[1]=radio1[1]/= n*surf;
  radio2[1]=sqrt((radio2[1]/(n*surf*surf)-( radio1[1]* radio1[1])) /(n-1.5));
  //Radiosite 1er ordre
  B1st[0]=B1st1[0]/= n*surf;
  B1st2[0]=sqrt((B1st2[0]/(n*surf*surf)-( B1st1[0]* B1st1[0])) /(n-1.5));
  B1st[1]=B1st1[1]/= n*surf;
  B1st2[1]=sqrt((B1st2[1]/(n*surf*surf)-( B1st1[1]* B1st1[1])) /(n-1.5));
 
}//DiffT::calc_stat()



//-*************** DiffO::calc_stat() ************************

void DiffO::calc_stat(int niter, double Stoit){
  raus((niter<1),"Diffuseur[calc_stat] appel sans specifier le nb d'iterations!");  
  // calcul de la moyenne et de la variance
  double n=(double)niter,surf;
  
  surf=prim->surface()/Stoit;
  //printf(" DiffO::calc_stat : Spoly = %g, surf = %g\n ",prim->surface(),surf);
 // Absorbee
  Eabs[0]=E1a[0]/= n/Stoit;
  Eabs[1]=E1a[1]/= n/Stoit;
  E2a[0]=sqrt((E2a[0]/(n/Stoit/Stoit)-( E1a[0]* E1a[0])) /(n-1.5));
  E2a[1]=sqrt((E2a[1]/(n/Stoit/Stoit)-( E1a[1]* E1a[1])) /(n-1.5));
  //Incident
  Einc=E1i/= n*surf;
  E2i=sqrt((E2i/(n*surf*surf)-( E1i* E1i)) /(n-1.5));
  //radiosite
  radio=radio1/= n*surf;
  radio2=sqrt((radio2/(n*surf*surf)-( radio1* radio1)) /(n-1.5));
  //radiositev 1er ordre
  B1st=B1st1/= n*surf;
  B1st2=sqrt((B1st2/(n*surf*surf)-( B1st1* B1st1)) /(n-1.5));
}//DiffO::calc_stat()


//-*************** DiffT::norm_surf() ************************

void DiffT::norm_surf(int nbr,double Stoit){
  //MC00: ???? pourquoi seult la moyenne est normalisée??
  // calcul de la moyenne et de la variance
  double surf;
  
  surf=prim->surface()*nbr/Stoit;
  Eabs[0]   *= Stoit;
  Eabs[1]   *= Stoit;
  
  Einc[0]/= surf;
  Einc[1]/= surf;
  // cout<<"<norm_surf> radio[0] = "<< radio[0]<<endl;
  radio[0]   /= surf;
  radio[1]   /= surf;
  B1st[0]   /= surf;
  B1st[1]   /= surf;
}//DiffT::norm_surf()



//-*************** DiffO::calc_stat() ************************

void DiffO::norm_surf(int nbr,double Stoit){ 
  // calcul de la moyenne et de la variance
  // calcul de la moyenne et de la variance
  double surf;
  
  surf=prim->surface()*nbr/Stoit;

  Eabs[0]   *= Stoit;
  Eabs[1]   *= Stoit;
  
  Einc/= surf;
  radio   /= surf;
  B1st   /= surf;
}//DiffO::norm_surf()


/**************************************************
 ******************     Diff8     *****************
 **************************************************/

//-*************** Diff8::cstructor *****************************
Diff8::Diff8(Diffuseur * pdif,Vecteur & delta) {
  register int i;
  
  diff=pdif;
  prim=new Polygone(*((Polygone*)&(diff->primi())));
  assert(prim);
  for (i=0; i<prim->nb_sommet(); i++){
    prim->sommets(i)+=delta;
    //prim->sommets(i).show();
    //printf("delta[%d]=%lf - ",i,(double)delta[i]);
  }
  prim->called(diff->primi().name());
  prim->calc_cst_eq();
  //printf("Diff8::Diff8() : called \n");
  
}//Diff8::Diff8()

/**************************************************
 ******************     DiffP     *****************
 **************************************************/

//-*************** DiffP::cstructor *****************************
DiffP::DiffP(Primitive* pprim) {
  prim=pprim;
  Z=prim->sommets(0)[2];
  CO(printf("DiffP(%5.2f) :",Z);prim->show();)
  fluxm.alloue(2,ORDRE_MAX,4,18); //inf/sup, ordre, azimut, teta 
  phim.alloue(2,ORDRE_MAX);//inf/sup,ordre
  fluxm.maj(0.);
  phim.maj(0.);
  pene=0.;
}//DiffP::DiffP()
DiffP::~DiffP(){
  fluxm.free();
  phim.free();
  delete prim;
}

Param_Inter DiffP::interact(Param_Inter &parag,bool inside,int ordres){
  double theta,phi,cs,cg,sinphi,weight=parag.poid();
  char jj,phid,order,tetad;
  // Vecteur & delta=parag.direct();
  double delta[3];

  inside=inside;
  ordres=ordres;

  delta[0]=parag.direct()[0];
  delta[1]=parag.direct()[1];
  delta[2]=parag.direct()[2];

  theta=Macos(fabs(delta[2]));
  cs=cos(theta)*sin(theta);
  if(cs==0)
    cg=1e12;
  else
    cg=72./cs;
  tetad=(char) (theta*180.0/M_PI)/5;
  // printf("; deltaZ=%g,Macos(z)=%g,theta=%.2lf, tetad=%d\n",fabs(delta[2]), Macos(fabs(delta[2])),theta ,tetad); fflush(stdout);
  jj=(delta[2]>0)?0:1;
  if(theta!=0.0) {
    phi=Macos(delta[0]/sqrt(delta[0]*delta[0]+(delta[1]*delta[1])));
    if(delta[1]<0.0)  phi=M_2PI-phi;
  }
  else
    phi=M_2PI*drand48();
  // par rapport au soleil
  phi-=phisun;
  sinphi=sin(phi);
  phi=Macos(cos(phi));
  phi=(sinphi<0)? M_2PI-phi : phi;
  phi*=180.0/M_PI;
  phid=4;
  if((phi<=5) || (phi>=355)) phid=0;
  if(fabs(phi-180.0)<=5)     phid=1;
  if(fabs(phi-90.0)<=5)      phid=2;
  if(fabs(phi-270.0)<=5)     phid=3;
 
  order=(parag.ordre<ORDRE_MAX)?parag.ordre-1:ORDRE_MAX-1;
  if(order==-1){
    pene+=weight;
   CO(printf("DiffP[interact] : pene(%5.2f)= %g\n",Z,pene);fflush(stdout);)
     //  printf("DiffP[interact] : jj=%d,order=%d,phid=%d,tetad=%d :: %.0g\n",(int)jj,(int) order,(int)phid,(int)tetad,prim->name());fflush(stdout);
     //printf("%d,",Z);
  }
  else{
    phim(jj,order)+= weight;
    if(phid<4) {//remplissage matrice flux directionnel
      // printf("DiffP[interact] : jj=%d,order=%d,phid=%d,tetad=%d :: %.0g\n",(int)jj,(int) order,(int)phid,(int)tetad,prim->name());fflush(stdout);
      fluxm(jj,order,phid,tetad)+= weight;// a voir *cg;
    }//if directionnel
  }//if ordre
  // parag.direct().show();printf(" DiffP::interact : FIN\n");
  return parag;
}//DiffP::interact()


/**************************************************
 ******************     DiffP2     *****************
 **************************************************/

//-*************** DiffP2::cstructor *****************************
DiffP2::DiffP2(Primitive* pprim) {
  prim=pprim;
  No=1;
  E = new double[2];
  E[0]=E[1]=0.;
}//DiffP2::DiffP2()

DiffP2::DiffP2(Primitive* pprim,int nbo) {
  int i;

  prim=pprim;
  No=nbo;
  E = new double[No];
  for(i=0;i<No;i++)
    E[i]=0.;
}//DiffP2::DiffP2()

DiffP2::~DiffP2(){
  delete prim;
  delete E;
}

Param_Inter DiffP2::interact(Param_Inter &parag,bool inside,int ordres){
  double weight,scal;

  inside=inside;
  ordres=ordres;
  
  scal=parag.direct().prod_scalaire(prim->normal());
  weight=parag.poid();
  // cout<<"DiffP2[interact()]poid = "<<parag.poid();prim->qui();
  //  cout<<" Normale  : ";prim->normal().show();
  //  cout<<" Incident : ";parag.direct().show();
  //  inside=true;
  if(scal<0.0){
    //  inside=false;
  if(parag.ordre<No){
      E[parag.ordre]+=weight;
    }
    else{
      E[No-1]+=weight;
    }//if order >= No
  }//if Attaque par dessus
  // parag.direct().show();printf(" DiffP2::interact : FIN\n");
  return parag;
}//DiffP2::interact()












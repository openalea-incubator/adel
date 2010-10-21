//\ps -o comm,vsize | grep Radiox

#include <iostream>
using namespace std ;

#include <cstdio>
#include <cmath>
#include <ctime>


#include "Mmath.h"
#include "canopy.h"
#include "outils.h"

//#define IMGTEST
#ifdef IMGTEST
#include "image.h"
#endif

#define _NFF
#include "ff.h"
#ifdef _HD
#include "bzh.h"
#endif
/* Pour quoi ne pas faire une class NFF qui serait dans Canopy? */
//#define _KONTAC
#ifdef _KONTAC
#include "image.h"
static Tabdyn<char,2> tabnc;//nb de contact
#endif

//dimension Buffer
unsigned int NB=128; //plus simple si multiple de 4
static  unsigned int MB=128;
static  unsigned int SPEED_INC=25;
//Prototype et declaration de variable en static pour opti
/*
typedef struct {
  double A;
  double B;
  double C;
} MPE;
*/
// Chronometrage
time_t Tnff=0, Tproj=0, dTnff, dTproj,now;

//Variables globales

static double dFF,mpe[3][4],denv,neginvR;
static double *mpei;
static double &mpe_00=mpe[0][0],&mpe_01=mpe[0][1],&mpe_02=mpe[0][2],&mpe_03=mpe[0][3];
static double &mpe_10=mpe[1][0],&mpe_11=mpe[1][1],&mpe_12=mpe[1][2],&mpe_13=mpe[1][3];
static double &mpe_20=mpe[2][0],&mpe_21=mpe[2][1],&mpe_22=mpe[2][2],&mpe_23=mpe[2][3];

static double RE[4],REt[4];
//variable propre a diffR
static bool transp,acv;
static Vecteur u,v,w;
static Point Pp[3],Ps[3],O,P,mid[2],So[3];
static double M[3][3],b[3];
static Diffuseur * receiver;
static unsigned long int ff_cum,vfar_cum,far_cum,glop_cum;
static double nbFF=0;
//zbuff
//pd static  et template avec le cstructeur
static Tabdyn<Diffuseur *,2> DbuffS,DbuffI;
static Tabdyn<float,2> ZbuffS,ZbuffI;
static double* rhotab;//H : image carree NBxNB
static double* i2stab;//tabule le chgt de coord
static Tabdyn<double,2> Tenv;
static double dzc;//hauteur d'une couche de SAIL
static int Nc;//nombre de couches de SAIL
static double *Eclt;


//Declaration propre a la version HD
#ifdef _HD
// Bug en lecture binaire sur SGI 64
//static unsigned long int nb_prim,nb_face,diag_idx,ff_idx;
static int nb_prim,nb_face,diag_idx,ff_idx;

static Tabdyn<int,1>diag;
static Tabdyn<int,1>ligne;
static Tabdyn<double,2>bfc;
#endif
//fin des declarations

#define PAUSE(msg)  printf(msg);printf("- Taper la touche Any");getchar();

static inline void glob2loc(Point &P,Point &Pl) {
  //chgt de repere x,y,y -> u,v,w
  Vecteur T(P);
  Vecteur translat(O);
  
  T=T-translat;
  Pl[0]=T.prod_scalaire(u);
  Pl[1]=T.prod_scalaire(v);
  Pl[2]=T.prod_scalaire(w);
}//glob2loc

static inline void  loc2sph(Point &P,Point &Pl) {
  //chgt de repere x,y,y -> u,v,w
  Point Pt(P);
  Pl[2]=Msqrt(Pt[0]*Pt[0]+Pt[1]*Pt[1]+Pt[2]*Pt[2]);
  Pl[0]=Pt[0]/Pl[2];
  Pl[1]=Pt[1]/Pl[2];
}//glob2loc
static inline void  sph2loc(Point &P,Point &Pl) {
  //chgt de repere u,v,w -> x,y,z
  Point Pt(P);
  Pl[2]=Pt[2]*Msqrt(1.0-Pt[0]*Pt[0]-Pt[1]*Pt[1]);
  Pl[0]=Pt[0]*Pt[2];
  Pl[1]=Pt[1]*Pt[2];
}//glob2loc

static inline int sph2img(double u) {
  return (int)(((u+1)*NB-1)/2.0);
}
static inline double img2sph(int i) {
  return i2stab[i];
}
static inline double sqrtab(int i,int j) {
  int t,N_2;

  N_2=(double)(NB/2.0);
  i=(i<N_2)?i:NB-i-1;
  j=(j<N_2)?j:NB-j-1;
  // i=N_2-i-1;
  //j=N_2-j-1;
  if(j>i) {
    t=i;
    i=j;
    j=t;
  }
  return rhotab[(int)(i+j*(N_2-(j+1)/2.0))];   
}
static inline double env(int i,int j,double &denv,signed char tr) {
  if(Nc==0)
    return 0.0;
  else {
    register int t,N_2;
    Vecteur dir;
    Point M;
    N_2=(double)(NB/2.0);
    //calcul de la direction (i,j) dans le rep. local
    M[0]=img2sph(i);
    M[1]=img2sph(j);
    M[2]=sqrtab(i,j);
    //local -> global
    //verifier un jour le calcul de cette direction
    for(t=0;t<3;t++) 
      dir [t]=M[0]*u[t]+M[1]*v[t]+M[2]*w[t]*tr;
    // dir [t]=M[0]*u[t]*tr+M[1]*v[t]*tr+M[2]*w[t]*tr;
    M=O+dir*denv;
    //printf("=>env(%d,%d,%d) :",i,j,tr);
    //dir.show(); M.show();
    //calcul de la contrib envt de cette direction
    //printf(" : env() : %lf <? %lf <? %lf \n",Tenv[0][2], M[2],Tenv[Nc][2]); 
    if(M[2]<=Tenv(0,2) ||M[2]>Tenv(Nc,2))
      return 0.0;
    else {
      double layer,Lm,transm;
      unsigned char up,il;
      up=(dir[2]>0)? 0:1;
      layer=(M[2]-Tenv(0,2))/dzc;
      il=(int)layer;
      layer-=il;
      
      //printf("M = (%lf,%lf,%lf) - up=%d - il = %d - layer=%lf \n",M[0],M[1],M[2],up,il,layer);
      if(layer==0.0){
	transm=Tenv(il,up);
#ifdef _HD
	bfc(up+((tr<0)?2:0),il)+=1;
#endif
      }
      else {
	transm=((layer*Tenv(il,up))+(1-layer)*Tenv(il+1,up));
#ifdef _HD
	bfc(up+((tr<0)?2:0),il)+=layer;
	bfc(up+((tr<0)?2:0),il+1)+=(1.-layer);
#endif
      }
      //printf("... transm=%lf x \n",transm);
      // Correction en debug de caribu MCdec2003
      // MCsail produit des flux et non plus des valeurs relatives
      //return Eclt[0]*transm;//a PI pres (cf. NFF())
      return transm;//a PI pres (cf. NFF())
    }
  }//fielse Nc==0
}
static bool init_MPE(int &i,int&j) {
  double del[3];
  mpei=mpe[i];
  //cas A!=0
  del[0]=(Pp[i][1]*Pp[j][2])-(Pp[j][1]*Pp[i][2]);
  del[1]=-(Pp[j][0]*Pp[i][2])+(Pp[i][0]*Pp[j][2]);
  del[2]=(Pp[i][0]*Pp[j][1])-(Pp[j][0]*Pp[i][1]);
  //printf("init_MPE(%d,%d) : ",i,j);
  if(del[0]!=0.0) {
    mpei[0]=1;
    mpei[1]=-del[1]/del[0];
    mpei[2]=del[2]/del[0];
    //printf("\n Cas A=1\n");
    return false;
  }
  //cas ou B!=0
  if(del[1]!=0.0) {
    mpei[1]=1;
    mpei[0]=-del[0]/del[1];
    mpei[2]=-del[2]/del[1];
    //printf("\n Cas B=1\n"); 
    return false;
  }
  //cas ou C!=0
  if(del[2]!=0.0) {
    mpei[2]=1;
    mpei[0]=del[0]/del[2];
    mpei[1]=-del[1]/del[2];
    //printf("\n Cas C=1\n");
  return false;
  }
  //cas a la c..
  return true;
}//init_MPE()
static  bool solve3(double M[3][3],double *b) {
  //renvoie true ou false si le systeme a une solution, ecrite dans b
  double D=0.0,x[3]={0,0,0},t[3];
  register int i,j;
  
  D+=M[0][0]*(M[1][1]*M[2][2]-(M[2][1]*M[1][2]));
  D-=M[0][1]*(M[1][0]*M[2][2]-(M[2][0]*M[1][2]));
  D+=M[0][2]*(M[1][0]*M[2][1]-(M[2][0]*M[1][1]));
  if(D==0.0)//cas limite
    return false;
  else {
    for(j=0;j<3;j++) {
      for(i=0;i<3;i++) {
	t[i]=M[i][j];
	M[i][j]=b[i];
      }
      x[j]+=M[0][0]*(M[1][1]*M[2][2]-(M[2][1]*M[1][2]));
      x[j]-=M[0][1]*(M[1][0]*M[2][2]-(M[2][0]*M[1][2]));
      x[j]+=M[0][2]*(M[1][0]*M[2][1]-(M[2][0]*M[1][1]));
      x[j]/=D;
      for(i=0;i<3;i++) 
	M[i][j]=t[i];
    }//for j
    for(i=0;i<3;i++) 
      b[i]=x[i];
    return true;
  }
}//solve3

static inline void Pz0(Point &M, Point &A, Point &B) {
  double k;

  k = B[2]/(B[2]-A[2]);
  M[0] = k*A[0] + (1.-k)*B[0];
  M[1] = k*A[1] + (1.-k)*B[1];
  M[2] = 0.0;
}//Pz0()
static inline void centre(Point &M, Point &A, Point &B, Point &C) {

  M[0] = (A[0] + B[0] + C[0])/3.;
  M[1] = (A[1] + B[1] + C[1])/3.;
  M[2] = (A[2] + B[2] + C[2])/3.; 
}//Pz0()

/****************************************************/
/**********    Fonctions exteriorisables    *********/
/****************************************************/

//init_NFF()
#ifdef _HD
void   init_NFF(char * EnvName,double *Esource,int nbpr,int nbf,double &Rsph,bool bias) {
#else
void   init_NFF(char * EnvName,double *Esource,double &Rsph,bool bias) {
#endif
  //precalcul pour les dist. et les dFF
  int a;
  register int i,ii,j,N_2;
  double v2,u2;

  denv=Rsph;
  neginvR=-1./denv;
  acv=bias;
  //alloc des Tabdyn pas possible a la declaration cause ?
  if(verbose)
    Ferr <<"*** Resolution du disque de projection :  "  
	 << NB<<"x"  << NB<<'\n' ;
  MB=NB;
  DbuffS.alloue(NB,MB);
  DbuffI.alloue(NB,MB);
  ZbuffS.alloue(NB,MB);
  ZbuffI.alloue(NB,MB);
#ifdef _KONTAC
   tabnc.alloue(NB,MB);
   tabnc.maj(0);
#endif
  int dimension;
  dimension=(int)(NB/4.0*(NB/2.0+1));
  rhotab=new double[dimension];
  i2stab=new double[NB];
  ff_cum=vfar_cum=far_cum=glop_cum=0;
#ifdef _HD
  nb_prim=nbpr;
  nb_face=nbf;
  diag_idx=0;
  diag.alloue(nb_prim+1);
  diag(0)=1;
  ligne.alloue(nb_face);
#endif
  for(i=0;i<NB;i++) {
    i2stab[i]=(2*i+1)/(double)NB-1.0;
    // printf("\t i2stab[%d]=%lf\n",i,i2stab[i]);
  }
  N_2=(int)(NB/2.0);
  ii=0;
  for(j=0;j<N_2;j++) {
    v2=img2sph(j);
    v2*=v2;
    for(i=j;i<N_2;i++) {
      u2=img2sph(i);
      u2*=u2;
      //a=i+j*(N_2-(j+1)/2.0);
      //if(a!=ii) printf("init rhotab ya 1cdp!\jn");
      if((u2+v2)<1.0)
	rhotab[ii]=Msqrt(1-u2-v2);
      else
	rhotab[ii]=-1;
      ii++;
    }//for i
  }//for j
  //Calcul de dFF
  dFF=4/(double)(NB*MB)/M_PI;
  /*for(i=0;i<NB;i++) 
    for(j=0;j<NB;j++) 
      if(sqrtab(i,j)>=0)//cas ou je suis dans le disque de projection
	dFF+=1;
  dFF=1./dFF;//moins bon sur le cube
  */
  //Chargement du flux moyen
  Eclt=Esource;
  if(EnvName==NULL)
    Nc=0;
  else {
    FILE *fenv;
    fenv=fopen(EnvName,"r");
    if(fenv==NULL){
      fflush(stdout);
      Ferr <<"<!> FF.C <init_NFF> Unable to open the file "  << EnvName<<'\n' ;
      fflush(stderr);
      exit(15);
    }
    fscanf(fenv,"%d %lf",&Nc,&dzc);
    if(verbose>3) printf("Nc = %d - dz = %lf\n",Nc,dzc);
#ifdef _HD    
    bfc.alloue(4,Nc+1);
    FILE* ffb;
    if(verbose>5) 
      Ferr<<"FF.cpp: init_NFF(): file "<<pcBfName<<" open for writing\n"; 
    ffb=fopen(pcBfName,"ab");
    fwrite(&Nc,sizeof(int),1,ffb);
    fwrite(&nb_prim,sizeof(int),1,ffb); 
    fclose(ffb);
#endif
    Tenv.alloue(Nc+1,3);
    for(i=0;i<=Nc;i++) {
      fscanf(fenv,"%lf %lf %lf ",&(Tenv(i,2)),&(Tenv(i,0)),&(Tenv(i,1)));
      if(verbose>1) printf("-> z= %lf : trans[%d] = %lf - ref[%d] = %lf\n",Tenv(i,2),i,Tenv(i,0),i,Tenv(i,1));
    }
    fclose(fenv);
  }
  //printf("***  SPEED_INC = ");scanf("%d",&SPEED_INC);
  //printf(" speed_inc =%d\n",SPEED_INC);
}//init_NFF
//init_proj (en fct de diffR)

#ifdef _HD
//Anti Doublon System
int ADS_val(int n) {
  return ligne(n);
}
void ADS_maj(int n,int val) {
  ligne(n)=val;
}
#endif

void init_proj(Diffuseur * diffR ) {
#ifdef _HD
  ligne.maj(0);
#endif
  transp=!diffR->isopaque();
  //creation du repere local (u,v,w)
  //diffR->togle_face();
  O=diffR->centre();
  //printf("FF:init_proj() Oz=%lf\n",O[2]);
  u.formation_vecteur(diffR->primi().sommets(1),diffR->primi().sommets(0));
  u.normalise();
  w=diffR->normal();
  v=w.prod_vectoriel(u);
  //debug Bfar
  //printf("FF:init_proj() (u,v,w) = \n");
  //u.show(); v.show();w.show();
  receiver=diffR;
  //init des Buffers
  DbuffS.maj(NULL);
  ZbuffS.maj(99999999999.9);
  if(transp) {
    DbuffI.maj(NULL);
    ZbuffI.maj(999999999999.9);
  }
}//init_proj()

//Projection orthospherique (cas du triangle)
//variable locale (globale ) a proj_ortho
static double x,y,extrem,div_extr,den_extr,pos_extr,speed_B[4];
void proj_ortho(Diffuseur* E,reel *inc) {
  int i,j,i1,i2,jm,j0;
  signed char pasglop=0,in=0,ins=0,nbd,face,iz0[2]={-1,-1},izm[2]={-1,-1},izp[2]={-1,-1},t;
  //evite les calcul de l'opaque vu par derriee (tige)
  //printf(" Ok guy, I'm in proj_ortho()\n");

  //Chrono
  time(&dTproj);
  
  if(E->isopaque()) {
    Vecteur er;
    er.formation_vecteur(O,E->centre());
    if(er.prod_scalaire(E->normal())>0) {
      if(verbose>1) printf("*** pas vu : opak par derriere\n");
      return;
    }
  }
  Point B;
  //projette E sur R
  //printf("E = %ld\n",E->num());
  for(i=0;i<3;i++) {
    P=E->primi().sommets(i);
    P[0]+=inc[0];
    P[1]+=inc[1];
    glob2loc(P,Pp[i]);
    // cout<<" Pp["<<i<<"] :";Pp[i].show();
    if(Pp[i][2]>1e-6) {
      pasglop++;
      if(izp[0]==-1)
	izp[0]=i;
      else
	izp[1]=i;
    }
    else
      if(Pp[i][2]<-1e-6){
	pasglop--;
	if(izm[0]==-1)
	  izm[0]=i;
	else
	  izm[1]=i;
      }
      else {
	pasglop+=4;
	if(iz0[0]==-1)
	  iz0[0]=i;
	else
	  iz0[1]=i;
      }
     }//for i sommet
  /*  B[i]=(Pp[0][i]+Pp[1][i]+Pp[2][i])/3.0;
      }//for i sommet
      loc2sph(B,B);
      B.show();
      i=sph2img(B[0]);
      j=sph2img(B[1]);
      printf("Bi=(%d,%d) - Bl=(%lf,%lf)\n",i,j,img2sph(i),img2sph(j));
      */
  //Traitement des triangles a` cheval sur le plan Z=0
  // cout<<" pasglop = "<<(int)pasglop<<endl;
  switch (pasglop) {
  case 4 : // subdiv en 2 triangles
    nbd=2;
    mid[0]=Pp[iz0[0]];
    So[0]=Pp[izp[0]];
    So[1]=Pp[izm[0]];
    Pz0(mid[1], So[0], So[1]);
    izm[0]=0;//sert  pour la transpa : 0 Z+ , 1 : Z-
    break;
  case 1 : //  pointe en Z- : subdiv en 3 triangles
    nbd=3;
    So[0]=Pp[izm[0]];
    So[1]=Pp[izp[0]];
    So[2]=Pp[izp[1]];
    Pz0(mid[0], So[0], So[1]);
    Pz0(mid[1], So[0], So[2]);
    izm[0]=1;//sert  pour la transparence
    break;
  case -1 : //  pointe en Z+ : subdiv en 3 triangles
    nbd=3;
    So[0]=Pp[izp[0]];
    So[1]=Pp[izm[0]];
    So[2]=Pp[izm[1]];
    Pz0(mid[0], So[0], So[1]);
    Pz0(mid[1], So[0], So[2]);
    izm[0]=0;//sert  pour la transparence
    break;
  default : //cas ou pas besoin de subdiviser
    //printf("==> pas de subdiv \n");
    nbd=1;
    izm[0]=(pasglop==3 || pasglop==6 || pasglop==9 )?0 : 1;
    mid[0]= Pp[0]; 
    mid[1]= Pp[1]; 
    So[0] = Pp[2]; 
  }

   for(t=0;t<nbd;t++) {//boucle de subdivion du triangle
    switch (t) {
    case 0 :
      Pp[0]=mid[0]; Pp[1]=mid[1]; Pp[2]=So[0]; break;
    case 1 :
      Pp[0]=mid[1]; Pp[1]=mid[0]; Pp[2]=So[1]; izm[0]=(izm[0]+1)%2; break;
    case 2 :
      Pp[0]=mid[1]; Pp[1]=So[1]; Pp[2]=So[2]; break;
    }//switch t
    
    //Symetrie /plan Z=0 tq validite des eq.
    if(izm[0]==1) {
      Pp[0][2]=fabs(Pp[0][2]);
      Pp[1][2]=fabs(Pp[1][2]);
      Pp[2][2]=fabs(Pp[2][2]);
    }//triangle en dessous 
   
    //Traitement des cas-limite
    signed char caxe[2]={0,0};
    for(i=0;i<3;i++) {
      if((Pp[i][0]==0)&&(Pp[i][1]==0.))
	Pp[i][0]=Pp[i][1]=0.0000001;
      else
	for(j=0;j<2;j++) {
	  if(Pp[i][j]==0) 
	    if(caxe[j]!=0)
	      Pp[i][j]= Pp[caxe[j]-1][j]=0.0000001;
	    else
	      caxe[j]=i+1;
	}
      if(Pp[i][2]<1e-6) 
	Pp[i][2]+=0.0000001;
    //Remplissage de la mat. du plan du triangle
    //printf("Pp[%d] : ",i); 
    for(j=0;j<3;j++) {
      //printf("%lf  ",Pp[i][j]);
      M[i][j]=Pp[i][j];
    }
    //printf("\n");
      b[i]=-1;//b contiendra les coeff de l'equation du plan du triangle E
      //printf("PROJ_ORHTO() : Pp[%d] = %f, %f, %f\n",i,Pp[i][0],Pp[i][1],Pp[i][2]);
    }//for i

     //printf("PROJ_ORHTO() : diffE no. %d : sens = %d\n",E->num(),izm[0]);
    if(((izm[0]==0) ||(transp && (izm[0]==1))) && solve3(M,b)) {
      //cas ou E n'intersecte pas le plan de projection
      //projection ortho (C. Renaud- LIL)
      //Calcul du centre du triangle P
      centre(P,Pp[0],Pp[1],Pp[2]);
      //printf("PROJ_ORHTO() : P = %f, %f, %f\n",P[0],P[1],P[2]);
      loc2sph(P,P);
      //Warning this is a bug ! :extrem=sqrtab(sph2img(P[0]),sph2img(P[1]));
      extrem=Msqrt(1-P[0]*P[0]-P[1]*P[1]);
      //for(i=0;i<3;i++) cout<<"b["<<i<<"]="<<b[i]<<", ";
      for(i=0;i<3;i++) {//loop arete
	mpei=mpe[i];
	j=(i+1)%3;
	if(init_MPE(i,j)){
	  Vecteur er;
	  er.formation_vecteur(O,E->centre());
	  Ferr <<" Pas possible de calculer init_MPE : triangle no. "  
	       << E->num()<<"  biz? en le projetant sur "  << receiver->num()
	       <<", qui sont a une distance de "  << er.norme()<<'\n' ;
	  //exit(0);
	  return;
	}
	loc2sph(Pp[i],Ps[i]);
	//printf("coord sph : Ps[%d] = %f, %f, %f\n",i,Ps[i][0],Ps[i][1],Ps[i][2]);//valeur de MPE au centre du TRiangle
	mpei[3]=mpei[0]*P[0]+mpei[1]*P[1]+mpei[2]*extrem;
	//printf("MPE[%d][3] = %lf\n",i, mpei[3]);
      }//lopp arete
      //calcul du rectangle englobant
      // 0 left - 1 right - 2 down - 3 up (codage des RE)
      RE[0]=RE[2]=1.5;
      RE[1]=RE[3]=-1.5;
      double maxP,minP,tab_mpe2[4];
      signed char kk;
      for(i=0;i<3;i++) {//loop arete
	mpei=mpe[i];
	tab_mpe2[0]=mpei[0]*mpei[0];
	tab_mpe2[1]=mpei[1]*mpei[1];
	tab_mpe2[2]=mpei[2]*mpei[2];
	tab_mpe2[3]=mpei[0]*mpei[1];
	//div_extr=mpei[0]*mpei[0]+mpei[1]*mpei[1]+mpei[2]*mpei[2];
	div_extr=tab_mpe2[0] +  tab_mpe2[1] + tab_mpe2[2];
	// cout<<"a"<<i;
	for(j=0;j<2;j++) {
	  // cout<<" - axe "<<j;
	  //den_extr=mpei[j]*mpei[j]+mpei[2]*mpei[2];
	  den_extr =  tab_mpe2[j] + tab_mpe2[2];
	  extrem=Msqrt(den_extr/div_extr);
	  if(den_extr==0.) {
	    Vecteur er;
	    er.formation_vecteur(O,E->centre());
	    Ferr <<"=> ARRET proj_ortho() --->arete no. "  << i
		 <<" -  den_extr=0!\n" ;
	    Ferr <<"   triangle no. "<< E->num()<<"  biz? en le projetant sur "
		 << receiver->num()<<", qui sont a une distance de "  
		 << er.norme()<<'\n' ;
	    for(char jj=0;jj<3;jj++)
	      Ferr <<"mpe["  << jj<<"]= "  << mpei[jj]<<" - " ;
	    Ferr <<"\n-----------------------------------------"<<'\n' ;
	    //exit(0);
	    return;
	  }
	  pos_extr=(-tab_mpe2[3]*extrem)/den_extr;
	  jm=3-j*2;
	  if(Ps[(i+1)%3][j]>Ps[i][j]) {
	    i1=i;i2=(i+1)%3;
	  }
	  else {
	    i2=i;i1=(i+1)%3;
	  }
	  
	  kk=(j+1)%2;
	  maxP=max(Ps[i1][kk],Ps[i2][kk]);
	  minP=min(Ps[i1][kk],Ps[i2][kk]);
	  //printf("\nABSI :  min=%f, extr=%lf, max=%f\n",Ps[i1][j],pos_extr,Ps[i2][j]);
	  //printf("ORDO :  min=%lf, extr=%lf, max=%lf\n",Ps[i1][kk],extrem,Ps[i2][kk]);
	  
	  if((Ps[i1][j]>0)^(Ps[i2][j]>0) == 1) {
	    double a,ya;
	    a=(pos_extr-Ps[i1][j])/(Ps[i2][j]-Ps[i1][j]);
	    ya=(1-a)*Ps[i1][kk]+a*Ps[i2][kk];
	    if((ya>0) ^ (extrem>0) == 1) {
	      //printf(" pasglop: extrem a la mmln!\n");
	      pos_extr=-pos_extr;
	      extrem=-extrem;
	    }
	  }
	  if((pos_extr>Ps[i1][j])&&(pos_extr<Ps[i2][j])) {
	    //cas ou le max est dans l'interval
	    REt[jm]  =(extrem>maxP)? extrem : maxP;
	    REt[jm-1]=(extrem<minP)? extrem : minP;
	  }
	  else {
	    REt[jm]=maxP;
	    REt[jm-1]=minP;
	  }
	  RE[jm]=max( REt[jm],RE[jm]);
	  jm--;
	  RE[jm]=min( REt[jm],RE[jm]);
	  //jm++;cout<<jm<<": ordo(min)="<<Ps[i1][kk]<<", ordo(max)="<<Ps[i2][kk]<<", ordo(tgte)="<< extrem<<", max="<< REt[jm]<<", min="<<REt[jm-1]<<endl;
	}
      }//loop arete 2
      
      //Remplisage des Buffers proxel a proxel
      i1=sph2img(RE[0]); i1=(i1==0)? 0 : i1-1; 
      i2=sph2img(RE[1]); i2=(i2<NB-1)? i2+1 : i2; 
      j0=sph2img(RE[2]); j0=(j0==0)? 0 : j0-1; 
      jm=sph2img(RE[3]); jm=(jm<NB-1)? jm+1 : jm; 

      //printf("E.Gsph: %lf   %lf   %lf \n",P[0],P[1],P[2]);
      //printf("REsph : %lf ; %lf - %lf ; %lf \n",RE[0],RE[1],RE[2],RE[3]) ;
      //printf("REimg : %d ; %d - %d ; %d \n",i1,i2,j0,jm);
      for(i=i1;i<=i2;i++) {
	x=img2sph(i);
	for(j=j0;j<=jm;j++) {//loop sur les proxels du RE
	  if(sqrtab(i,j)!=-1) {//dans le cercle
	    y=img2sph(j);
	    //debug pour voir le carre englobant : DbuffS(i,j)=receiver;
	    // MC debug 151203 in=3;
	    int in=3;
	    extrem=sqrtab(i,j);//extrem var. tempo pour economiser les appel a sqrtab()
	    // Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
	    in -= (x + mpe_01*y + mpe_02*extrem > 0) ^ (mpe_03 > 0); //eq 2.5
	    //Ferr << __FILE__ " : "<< __LINE__ <<"; "<<(int)(x + mpe_11*y + mpe_12*extrem > 0)<<" "<<(int)(mpe_13 > 0)<<" "<<in<< '\n' ;	  	   
	    in -= (x + mpe_11*y + mpe_12*extrem > 0) ^ (mpe_13 > 0);//mpe_i0=1
	    //Ferr << __FILE__ " : "<< __LINE__ << '\n' ;	  
	    in -= (x + mpe_21*y + mpe_22*extrem > 0) ^ (mpe_23 > 0);
	    //Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
	    if(in==3){//dans le mille!	      
	      //calcul de la distance rho
	      extrem=b[0]*x+b[1]*y+b[2]*sqrtab(i,j);
	      //printf("d2(%d,%d) = %g, neginvR=%g\n",i,j,extrem,neginvR );
	      //maj des buffer
	      if(!(!acv && (extrem > neginvR) )){ //cf. PhD Renaud p.41
		if(izm[0]==0){
		  //printf("rho(%d,%d)=%lf\n",i,j,extrem);
		  if(extrem<ZbuffS(i,j)) {
		    //cas reflechi
		    ZbuffS(i,j)=extrem;
		    DbuffS(i,j)=E;
		  }
		}
		else 
		  if(extrem<ZbuffI(i,j)) {
		    //cas transmis
		    ZbuffI(i,j)=extrem;
		    DbuffI(i,j)=E;
		    //printf("ZBUFF(%d,%d) = %lf\n",i,j,ZbuffS(i,j));
		}
	      }// if a cheval (traitt special si acv==false)
	    }//dans le mille
	  }//dans le cercle 
	}//fin loop RE j
      }//fin loop RE i
    }//if glop
  }//for t (subdiv triangle)
  //printf("proj_ortho() : FIN\n");
  //Chrono stop
  
  //Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
  time(&now);
  dTproj-=now;
  Tproj-=dTproj;
}//proj ortho  

//Nusselt Form-Factor (Renaud, LIFL)

//variable de NFF
static   Diffuseur *E;
static  SPROW	*r_sup,*r_inf;
static   double FFenv,rho[2],tau[2];


#ifdef _HD
void NFF(int i_sup,int i_inf,VEC **Cfar){
  int ii,nnz=0;
#else
void NFF(SPMAT*FF,int &i_sup,int &i_inf,VEC **Cfar){
#endif
  //Rq : 'cause envt on ne peut pas profiter de Ai.Fij=Aj.Fji
  register int i,j,n,idx;
  //printf("NFF() : DEBUT\n");
  //Chrono
  time(&dTnff);
  //init des prp optiques
  rho[0]=receiver->rho();
  if (transp) {
    tau[0]=receiver->tau();
    receiver->togle_face();
    tau[1]=receiver->tau();
    rho[1]=receiver->rho();
    receiver->togle_face();
  }
  //init des lignes
#ifdef _HD
  ligne.maj(0);
  if(Nc>0)
    bfc.maj(0.);
#else
  r_sup = FF->row+ i_sup;
  r_sup->len=0;
  if(transp) {
    r_inf = FF->row+ i_inf;
    r_inf->len=0;
  }
  //Diagonale FF = 1
  sprow_set_val(r_sup,i_sup,1.);
  if(transp)
    sprow_set_val(r_inf,i_inf,1.);
#endif
  //parcours des Buffers
  for(i=0;i<NB;i++) {
    for(j=0;j<NB;j++) {
      if(sqrtab(i,j)>=0){//cas ou je suis dans le disque de projection
	E=DbuffS(i,j);
	//printf("==> NFF(%d,%d) = %d : ",i,j,(E==NULL)?0:1);
	if(E!=NULL) {
	  n=E->num();
	  //face sup en reflectance
	  ff_cum++;
#ifdef _HD
	  ligne(n)++;
#else
	  idx=sprow_idx(r_sup,n);
	  if( idx<0){//on init FF(i,j)
	    sprow_set_val(r_sup,n,dFF);
	    nbFF++;
	  } else//on incremente FF
	    r_sup->elt[idx].val+=dFF;
#endif
	}
	else {
	  FFenv=env(i,j,denv,1)*dFF;//xier par M_PI si luminance
	  //printf("FFenv(sup)=%lf\n",FFenv);
	  Cfar[0]->ve[i_sup]+=FFenv*rho[0]; //face sup en refl
	  far_cum++;
	  if(FFenv>0) vfar_cum++;
	  //printf("\t (%d,%d) en far\n",i,j);
	}
	if(transp) {
	  if(E!=NULL) {
	    //face inf en transmittance
	    ff_cum++;
#ifndef _HD
	    idx=sprow_idx(r_inf,n);
	    if( idx<0){//on init FF(i,j)
	      sprow_set_val(r_inf,n,-dFF);
	      nbFF++;
	    } else//on incremente FF
	      r_inf->elt[idx].val-=dFF;
#endif
	  }
	  else {
	    Cfar[0]->ve[i_inf]+=FFenv*tau[0];//face inf en trans
	    far_cum++; if(FFenv>0) vfar_cum++;
	  }
	  //autre hemisphere
	  E=DbuffI(i,j);
	  if(E!=NULL) {
	    n=E->num();
	    //face inf en reflectance
	    ff_cum++;
	    ff_cum++;
#ifdef _HD
	    ligne(n)--;
#else
	    idx=sprow_idx(r_inf,n);
	    if( idx<0){//on init FF(i,j)
	      sprow_set_val(r_inf,n,dFF);
	      nbFF++;
	    }else//on incremente FF
	      r_inf->elt[idx].val+=dFF;
	    //face sup en transmittance
	    idx=sprow_idx(r_sup,n);
	    if( idx<0){//on init FF(i,j)
	      sprow_set_val(r_sup,n,-dFF);
	      nbFF++;
	    }else//on incremente FF
	      r_sup->elt[idx].val-=dFF;
#endif
	  }
	  else {
	    //contribution du milieu
	    FFenv=env(i,j,denv,-1)*dFF;//pense a xier par rho ou tau
	    Cfar[0]->ve[i_inf]+=FFenv*rho[1]; //face inf en refl
	    Cfar[0]->ve[i_sup]+=FFenv*tau[1]; //face sup en trans 
	    far_cum+=2; if(FFenv>0) vfar_cum+=2; }
	}//if transp
	
      } //if dans le disque
    }//for j
  }//for i
#ifdef _HD
  FILE *ffb;
  //  if(verbose>5)  Ferr<<"FF.cpp: NFF(): file "<<pcNzName<<" open for writing\n"; 
  ffb=fopen(pcNzName,"ab");
  int tamp[2];
  for(i=0;i<nb_face;i++) 
    if(ligne(i)!=0) {
      nnz++;
      tamp[0]=i;
      tamp[1]=ligne(i);
      // if(i_sup/2==29) printf(" rec %d  - emit = %d - iff = %d\n",i_sup/2,i,ligne(i));
      fwrite(tamp,sizeof(int),2,ffb);
      //fprintf(ffa,"%d %d  ",tamp[0],tamp[1]);
    }
  fclose(ffb);
  diag(diag_idx+1)=(int)fabs(diag(diag_idx))+nnz;
  diag_idx++;
  nbFF+=nnz;
  if(transp){
    nbFF+=nnz;
    diag(diag_idx)*=-1;
  }
  //printf("diag(%d) = %d, nnz=%d\n",diag_idx,diag(diag_idx),nnz); fflush(stdout);//debug
  
  //remplissage en append du fichier des coeff de la CL des Bfar
  if(Nc>0){
    ffb=fopen(pcBfName,"ab");
    float clc[2];
    for(i=0;i<=Nc;i++){
      clc[0]=bfc(0,i)*dFF;
      clc[1]=bfc(1,i)*dFF;
      fwrite(clc,sizeof(float),2,ffb);
    }
    if(transp)
      for(i=0;i<=Nc;i++){
	clc[0]=bfc(2,i)*dFF;
	clc[1]=bfc(3,i)*dFF;
	fwrite(clc,sizeof(float),2,ffb);
      } 
    fclose(ffb);
  }
#endif
  time(&now);
  dTnff-=now;
  Tnff-=dTnff;
}//NFF()

void stat_NFF() {
  //affichage des stats
#define EPSILON 1E-6
#define NONZERO(A) if ((A<EPSILON)&&(A>-EPSILON)){A= A>0 ? EPSILON: -EPSILON;}
  double dummy ;

  dummy= (double)(ff_cum+far_cum) ;
  //Ferr<<"dummy= "<<dummy<<'\n';
  NONZERO(dummy);
  //Ferr<<"NONZERO(dummy) = "<<dummy<<'\n';
  dummy= ff_cum/dummy ;
  printf("\n\t%%age interaction proche  = %lf\n", dummy);

  dummy= (double)(ff_cum+vfar_cum) ;
  NONZERO(dummy);
  dummy= vfar_cum/dummy ;
  printf("\n\t%%age de SAIL non nul     = %lf\n", dummy );

  dummy = (double)(ff_cum+far_cum+glop_cum) ;
  NONZERO(dummy);
  dummy =  glop_cum/dummy ;
  printf(  "\t%%age de triangle pasglop = %lf\n", dummy );

  printf( "\tNombre de FF calcules    = %.0lf\n", nbFF);
  printf("\n@ Temps total de proj_ortho() : %d:%d:%d (%ld s)\n",(int)(Tproj/3600),(int)((Tproj%3600)/60),(int)(Tproj%60),Tproj);
  printf("@ Temps total de NFF()        : %d:%d:%d (%ld s)\n\n",(int)(Tnff/3600),(int)((Tnff%3600)/60),(int)(Tnff%60),Tnff);
  fflush(stdout);
#ifdef _HD
  register int i;
  //Ecriture du fichier binaire diag.bzh
  FILE* diagb;
  Ferr<<"FF.cpp: stat_NFF(): file "<<pcDgName<<" open for writing\n"; 
  diagb=fopen(pcDgName,"wb");
  fwrite(&nb_face,sizeof(int),1,diagb);
  fwrite(&nb_prim,sizeof(int),1,diagb);
  fwrite(&dFF,sizeof(double),1,diagb);	  
  for(i=0;i<=diag_idx;i++) 
    fwrite(&(diag(i)),sizeof(int),1,diagb); 
  fclose(diagb);
  //liberez la memoire!
  diag.free();
  ligne.free();
  ZbuffI.free(); ZbuffS.free();
  DbuffI.free(); DbuffS.free();
  delete rhotab; delete i2stab;  
#endif
}//stat_NFF()

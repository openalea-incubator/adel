/*******************************************************************
*                 canopy.C - mike 95                               *
*            classe scene de la radioxity                          *
*              Precalcul pour le solver                            *
********************************************************************/
#include <iostream>	// pour utiliser namespace (compile mieux)
using namespace std ;

#include <cmath>
#include <cstdio>

#include "outils.h"
#include "chrono.h"

#include "canopy.h"
#include "ff.h"
#include "Mmath.h"

//#define FFs
#ifdef FFs
//lib FF de Schroder pour polygone ss occultation
extern "C"{
#include "ff.h"
}
#endif
#define REELLE float

#include "infini.h"
#include "image.h"

#define DEBOG
#include "debug.h"

#define A Pp[i]
#define B Pp[j]
#define C Pp[k]
#define D Pp[3]

typedef struct{
  BSP * box;
  reel inc[2];
}Boxi;

static Chrono chron;
static Chrono tps;
static int proj_cpt;
#define PAUSE(msg)  printf(msg);printf("- Taper la touche Any");getchar();


// Cas de SAIL pur (Dsph==0)
void Canopy::sail_pur(VEC **Cfar,double *Eclt,char* envname){
  int Nc;
  long i_sup=0,i_inf=0;
  Point G;
  Tabdyn<double,2> Tenv;
  FILE *fenv;
  double dzc,layer,transm,nz,Ei[2],E[2],costl,rho[2],tau[2];
  unsigned char up,down,il;
  Diffuseur *diffR;
  
  // Chargement des profils d'eclairement (E+,E-) calcule par SAIL    
  fenv=fopen(envname,"r");
  if(fenv==NULL){
    fflush(stdout);
    Ferr <<"<!> FF.C <init_NFF> Unable to open the file " << envname<<'\n' ;
    exit(4);
  }
  fscanf(fenv,"%d %lf",&Nc,&dzc);
  printf("Canopy::sail_pur(): Nc = %d - dz = %lf\n",Nc,dzc);
  Tenv.alloue(Nc+1,3);
  for(short i=0;i<=Nc;i++) {
    fscanf(fenv,"%lf %lf %lf ",&(Tenv(i,2)),&(Tenv(i,0)),&(Tenv(i,1)));
    if(verbose>1)
      printf("-> z= %lf : trans[%d] = %lf - ref[%d] = %lf\n",Tenv(i,2),i,Tenv(i,0),i,Tenv(i,1));
  }
  fclose(fenv);
  
  // Calcul de Bfar (Boucle sur tousles diffuseurs => maj de Cfar)
  for(Ldiff.debut();! Ldiff.finito();Ldiff.suivant()){
    diffR=Ldiff.contenu();
    diffR->active(0);//active la face sup
    nz = diffR->primi().normal()[2];
    costl=fabs(nz);
    G=diffR->centre(); 
    i_sup=diffR->num();
    rho[0]=diffR->rho();
    if(!diffR->isopaque()) {
      tau[0]=diffR->tau();
      diffR->togle_face();
      i_inf=diffR->num();
      tau[1]=diffR->tau();
      rho[1]=diffR->rho();
      diffR->togle_face();
    }
    layer=(G[2]-Tenv(0,2))/dzc;
    il=(int)layer;
    layer-=il;
    up=(nz>0)? 0:1;
    down=(up+1)%2;
    if(layer==0.0){//eg sol
      Ei[0]=Tenv(il,0);//E-(0)
      Ei[1]=Tenv(il,1);//E+(0)    
    }
    else {
      Ei[0] = layer*Tenv(il,0) + (1-layer)*Tenv(il+1,0); //E-(z)
      Ei[1] = layer*Tenv(il,1) + (1-layer)*Tenv(il+1,1); //E+(z)
    }
    // printf("(%d,%d) layer=%lf, Ei[0]=%lf, Ei[1]=%lf\n",(int)il,(int)(il+1),layer,Ei[0], Ei[1]);
    //Eclairement face sup
    E[0] = .5*(1+costl)*Ei[up] + .5*(1-costl)*Ei[down];// Formules dans (Ross,1981) p.131
    Cfar[0]->ve[i_sup]=rho[0]*E[0];
    if(!diffR->isopaque()) {
      Cfar[0]->ve[i_inf]=tau[0]*E[0];
      E[1]=.5*(1-costl)*Ei[up] + .5*(1+costl)*Ei[down];// Formules dans (Ross,1981) p.131
      Cfar[0]->ve[i_sup]+=tau[1]*E[1];;
      Cfar[0]->ve[i_inf]+=rho[1]*E[1];;
    }// if transp
    // printf("costl=%lf => E[0]=%g,rho[0]=%g,  E[1]=%g\n",costl,E[0],rho[0], E[1]);
    // printf("==> Bsup(%ld)=%g, Binf(%ld)=%g\n",i_sup,Cfar[0]->ve[i_sup],i_inf,Cfar[0]->ve[i_inf]);
  }// for Ldiff  
}//sail_pur()


//++ Prototype pour la projection

void colorie_triangle( void *,void *,REELLE **, const Punkt,const Punkt,const Punkt, int,int, double,double, void (*f)(void *, int, int, void*));
void colorie_capteur(REELLE **Zbuf,Punkt a,Punkt b,Punkt c, int Tx, int Ty,double tx, double ty,int& pB0);

void pt2pkt(Point &, Punkt);
void zproj(void *, int, int, void*);
void zFF(void *, int, int, void*);



//+************ zproj()
void zproj(void * Zprim, int i, int j, void* prim) {
  ((Diffuseur***) Zprim)[i][j] = (Diffuseur *) prim;
}//zproj()

//+************ zFF()
void zFF(void *, int, int, void*) {

}//zFF()

static  int addbox(BSP * box,reel dx,reel dy, Boxi *Tabox,int ind) {
  register int i;
  if(ind!=0) {
    if(box->nb_diffuseur()==0)
       return ind;
    for(i=0;i<ind;i++)
      if(box==Tabox[i].box && Tabox[i].inc[0]==dx && Tabox[i].inc[1]==dy)
      // BUG : if(box==Tabox[i].box && Tabox[ind].inc[0]==dx && Tabox[ind].inc[1]==dy)
      return ind;
  }
  Tabox[ind].box=box;
  Tabox[ind].inc[0]=dx;
  Tabox[ind].inc[1]=dy;
  return ++ind;
}//addbox()

#ifdef _HD
void Canopy::calc_FF_Bfar(
#else
void Canopy::calc_FF_Bfar(SPMAT *FF,
#endif
			  VEC **Cfar,
			  double *Esource,
			  char * envname,
			  bool bias,
			  double denv,
			  int nbsim) {
  register int p,t,b,nb_rec=0,nb_emi=0,nb_test=0,cum_box=0,nb_error=0,nb_ff,nop=0;;
  register unsigned char i,nb_box,pr=0;
  int pmax,tmax,inc[3],Gi[3];
  int n,i_sup,i_inf=0,idx,idxn;
  Diffuseur *diffR, *diffE;
  BSP *box;
  Boxi Tabox[8];
  Point G,T;
  reel move[2];
  double mid,S[3],dGS[3],d2env=denv*denv,ff,dER2; //attention aux tests entre d et d^2
  bool rechauffe,sol,select;
  
  if(verbose>1) {
    Ferr <<"Canopy::calc_FF_Bfar() denv=Rsph="  << denv
	 <<" ; taille_vox = "  << mesh.taille()<<'\n' ;
  }
#ifdef _HD
  init_NFF(envname,Esource,Ldiff.card(),radim,denv,bias);
#else
  SPROW	*r_sup,*r_inf;
  init_NFF(envname,Esource,denv,bias);
#endif
  for(Ldiff.debut();! Ldiff.finito();Ldiff.suivant()){
    nb_ff=0;
    diffR=Ldiff.contenu();
    sol = (diffR->primi().name()==0)? true : false;
    //pmax=diffR->nb_patch();
    // printf("[calc_FF_Bfar] diffR = %d\n",diffR->num());
    //printf("%d ",diffR->num());
    //pr++;
    //if((pr%10)==0) printf("\n");
    //chron.Start();
    //proj_cpt=0;
    //printf("* pmax = %d\n",pmax); 
    //for(p=0;p<pmax;p++) {//cas du diffuseur patche
    //diffR->select_patch(p);
    diffR->active(0);//active la face sup
    init_proj(diffR);
    G=diffR->centre();
    i_sup=diffR->num();
    //printf("\n %d -",i_sup);
#ifndef _HD
    r_sup = FF->row+ i_sup;
#endif
    if(!diffR->isopaque()) {
      diffR->togle_face();
      i_inf=diffR->num();
      // printf(" %d ",i_inf);
#ifndef _HD
      r_inf = FF->row+ i_inf;
#endif
      diffR->togle_face();
    }
    nb_box=0;
    for(i=0;i<3;i++) {
      //Gi[i]=mesh.coord(i,G[i]);
      Gi[i]=(int) (G[i]-mesh.origine()[i])/mesh.taille();
      //printf("i=%d - G[i] = %lf - Gi[i]=%d\n",(int)i,G[i],Gi[i]);
      mid=mesh.milieu_vox(i,Gi[i]);
      if(G[i]<mid && (Gi[i]!=0 || (infty &&(i!=2)))) {
	if(Gi[i]!=0) {
	  inc[i]=-1;
	  if(i!=2) move[i]=0.0;
	}
	else{
	  inc[i]=mesh.nb_voxel(i)-1;
	   if(i!=2)move[i]=-delta[i];
	}
	S[i]=mesh.sommet(i,Gi[i]);
      }
      else
	if(G[i]>mid && (Gi[i]!=mesh.nb_voxel(i)-1|| (infty &&(i!=2))) ) {
	  if(Gi[i]!=mesh.nb_voxel(i)-1) {
	  inc[i]=1;
	  if(i!=2) move[i]=0.0;
	}
	else{
	  inc[i]=-mesh.nb_voxel(i)+1;
	  if(i!=2) move[i]=delta[i];
	}
	  S[i]=mesh.sommet(i,Gi[i]+1);
	}
	else {
	  inc[i]=0;
	  if(i!=2) move[i]=0.0;
	  S[i]=0.0;
	}
      if((inc[i]!=0 ||(infty && mesh.nb_voxel(i)==1 && i!=2 &&move[i]!=0))&& fabs(G[i]-S[i])>denv) {
	 inc[i]=0;
	 if(i!=2) move[i]=0.0;
	//cout<<(int)i<<" - "<<nb_rec<<" - "<<pmax<<"bordel acqueux!\n";
      }
    }//loop sur les axes

    //Ferr << __FILE__ " : "<< __LINE__ << '\n' ;

    nb_box=addbox(mesh(Gi[0],Gi[1],Gi[2]),0,0,Tabox,nb_box);
    char axe[3]= {2,2,2},diag=3;
    //axe -> 0 : x-y, 1 : y-z, 2 : z-x
    if(inc[0]!=0 ||(infty && mesh.nb_voxel(0)==1 &&move[0]!=0 )) {
      nb_box=addbox(mesh(Gi[0]+inc[0],Gi[1],Gi[2]),move[0],0,Tabox,nb_box);
      axe[0]--; axe[2]--; diag--;
    }
    if(inc[1]!=0 ||(infty && mesh.nb_voxel(1)==1  &&move[1]!=0 )) {
      nb_box=addbox(mesh(Gi[0],Gi[1]+inc[1],Gi[2]),0,move[1],Tabox,nb_box);
      axe[0]--; axe[1]--; diag--;
    }
    if(inc[2]!=0) {
      nb_box=addbox(mesh(Gi[0],Gi[1],Gi[2]+inc[2]),0,0,Tabox,nb_box);
      axe[1]--; axe[2]--; diag--;
    }
    // precacul distance d(G,S)
    dGS[0]=(G[0]-S[0])*(G[0]-S[0]);
    dGS[1]=(G[1]-S[1])*(G[1]-S[1]);
    dGS[2]=(G[2]-S[2])*(G[2]-S[2]);
    //test des boites jointives par un axe 
    if(axe[0]==0) {
      //test/ dist(G, axe x-y ie z)
      if( (dGS[0]+dGS[1]) < d2env)
	nb_box=addbox(mesh(Gi[0]+inc[0],Gi[1]+inc[1],Gi[2]),move[0],move[1],Tabox,nb_box);
    }
    if(axe[1]==0) {
      //test/ dist(G, axe y-z ie x)
      if( (dGS[1]+dGS[2]) < d2env)
	nb_box=addbox(mesh(Gi[0],Gi[1]+inc[1],Gi[2]+inc[2]),0,move[1],Tabox,nb_box);
    }
    if(axe[2]==0) {
      //test/ dist(G, axe z-x ie y)
      if( (dGS[2]+dGS[0]) < d2env)
	nb_box=addbox(mesh(Gi[0]+inc[0],Gi[1],Gi[2]+inc[2]),move[0],0,Tabox,nb_box);
    }
    //test de la  boite jointive par le sommet S
    if(diag==0) {
      //test/ dist(G, axe y-z ie x)
      //printf("%lf %lf %lf \n",dGS[0],dGS[1],dGS[2]);
      if( (dGS[0]+dGS[1]+dGS[2]) < d2env){
	//printf("mesh(%d+%d, %d+%d, %d+%d), move(%lf,%lf), mesh[%d,%d,%d] \n",Gi[0],inc[0],Gi[1],inc[1],Gi[2],inc[2],move[0],move[1], mesh.nb_voxel(0), mesh.nb_voxel(1), mesh.nb_voxel(2));fflush(stdout);
	nb_box=addbox(mesh(Gi[0]+inc[0],Gi[1]+inc[1],Gi[2]+inc[2]),move[0],move[1],Tabox,nb_box);
      }
    }
    cum_box+=nb_box;
 
   
    //allez zou, on calcule les FF et les Bfar par projection
    if(verbose==2) {
      tps.Start(); 
      //Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
    }
    for(b=0;b<nb_box;b++) {//loop sur les boites dans denv
      box=Tabox[b].box;
      //printf("Boite no %d/%d  = dx = %f, dy = %f\n",b,nb_box,Tabox[b].inc[0],Tabox[b].inc[1]);
      // Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
      //Ferr << "Boite no "<<b<<'\n';
      //loop sur les diffuseurs de la boite
      for(box->Ldiff.debut();! box->Ldiff.finito();box->Ldiff.suivant()){
	//Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
	diffE=box->Ldiff.contenu();
	if (diffE->isreal() && diffE!=diffR) 
	  //si diffE n'est pas un capteur virtuel et si pas cas Diagonale : FF=1
	  if(!(sol && diffE->primi().name()==0)) {
	    //interaction sol-sol (H sol plan)
	    //cas du diffuseur patche
	    //tmax=diffE->nb_patch();
	    //for(t=0;t<tmax;t++) {
	    //diffE->select_patch(t);
	    
	    nb_test++;
	    //OLD si bias==false, test de distance pas fait et triangle traite
	     //OLDif(!(bias && G.dist2(T)> d2env)) {
	    // si bias==false, test de distance point-triangle
	    // si bias==true , triangle teste que si  G.dist2(T)>d2env est faux
	    T=diffE->centre();
	    T[0]+=Tabox[b].inc[0];
	    T[1]+=Tabox[b].inc[1]; 
	    if(bias){
	      dER2= G.dist2(T);
	      select = dER2<d2env;
	    }else{
              /* version exacte, mais trop lent 
	      // a optimiser !!!
	      char polystr[255],str[100];
	      Polygone *E; // Perte du polymorphisme(a voir)
     
	      sprintf(polystr,"%d ",diffE->primi().nb_sommet());
	      for(signed char p=0;p<diffE->primi().nb_sommet();p++){
		sprintf(str," %g %g %g ",
			diffE->primi().sommets(p)[0]+Tabox[b].inc[0],
			diffE->primi().sommets(p)[1]+Tabox[b].inc[1],
			diffE->primi().sommets(p)[2]);
		strcat(polystr,str);
	      } 
	      //printf("polystr = %s\n",polystr);
	      E = new Polygone(polystr,1,NULL,NULL,false);
	      T=E->centre();
	      dER2=E->distance2_point(G);
	      delete E;
	      */
	      // calcul un sous estimateur de la dist au triangle pour selectionner les projetes
	      Point I;
	      I=G;
	      I[0]-=Tabox[b].inc[0];
	      I[1]-=Tabox[b].inc[1];
	      select=diffE->primi().appart_sphere(I,denv);
	      //printf("denv=%g, select=%d\n",denv,(int) select);//I.show();diffE->show();
	    }
	    if(select){
	      Vecteur dir(G,T);
	      if(dir.prod_scalaire(diffE->normal())!=0){
		nb_ff++;
		diffE->active(dir);
		n=diffE->num();
		//Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
#ifdef _HD
		idx= ADS_val(n);
#else
		idx=(int) sp_get_val(FF,i_sup,n);
#endif
		idxn=(Tabox[b].inc[0]==0 && Tabox[b].inc[1]==0)?2:3;
		//if( sprow_idx(r_sup,n)<0) {
		if(((idx*idxn)/2)%3 == 0){
		  //Anti-Doublon System
		  if(idx+idxn==5) nb_error++;
#ifdef _HD
		  ADS_maj(n,idx+idxn);
#else
		  sp_set_val(FF,i_sup,n,idx+idxn);
#endif
		  // Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
		  proj_ortho(diffE,Tabox[b].inc);
		  //Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
		  nb_emi++; proj_cpt++;
		  
		}//if pas deja traite
	      }//if pas parllele a la direction
	    }//if diffE dans l'envt
	    //}//for nb patch diffE (t)	 
	  }//if pas sol-sol //if diffE != diffR
      }//for box->Ldiff 
    }//jqa nb_box(b) (ou if rechauffe si AntiDoublon)
    nb_rec++;
    //Calcul des FF en fct des Buffers
    if(verbose>3){
      tps.Stop();
      cout<<++nop<<" : proj ortho-sph de "<<nb_ff<<" T en "<<tps<<endl;
      //Bug MC2006: ++nop est fait deux fois!!
      //    Ferr <<++nop<<" : proj ortho-sph de "<<nb_ff<<" T en "<<tps<<'\n';
      Ferr <<nop<<" : proj ortho-sph de "<<nb_ff<<" T en "<<tps<<'\n';
    }
    
#ifdef _HD
    NFF(i_sup,i_inf,Cfar);
#else
    NFF(FF,i_sup,i_inf,Cfar);
#endif    //chron.Stop();

    if(verbose>4){
      Ferr <<__FILE__ <<" : "<< __LINE__<< " : ";
      Ferr<<proj_cpt<<'\n'<<" NFF caculées "<<'\n';
    }

    //}//for nb patch diffR
  }//for Ldiff

  stat_NFF(); 

#define EPSILON 1E-6
#define NONZERO(A) if ((A<EPSILON)&&(A>-EPSILON)){A= A>0 ? EPSILON: -EPSILON;}
  double dummy = (double)nb_rec ;
  NONZERO(dummy) ;
  cout<<"\tmean(diff/patch)     = "<<nb_test/dummy<<"\n";
  cout<<"\tmean(diff_proj/patch)= "<<nb_emi/dummy<<"\n";
  cout<<"\tmean(box/patch)      = "<<cum_box/dummy<<"\n";
  cout<<endl;
  // ,nb_test/(double)nb_rec,nb_emi/(double)nb_rec,cum_box/(double)nb_rec);
  if(infty && denv>=min(delta[0],delta[1]))
    Ferr <<"####>> Nb Erreurs dues a active_face et infini = "
	 << nb_error<<"\n" ;
}//Canopy::calc_FF_Bfar()
#ifndef _HD
#ifdef FFs

void Canopy::calc_FF(SPMAT *FF){;
  Point G,T;
  register unsigned char i,j;
  int n,i_sup,i_inf,idx;
  Diffuseur *diffR, *diffE;
  SPROW	*r_sup,*r_inf;
  double R[3][3], E[3][3],ff;
  for(Ldiff.debut();! Ldiff.finito();Ldiff.suivant()){
    diffR=Ldiff.contenu();
    diffR->active(0);//active la face sup
    i_sup=diffR->num();
    r_sup = FF->row+ i_sup;
    if(!diffR->isopaque()) {
      diffR->togle_face();
      i_inf=diffR->num();
      r_inf = FF->row+ i_inf;
      diffR->togle_face();
    }
     G=T=diffR->centre();
     for(i=0;i<3;i++)
        for(j=00;j<3;j++)
	  R[i][j]=diffR->primi().sommets(i)[j];
    //allez zou, on calcule les FF par la lib Schroeder
    Ldiff.deb_pause();
    for(Ldiff.debut();! Ldiff.finito();Ldiff.suivant()){
      diffE=Ldiff.contenu();
      if (diffE!=diffR) {//si pas cas Diagonale : FF=1
	T=diffE->centre();
	Vecteur dir(G,T);
	diffE->active(dir);
	n=diffE->num();
	for(i=0;i<3;i++)
	  for(j=00;j<3;j++)
	    R[i][j]=diffE->primi().sommets(i)[j];
	ff=FormFactor(R,3, E,3);
	sp_set_val(FF,i_sup,n,ff);
      }//if pas deja traite
      else
	sp_set_val(FF,i_sup,i_sup,1);
    }//for Ldiff(E)   
     Ldiff.fin_pause();
  }//for Ldiff(R)
}//Canopy::calc_FF()
#else
void Canopy::calc_FF(SPMAT *FF){
  Ferr <<"Sorry, but your version desn't have the Schoeder Form Factor library"
    " \n => contact chelle@bcgn.grignon.inra.fr"  << (char)7<<"\n" ;
}//Canopy::calc_FF()
#endif

#endif

void Canopy::projplan(Vecteur &visee,bool infty, double* Bo) {
  register int i,j,k,l,img_surf;
  Point roof[4];
  REELLE **Zbuf,**pZbuf,*ptZ;
  double tx,ty,costeta,Apix;
  Diffuseur ***Zprim,***pZprim;
  l=0;
  //   Tabdyn<REELLE, 2> Zbuf(img->taille(0),img->taille(1));
  //Tabdyn<Diffuseur *, 2> Zprim(img->taille(0),img->taille(1));
  if(verbose>2) printf("%c => projplan() DEBUT res. %d x %d\n%c",7, Timg,Timg);
  Zbuf= new REELLE*[Timg];
  pZbuf=Zbuf;
  img_surf=Timg*Timg;
  for(i=0;i<Timg;i++,pZbuf++) {
    (*pZbuf)=new REELLE[Timg];
  }
  pZbuf=Zbuf;
  for(i=0;i<Timg;i++,pZbuf++) {
    ptZ=*pZbuf;
    for(j=0;j<Timg;j++,ptZ++)
      *ptZ=99999999999.9;
  }
  Zprim=new Diffuseur**[Timg];
  for(i=0,pZprim=Zprim;i<Timg;i++,pZprim++) {
    (*pZprim)=new Diffuseur*[Timg];
    for(j=0;j<Timg;j++,ptZ++)
    (*pZprim)[j]=NULL;
  }
  //&&&&&& ProjPlan() &&&&&&&&&
  //calcul de laposition de l'ecran en fonction des bornes de la scene
  Point Ecran[4];
  double Pts[8][2],du,dv;
  Vecteur u,v,w;
  visee.normalise();
  costeta=-visee[2];
  
  if(fabs(visee[2]+1.0)<1e-6) {  
    Ferr<<"Projplan(): cas de la visee verticale\n";
    Ecran[0][0]=vmin[0];
    Ecran[0][1]=vmax[1];
    Ecran[1][0]=vmin[0];
    Ecran[1][1]=vmin[1];
    Ecran[2][0]=vmax[0];
    Ecran[2][1]=vmin[1];
    Ecran[3][0]=vmax[0];
    Ecran[3][1]=vmax[1];
    Ecran[0][2]=Ecran[1][2]=Ecran[2][2]=Ecran[3][2]=vmax[2]+1;
    du=vmax[0]-vmin[0];
    dv=vmax[1]-vmin[1];
    u=visee;
    v[0]=1; v[1]=0; v[2]=0;
    w[0]=0; w[1]=-1; w[2]=0;
    //cout<<"Base (u,v,w) :\n";u.show(); v.show(); w.show();
  }
  else { //cas  : visee non verticale
    double M[2][2],tmp,maxi[2]={-999999.9,-999999.9},mini[2]={999999.9,999999.9},dz=vmax[2]-vmin[2];
    
    //points a Z=0   
    Ecran[0][0]=Pts[0][0]=vmin[0];
    Ecran[0][1]=Pts[0][1]=vmin[1];
    Ecran[1][0]=Pts[1][0]=vmax[0];
    Ecran[1][1]=Pts[1][1]=vmin[1];
    Ecran[2][0]=Pts[2][0]=vmax[0];
    Ecran[2][1]=Pts[2][1]=vmax[1];
    Ecran[3][0]=Pts[3][0]=vmin[0];
    Ecran[3][1]=Pts[3][1]=vmax[1];
    //projection selon visee des points a z = vmax[2]
    for(i=0;i<4;i++) {
      roof[i]=Ecran[i];
      Ecran[i][2]=vmax[2];
      Ecran[i]=Ecran[i]+visee*fabs(dz/visee[2]);
      Pts[i+4][0]=Ecran[i][0];
      Pts[i+4][1]=Ecran[i][1];
      //printf("projplan() : Ecran[%d][2] = %g\n",i,Ecran[i][2]); 
      Ecran[i][2]=vmin[2];
    }
    roof[0][2]=roof[1][2]=roof[2][2]=roof[3][2]=vmax[2];
    //chgt de repere x,y -> u,v && calculer min et max
	 u=visee;
     u[2]=0;
     u.normalise();
	 v[0]=-u[1];
     v[1]=u[0];
     v[2]=0;
     w=visee.prod_vectoriel(v); 
     w.normalise();
     if(w[2]<0) {
       w=-w;
       v=-v;
     }
     if(verbose>1){cout<<"Base (u,v,w) :\n";u.show(); v.show(); w.show();}
    M[0][0]=u[0] ;
    M[1][0]=v[0]; ;
    M[0][1]=u[1] ;
    M[1][1]=v[1] ;
    for(i=0;i<8;i++) {
      tmp=Pts[i][0];
      for(j=0;j<2;j++) {
	Pts[i][j]=M[j][0]*tmp+M[j][1]*Pts[i][1];
	mini[j]=min(mini[j],Pts[i][j]);
	maxi[j]=max(maxi[j],Pts[i][j]);      
      }
    }
    //coord de la projection de l'ecran sur le sol (Ruv)
    Ecran[0][0]=mini[0];
    Ecran[0][1]=mini[1];
    Ecran[1][0]=maxi[0];
    Ecran[1][1]=mini[1];
    Ecran[2][0]=maxi[0];
    Ecran[2][1]=maxi[1];
    Ecran[3][0]=mini[0];
    Ecran[3][1]=maxi[1];
    //chgt repere uv -> xy
    M[0][0]=u[0] ;
    M[1][0]=u[1]; ;
    M[0][1]=v[0] ;
    M[1][1]=v[1] ;
    for(i=0;i<4;i++) {
      tmp=Ecran[i][0];
      //Ecran[i].show();
      for(j=0;j<2;j++) {
	Ecran[i][j]=M[j][0]*tmp+M[j][1]*Ecran[i][1];
      }
      //Ecran[i].show();
    }
    //redressement de l'ecran tq ortho a visee (Ecran 0 et 3 : u min)
    dv=tmp=(maxi[0]-mini[0])*u.prod_scalaire(w);
    //printf("mini[0] = %g - mini[1] = %g ; maxi[0]= %g - maxi[1] = %g\n",mini[0],mini[1],maxi[0], maxi[1]);
    //printf(" translation de d de la base Ecran : d = %g\n",tmp);
    Ecran[1]=Ecran[0]+w*tmp;
    Ecran[2]=Ecran[3]+w*tmp;
    //translation de l'ecran tq tota la scene soit vue
    tmp=(dz+1)/visee[2];
    for(i=0;i<4;i++)
      Ecran[i]= Ecran[i] + visee*tmp;

    du = maxi[1]-mini[1];
    u=visee;
  }//if visee verticale
 
/* validations   
  //validation geom
  if(false){
  FILE *fcan;
  fcan=fopen("proj.can","w");
  // Bounding Box
  fprintf(fcan,"p 1 1 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",vmin[0],vmin[1],vmin[2], vmax[0],vmin[1],vmin[2], vmax[0],vmin[1],vmax[2], vmin[0],vmin[1],vmax[2]);
  fprintf(fcan,"p 1 2 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",vmax[0],vmin[1],vmin[2], vmax[0],vmax[1],vmin[2], vmax[0],vmax[1],vmax[2], vmax[0],vmin[1],vmax[2]);
  fprintf(fcan,"p 1 3 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",vmax[0],vmax[1],vmin[2], vmin[0],vmax[1],vmin[2], vmin[0],vmax[1],vmax[2], vmax[0],vmax[1],vmax[2]);
  fprintf(fcan,"p 1 4 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",vmin[0],vmax[1],vmin[2], vmin[0],vmin[1],vmin[2], vmin[0],vmin[1],vmax[2], vmin[0],vmax[1],vmax[2]);
  fprintf(fcan,"p 1 5 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",vmin[0],vmin[1],vmax[2], vmax[0],vmin[1],vmax[2], vmax[0],vmax[1],vmax[2], vmin[0],vmax[1],vmax[2]);
  //Ecran de projection
 fprintf(fcan,"p 1 0 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",Ecran[0][0],Ecran[0][1],Ecran[0][2], Ecran[1][0],Ecran[1][1],Ecran[1][2], Ecran[2][0],Ecran[2][1],Ecran[2][2], Ecran[3][0],Ecran[3][1],Ecran[3][2]);
	fclose(fcan);
  }
*/

  //Projection prallele a visee sur Ecran
  Point Pp[4];
  Punkt a,b,c;
  Diffuseur *pdiff;
  double distZ,pente;
  Vecteur SvE=Ecran[0]; // SvE : Scene vers Ecran
  bool up,down,pastoutvu;
  //taille ecran stocke en tmp et dz
  //variable pr periodic infini
  register signed char acv_idx,acv_fin=0;
  Vecteur delta[3];
  
  delta[0][0]=delta[2][0]=0.0;
  delta[0][1]=delta[1][1]=0.0;
  delta[0][2]=delta[1][2]=delta[1][2]=0.0;
  delta[1][0]=vmax[0]-vmin[0];
  delta[2][1]=vmax[1]-vmin[1];
  
  // Cas des primitives (non capteurs virtuels)
  // int comptr;
  //comptr=0;
  for(Ldiff.debut();(! Ldiff.finito()) && Ldiff.contenu()->isreal();Ldiff.suivant()){
    //Ferr<<"* diff no. "<<++comptr<<endl;
    pdiff=Ldiff.contenu();
    pastoutvu=false;
    //cout <<"Canopy[projplan] primitive = "<<pdiff->primi().name()<<endl;
    //cout <<"Canopy[projplan] P{Re} = ";Ecran[i+1].show();
    /*if( (Ecran[i+1][2]<0.0) || pastoutvu){  // Prim  PARTIELLEMENT pas vue
      if( Ecran[i+1][2]<0.0) {
      pastoutvu=true; 
      //cout <<"Camera[calc_visi] Z neg\n";
      }
      }
      else{
      */
    //cout <<"Camera[projplan] Pp{Rimage}["<<i<<"]  = ";Pp[i].show();
    switch(pdiff->acv) {
    case 0: acv_fin=0; break;
    case 1: acv_fin=1; break;
    case 2: 
    case 4: acv_fin=2; break;
    }
    if(!pastoutvu) 
      for(acv_idx=0;acv_idx<=acv_fin;acv_idx++) {
	if(acv_idx==1 && pdiff->acv==2)
	  acv_idx++;
	for(i=0;i<3;i++) { // Cas des triangles
	  Ecran[i+1]=pdiff->primi()[i];
	  Ecran[i+1]+=delta[acv_idx];
	  Ecran[i+1]-=SvE;
	  Ecran[i+1]=Ecran[i+1].chgt_base(v,w,u);
	  Pp[i][0]=Ecran[i+1][0];
	  Pp[i][1]=Ecran[i+1][1];
	  Pp[i][2]=Ecran[i+1][2];//distZ;
	}//for triangle
	// Tri sommets tq Pp[i][1]<<Pp[j][1]<<Pp[k][1] ie A[1] < B[1] < C[1]
	j = (Pp[1][1]>Pp[2][1])? 1: 2; // calc intermed
	k = (Pp[0][1]>Pp[j][1])? 0: j; // indice max pour coord y
      i = (k+1)%3; j= (i+1)%3;
      i = (Pp[i][1]<Pp[j][1])? i : j; // indice min pour coord y
      j = 3- i-k;
      //       cout<<" (i,j,k) = "<<i<<j<<k<<endl;
      if ((i!=k)&& !((A[0]==B[0])&&(B[0]==C[0]))&& !((A[1]==B[1])&&(B[1]==C[1]))){
	// Pts A,B,Cpas  alignes selon les axes Xou Y
	// tri Ok
	up=down=false;
	if (A[1]==B[1]) // up 
	{ i = (A[0] <B[0])?i:j;
	j = 3-i-k; 
	up=(A[0]==B[0])?false: true;
	}//if up
	else{
	  if(B[1] == C[1]){ // down 
	    k = (B[0] <C[0])?k:j; 
	    j = 3-i-k;
	    down =(B[0]==C[0])?false: true;
	  }//if down
	  else { 
	    D[1] = B[1];
	    pente=(D[1]-A[1])/(C[1]-A[1]);
	    D[0] = pente*(C[0]-A[0])+A[0];
	    // D[2] = (A[2]*(C[1]-D[1]) + C[2]*(D[1]-A[1]))/(C[1]-A[1]);
	    D[2] = pente*(C[2]-A[2])+A[2];
	    up=down=(B[0]==D[0])?false: true;
	    if(D[0]>B[0]) { l=i; i=j; j=3;}
	    else          { l=i; i=3;     }
	  }//else cas quelconque, ni up , ni down 
	}
	//     cout<<"2 (i,j,k) = "<<i<<j<<k<<endl;
	// rappel syntaxe colorie_triangle()
	//void colorie_triangle( void * tria,void *Zprim, double **Zbuf,Punkt a,Punkt b,Punkt c, int Tx, int Ty,double tx, double ty,void (*f)(void * Zprim, int i, int j, void* tria)){

	if(up) {
	  pt2pkt(A,a);
	  pt2pkt(B,b);
	  pt2pkt(C,c);
	  /* printf(" UP : \tA[0]= %g, A[1]=%g, A[2]=%g\n",A[0], A[1],A[2]);
	  printf(" \t\tB[0]= %g, B[1]=%g, B[2]=%g\n",B[0], B[1],B[2]);
	  printf(" \t\tC[0]= %g, C[1]=%g, C[2]=%g\n",C[0], C[1],C[2]);
	  */
	  colorie_triangle(pdiff,Zprim,Zbuf, c,a,b,Timg,Timg,du,dv,zproj);
	}
	if(down) {
	  if (up) { k=j; j=i; i=l; }
	  pt2pkt(A,a);
	  pt2pkt(B,b);
	  pt2pkt(C,c);
	  /* printf(" DOWN : \tA[0]= %g, A[1]=%g, A[2]=%g\n",A[0], A[1],A[2]);
	  printf(" \t\tB[0]= %g, B[1]=%g, B[2]=%g\n",B[0], B[1],B[2]);
	  printf(" \t\tC[0]= %g, C[1]=%g, C[2]=%g\n",C[0], C[1],C[2]);
	  */
	  colorie_triangle(pdiff,Zprim,Zbuf, a,b,c, Timg,Timg,du,dv,zproj);
	}// if down
      }//if pas un triangle plat
    }//if !pastoutvu 
  }// for liste diffuseurs

  //Infinitisation
  if(infty && visee[2]>-1+1e-6) {
    int **roofi;
    double cdist;//cste de distance ne depend que de l'inclinaison de la visee
    roofi=new int*[4];
    for(i=0;i<4;i++) {
      roofi[i]=new int[2];
      roof[i]-=SvE;
      roof[i]=roof[i].chgt_base(v,w,u);
      roofi[i][0] = (int)(roof[i][0]*(Timg)/du);
      roofi[i][1] = (int)(roof[i][1]*(Timg)/dv);
      }
    cdist=tan(Macos(-visee[2]))*dv/(double)Timg;
    //infinitisation sans duplication des primi (juste ombrage)
    infinitise((void ***)Zprim,Zbuf,cdist,roofi,Timg,Timg,false);
    for(i=0;i<4;i++)
      delete [] roofi[i];
    delete [] roofi;
  }//if infty
  
//Valeur d'un pixel 
  Apix=du*dv/(double)(Timg*Timg)/costeta;

//Ferr <<"==> Traitement des "  << nbcell<<" Capt\n" ;

  // Cas des capteurs virtuels 
  if(nbcell>0){  
    int nbpix;
    for( ;! Ldiff.finito();Ldiff.suivant()){
      pdiff=Ldiff.contenu();
      if(pdiff->isreal()){
	Ferr <<"<!> Attention triangle reel dans la liste capteur virtuel!!\n";
	continue;
      }
      // Ferr<<"\t prodscal="<<visee.prod_scalaire(pdiff->normal())<<endl;
      if(visee.prod_scalaire(pdiff->normal())<0){//capteur bien vu par dessus
	nbpix=0;
	switch(pdiff->acv) {
	case 0: acv_fin=0; break;
	case 1: acv_fin=1; break;
	case 2: 
	case 4: acv_fin=2; break;
	}
	for(acv_idx=0;acv_idx<=acv_fin;acv_idx++) {
	  if(acv_idx==1 && pdiff->acv==2)
	    acv_idx++;
	  for(i=0;i<3;i++) { // Cas des triangles
	    Ecran[i+1]=pdiff->primi()[i];
	    Ecran[i+1]+=delta[acv_idx];
	    Ecran[i+1]-=SvE;
	    Ecran[i+1]=Ecran[i+1].chgt_base(v,w,u);
	    Pp[i][0]=Ecran[i+1][0];
	    Pp[i][1]=Ecran[i+1][1];
	    Pp[i][2]=Ecran[i+1][2];//distZ;
	  }//for triangle
	  // Tri sommets tq Pp[i][1]<<Pp[j][1]<<Pp[k][1] ie A[1] < B[1] < C[1]
	  j = (Pp[1][1]>Pp[2][1])? 1: 2; // calc intermed
	  k = (Pp[0][1]>Pp[j][1])? 0: j; // indice max pour coord y
	  i = (k+1)%3; j= (i+1)%3;
	  i = (Pp[i][1]<Pp[j][1])? i : j; // indice min pour coord y
	  j = 3- i-k;
	  if ((i!=k)&& !((A[0]==B[0])&&(B[0]==C[0]))&& !((A[1]==B[1])&&(B[1]==C[1]))){
	    // Pts A,B,Cpas  alignes selon les axes Xou Y
	    up=down=false;
	    if (A[1]==B[1]) // up 
	      { i = (A[0] <B[0])?i:j;
	      j = 3-i-k; 
	      up=(A[0]==B[0])?false: true;
	      }//if up
	    else{
	      if(B[1] == C[1]){ // down 
		k = (B[0] <C[0])?k:j; 
		j = 3-i-k;
		down =(B[0]==C[0])?false: true;
	      }//if down
	      else { 
		D[1] = B[1];
		pente=(D[1]-A[1])/(C[1]-A[1]);
		D[0] = pente*(C[0]-A[0])+A[0];
		D[2] = pente*(C[2]-A[2])+A[2];
		up=down=(B[0]==D[0])?false: true;
		if(D[0]>B[0]) { l=i; i=j; j=3;}
		else          { l=i; i=3;     }
	      }//else cas quelconque, ni up , ni down 
	    }
	    // rappel syntaxe colorie_capteur()
	    //void colorie_capteur(  double **Zbuf,Punkt a,Punkt b,Punkt c, int Tx, int Ty,double tx, double ty,double nbpix)){
	    
	    if(up) {
	      pt2pkt(A,a);
	      pt2pkt(B,b);
	      pt2pkt(C,c);
	      colorie_capteur(Zbuf, c,a,b,Timg,Timg,du,dv,nbpix);
	    }
	    if(down) {
	      if (up) { k=j; j=i; i=l; }
	      pt2pkt(A,a);
	      pt2pkt(B,b);
	      pt2pkt(C,c);
	      colorie_capteur(Zbuf, a,b,c, Timg,Timg,du,dv,nbpix) ;
	    }// if down
	  }//if pas un triangle plat
	  //maj de l'eclairement 
	  //Ferr <<"Capt no. "  << (int)(pdiff->num())<<" => nbpix=" 
	  //     << nbpix<<"\n" ;

	  Bo[pdiff->num()]+=nbpix*Apix;
	}//for acv
      }//if vu par au-dessus
    }// for liste diffuseurs
  }  // Fin traitt Capteurs virtuels 
  
  //maj de l'image en fonction de Zprim
  double cocnomen,alt;
  int nummer;
  double cocmax=0, cocmin=99999999999.,altmin=99999999.,altmax=0;
  //calcul de l'eclairage direct
  if(verbose>1) printf("projplan() : du=%lf - dv=%lf =>  Apix= %lf\n",du,dv,Apix);
  for(i=0;i<Timg;i++)
    for(j=0;j<Timg;j++) {
      pdiff=Zprim[i][j];
      if(0){
	cocnomen=(pdiff==NULL)? 0:pdiff->primi().name()/1e7;
	alt=(pdiff==NULL)? 0: 100*Zbuf[i][j];
	cocmax=(cocnomen>cocmax) ? cocnomen:cocmax;
	cocmin=(cocnomen<cocmin) ? cocnomen:cocmin;
	altmax=(alt>altmax) ? alt:altmax;
	altmin=(alt<altmin) ? alt:altmin;
      }
      //calcul de la visibilite'
      if(pdiff!=NULL) {
	nummer=pdiff->num();
	pdiff->active(visee);
        //printf("projplan : img(%d,%d)=%d\n",i,j,pdiff->num());
	Bo[pdiff->num()]+=Apix;
	if(!pdiff->isopaque()) {
	  pdiff->togle_face();
	  Bo[pdiff->num()]-=Apix;
	  pdiff->togle_face();
	}//si transparent
	pdiff->activ_num(nummer);
        if(pdiff->num()!=nummer){
	  Ferr <<" pdiff->num()!=nummer\n" ;
          exit(5);
        }
      }
    }

 //ecriture de projplan.ppm et de prozplan.ppm
 ////////////////////// GAFFE !! ecriture par fprintf ET fwrite 
 //                              dans le même fichier
 //////////////////////          Ne fonctionnera pas sous Win
 /*FILE*fz,*fprim;
 unsigned char bit;
 fz=fopen("proz.ppm","w");
 fprim=fopen("proj.ppm","w");
 fprintf(fz,"P5\n# Radioxity - Image des z \n# BitMin=%lf - BitMax=%lf\n %d %d\n255\n",altmin,altmax,Timg,Timg);
 fprintf(fprim,"P5\n# Radioxity - Image des labels \n# BitMin=%lf - BitMax=%lf\n %d %d\n255\n",cocmin,cocmax,Timg,Timg);
 altmax+=(altmax-altmin)/10.;
 cocmin-=((cocmax-cocmin>0?cocmax-cocmin:1))/10.;
 for(i=Timg-1;i>=0;i--)
   for(j=Timg-1;j>=0;j--) {
     pdiff=Zprim[i][j];
     alt=(pdiff==NULL)? altmax: 100*Zbuf[i][j];
     bit= (unsigned char) ((alt-altmin)/(altmax-altmin)*255);
     bit=(bit>255|| bit<0)?255:bit;
     fwrite(&bit,sizeof(char),1,fz);
     cocnomen=(pdiff==NULL)? cocmin:pdiff->primi().name()/1e10;
     bit= (unsigned char) ((cocnomen-cocmin)/(cocmax-cocmin)*255);
     bit=(bit>255||bit<0)?255:bit;
     fwrite(&bit,sizeof(unsigned char),1,fprim);
   }
 fclose(fz);
 fclose(fprim);
 */
 // libere les Xbuff
 pZbuf=Zbuf;
 for(i=0;i<Timg;i++,pZbuf++) 
   delete [] (*pZbuf);
 delete []  Zbuf;
 for(i=0,pZprim=Zprim;i<Timg;i++,pZprim++) 
   delete [] (*pZprim);
 delete Zprim;
 if(verbose>2) printf("<= projplan() FIN\n%c",7);
}//Canopy::projplan()



//-*************** colorie_triangle() ************************

void colorie_triangle( void * tria,void *Zprim, REELLE **Zbuf,const Punkt a,const Punkt b,const Punkt c, int Tx, int Ty,double tx, double ty,void (*f)(void *, int, int, void*)){
  double penteL,penteR,xL,xR,zL,yrel,vab,vac,vbc,z,dz; // L Left, R  Right
  register int i,j, iR,iL, deby,finy;
  register signed char sens;
  double K[2];
  
  //            a           b ______ c 
  //           /\             \    /
  //          /  \      ou     \  /
  //        b/____\c            \/ a
  /*                                   */
  K[0]=(Tx-1)/tx;
  K[1]=(Ty-1)/ty;
  sens   = (b[1] > a[1])? 1 : -1;         
  penteL = (b[0] - a[0]) / (b[1] - a[1]);   
  penteR = (c[0] - a[0]) / (c[1] - a[1]);  
  
  
  if(sens>0){
    deby   = int(a[1]*K[1]+ 0.5);
    finy   = int(b[1]*K[1]- 0.5);
  }
  else{
    deby   = int(b[1]*K[1]+ 0.5);
    finy   = int(a[1]*K[1]- 0.5);
  }

  //printf("Canopy::colorie_triangle(): deby=%d, finy=%d \n",deby,finy); 
  
  if( finy<0 || deby>=Ty || deby>finy )
    return;  
  
  deby=max(0,deby);
  finy=min(Ty-1,finy);
  
  //	calcul de grandeurs fixes utilisees pour obtenir l'altitude des pixels
  vab = (b[2] - a[2])/(b[1] - a[1]);
  vbc = (c[2] - b[2])/(b[1] - a[1]);
  
  //	calcul de grandeurs qui seront incrementes dans la boucle sur les lignes
  yrel=((((double)deby+0.5)/K[1])-a[1]);
  xL=penteL*yrel+a[0];
  xR=penteR*yrel+a[0];
  zL= a[2] + yrel*vab + ((c[2]-b[2])/(c[0]-b[0])*0.5/K[0]);
  dz=vbc/(penteR-penteL)/K[0];
  penteL/=K[1];
  penteR/=K[1];
  vab/=K[1]; 
  //K[1]=-sens/K[1];
    //cout<<" Canopy[colorie_triangle] \n ";
  //cout<<"\t A : ";a.show();cout<<"\t B : "; b.show();cout<<"\t C : "; c.show();
  //cout<<"\t deby = "<<deby  <<" - finy = "<<finy<<" - yrel = "<<yrel<<endl;
  //cout<<"\t xL   = "<<xL    <<" - xR   = "<<xR<<" -zL="<<zL<<endl;
  //cout<<"\tpenteL= "<<penteL<<" -penteR= "<<penteR<<" - sens = "<<(int) sens<<" - vab= "<<vab<<" - vbc = "<<vbc<<endl;
  
  //	boucle sur les lignes 
  for(j=deby;j<=finy;j++){
    //iL=min(Tx-1,max(0,int(xL*K[0]+0.5)));    // centre de gravite 
    //iR=min(Tx-1,max(0,int(xR*K[0]-0.5)));    // dans le triangle
    iL=max(0,int(xL*K[0]+0.5));
    iR=min(Tx-1,int(xR*K[0]-0.5)); 
    //cout<<"LOOP\t xL= "<<xL    <<" - xR= "<<xR<< "(iL,iR)="<<iL<<", "<<iR<<" -j= "<<j<<" -zL="<<zL<< endl;
    if (iL<=iR) {
      //dz= yrel*vbc/(xR-xL);
      z = zL;
      for(i=iL;i<=iR;i++){
	//cout<<"Capteur[colorie-triangle] [i][j] = "<<i<<", "<<j<<endl;
	//cout <<"Camera [colorie_triangle][i][j] ("<<i<<","<<j<<") z = "<<z<<"- Z = "<<Zbuf[i][j]<<endl;
	if(z<Zbuf[i][j]){
	      //cout <<"Camera [colorie_triangle][i][j] ("<<i<<","<<j<<") z = "<<z<<"- Z = "<<Zbuf[i][j]<<endl;
	  Zbuf[i][j]=(REELLE) z;
	  //ndi[i][j]=idif; 
	  f(Zprim,i,j, tria);
	  //if(pimg!=NULL)pimg->maj(Tx-1-i,Ty-1-j,pdif->primi().name());
	} 
	z = z+dz;
      }//for i
    }//if ligne pleine
    //yrel+=K[1];
    xL+=penteL;  
    xR+=penteR;
    zL+=vab;
  }// for j	
  //   cout<<"FIN\t i = "<<i  <<" - j = "<<j<<" - yrel = "<<yrel<<endl;
  //   cout<<"\t xL   = "<<xL    <<" - xR   = "<<xR<<endl;
}// colorie_triangle()


/***** colorie_capteur() *****/
 void colorie_capteur(REELLE **Zbuf,Punkt a,Punkt b,Punkt c, int Tx, int Ty,double tx, double ty,int& pB0){
 double penteL,penteR,xL,xR,zL,yrel,vab,vac,vbc,z,dz; // L Left, R  Right
  register int i,j, iR,iL, deby,finy;
  register signed char sens;
  double K[2];
  
  //            a           b ______ c 
  //           /\             \    /
  //          /  \      ou     \  /
  //        b/____\c            \/ a
  /*                                   */
  K[0]=(Tx-1)/tx;
  K[1]=(Ty-1)/ty;
  sens   = (b[1] > a[1])? 1 : -1;         
  penteL = (b[0] - a[0]) / (b[1] - a[1]);   
  penteR = (c[0] - a[0]) / (c[1] - a[1]);  
  
  
  if(sens>0){
    deby   = int(a[1]*K[1]+ 0.5);
    finy   = int(b[1]*K[1]- 0.5);
  }else{
    deby   = int(b[1]*K[1]+ 0.5);
    finy   = int(a[1]*K[1]- 0.5);
  }
  if( finy<0 || deby>=Ty || deby>finy )
    return;  
  
  deby=max(0,deby);
  finy=min(Ty-1,finy);
  //	calcul de grandeurs fixes utilisees pour obtenir l'altitude des pixels
  vab = (b[2] - a[2])/(b[1] - a[1]);
  vbc = (c[2] - b[2])/(b[1] - a[1]);
  //	calcul de grandeurs qui seront incrementes dans la boucle sur les lignes
  yrel=((((double)deby+0.5)/K[1])-a[1]);
  xL=penteL*yrel+a[0];
  xR=penteR*yrel+a[0];
  zL= a[2] + yrel*vab + ((c[2]-b[2])/(c[0]-b[0])*0.5/K[0]);
  dz=vbc/(penteR-penteL)/K[0];
  penteL/=K[1];
  penteR/=K[1];
  vab/=K[1]; 
  
  //	boucle sur les lignes 
  for(j=deby;j<=finy;j++){
    iL=max(0,int(xL*K[0]+0.5));
    iR=min(Tx-1,int(xR*K[0]-0.5)); 
    if (iL<=iR) {
      z = zL;
      for(i=iL;i<=iR;i++){
	if(z<Zbuf[i][j]){
	  // printf("z=%.3f<%.3f=Z(%d,%d)\n",z,Zbuf[i][j],i,j);
	  pB0++;
	} 
	//else printf("z=%.3f>%.3f=Z(%d,%d)\n",z,Zbuf[i][j],i,j);
	z = z+dz;
      }//for i
    }//if ligne pleine
    xL+=penteL;  
    xR+=penteR;
    zL+=vab;
  }// for j	

}//colorie_capteur()
 




/*******************************************
**************   data3d     **************** 
********************************************/

void Canopy::data3d(int tx,int ty,Vecteur &visee,
		    bool infty,long int **&Zno){
  register int i,j,k,l,img_surf;
  Image imgZ (tx,ty,(char*)"proz.ppm");
  long int **pZno;
  REELLE **Zbuf,**pZbuf,*ptZ;
  float emax[3],emin[3],dz;
  Point roof[4];
  Diffuseur ***Zprim,***pZprim;
  //   Tabdyn<double, 2> Zbuf(img->taille(0),img->taille(1));
  //Tabdyn<Diffuseur *, 2> Zprim(img->taille(0),img->taille(1));
  l=0;
  Zno= new long*[imgZ.taille(0)];
  pZno=Zno;
  for(i=0;i<imgZ.taille(0);i++,pZno++) {
    (*pZno)=new long[imgZ.taille(1)];
    for(j=0;j<imgZ.taille(1);j++)
      (*pZno)[j]=-1;
  }
  
  Zbuf= new REELLE*[imgZ.taille(0)];
  pZbuf=Zbuf;
  img_surf=imgZ.taille(0)*imgZ.taille(1);
  for(i=0;i<imgZ.taille(0);i++,pZbuf++) {
    (*pZbuf)=new REELLE[imgZ.taille(1)];
  }
  pZbuf=Zbuf;
  for(i=0;i<imgZ.taille(0);i++,pZbuf++) {
    ptZ=*pZbuf;
    for(j=0;j<imgZ.taille(1);j++,ptZ++)
      *ptZ=99999999999.9;
  }
  Zprim=new Diffuseur**[imgZ.taille(0)];
  for(i=0,pZprim=Zprim;i<imgZ.taille(0);i++,pZprim++) {
    (*pZprim)=new Diffuseur*[imgZ.taille(1)];
    for(j=0;j<imgZ.taille(1);j++,ptZ++)
    (*pZprim)[j]=NULL;
  }

  //&&&&&& Data3d &&&&&&&&&
  //calcul de la position de l'ecran en fonction des bornes de la scene
  Point Ecran[4];
  double Pts[8][2],du,dv;
  int ofset,ei[4];
  Vecteur u,v,w;
 
  for (i=0;i<2;i++){
    emin[i]=vmin[i];
    emax[i]=vmax[i];    
  }
  visee.normalise();
  cout<<"Visee::";visee.show();
  if(visee[2]<=0){//visee par  dessus
    printf("=> Down: Visee par dessus...\n");
    ei[0]=0; ei[1]=1; ei[2]=2; ei[3]=3; 
    emin[2]=vmin[2];
    emax[2]=vmax[2];    
    ofset=1;
  }else{
    printf("=> Up: Visee par dessous...\n");
    ei[0]=3; ei[1]=2; ei[2]=1; ei[3]=0; 
    emin[2]=vmax[2];
    emax[2]=vmin[2];    
    ofset=-1;
  }
  
   if(fabs(visee[2])+1e-5>1) {// visee verticale
     Ecran[0][0]=emin[0];
     Ecran[0][1]=emax[1];
     Ecran[1][0]=emin[0];
     Ecran[1][1]=emin[1];
     Ecran[2][0]=emax[0];
     Ecran[2][1]=emin[1];
     Ecran[3][0]=emax[0];
     Ecran[3][1]=emax[1];
     Ecran[0][2]=Ecran[1][2]=Ecran[2][2]=Ecran[3][2]=emax[2]+ofset;
    du=emax[0]-emin[0];
    dv=emax[1]-emin[1];
    u=visee;
    v[0]=1; v[1]=0; v[2]=0;
    w[0]=0; w[1]=-1; w[2]=0;
    if(ofset==-1) w[1]=1;
    cout<<"Base (u,v,w) :\n";u.show(); v.show(); w.show();
  }
  else { //cas  : visee non verticale
    double M[2][2],tmp,maxi[2]={-999999.9,-999999.9},mini[2]={999999.9,999999.9};
    dz=emax[2]-emin[2];
    printf("==> dz=%.3g\n",dz);
    //points a Z=0   
    Ecran[0][0]=Pts[0][0]=emin[0];
    Ecran[0][1]=Pts[0][1]=emin[1];
    Ecran[1][0]=Pts[1][0]=emax[0];
    Ecran[1][1]=Pts[1][1]=emin[1];
    Ecran[2][0]=Pts[2][0]=emax[0];
    Ecran[2][1]=Pts[2][1]=emax[1];
    Ecran[3][0]=Pts[3][0]=emin[0];
    Ecran[3][1]=Pts[3][1]=emax[1];
    //projection selon visee des points a z = emax[2]
    for(i=0;i<4;i++) {
      roof[i]=Ecran[i];
      Ecran[i][2]=emax[2];
      Ecran[i]=Ecran[i]+visee*fabs(dz/visee[2]);
      Pts[i+4][0]=Ecran[i][0];
      Pts[i+4][1]=Ecran[i][1];
      printf("data3d() : roof[%d] = (%f,%f,%f)\n",i,roof[i][0],roof[i][1],roof[i][2] ) ; 
      Ecran[i][2]=emin[2];
    }
    roof[0][2]=roof[1][2]=roof[2][2]=roof[3][2]=emax[2];
    //chgt de repere x,y -> u,v && calculer min et max
     u=visee;
     u[2]=0;
     u.normalise();
     v[0]=-u[1];
     v[1]=u[0];
     v[2]=0;
     w=visee.prod_vectoriel(v); 
     w.normalise();
     if(w[2]<0) {
       w=-w;
       v=-v;
     }
     //if(ofset==-1) w=-w;
     cout<<"Base (u,v,w) :\n";u.show(); v.show(); w.show();

    M[0][0]=u[0] ;
    M[1][0]=v[0]; 
    M[0][1]=u[1] ;
    M[1][1]=v[1] ;
    for(i=0;i<8;i++) {
      tmp=Pts[i][0];
      for(j=0;j<2;j++) {
	Pts[i][j]=M[j][0]*tmp+M[j][1]*Pts[i][1];
	mini[j]=min(mini[j],Pts[i][j]);
	maxi[j]=max(maxi[j],Pts[i][j]);      
      }
    }
    //coord de la projection de l'ecran sur le sol (Ruv)
    Ecran[0][0]=mini[0];
    Ecran[0][1]=mini[1];
    Ecran[1][0]=maxi[0];
    Ecran[1][1]=mini[1];
    Ecran[2][0]=maxi[0];
    Ecran[2][1]=maxi[1];
    Ecran[3][0]=mini[0];
    Ecran[3][1]=maxi[1];
    //chgt repere uv -> xy
    M[0][0]=u[0] ;
    M[1][0]=u[1]; ;
    M[0][1]=v[0] ;
    M[1][1]=v[1] ;
    for(i=0;i<4;i++) {
      tmp=Ecran[i][0];
      //Ecran[i].show();
      for(j=0;j<2;j++) {
	Ecran[i][j]=M[j][0]*tmp+M[j][1]*Ecran[i][1];
      }
      //Ecran[i].show();
    }
    //redressement de l'ecran tq ortho a visee (Ecran 0 et 3 : u min)
    dv=tmp=(maxi[0]-mini[0])*u.prod_scalaire(w);
    //printf("mini[0] = %g - mini[1] = %g ; maxi[0]= %g - maxi[1] = %g\n",mini[0],mini[1],maxi[0], maxi[1]);
    //printf(" translation de d de la base Ecran : d = %g\n",tmp);
    Ecran[1]=Ecran[0]+w*tmp;
    Ecran[2]=Ecran[3]+w*tmp;
    // cout<<"dv="<<dv<<" - dv2 = "<<Ecran[0].dist(Ecran[1])<<endl;
    //translation de l'ecran tq tota la scene soit vue
    tmp=(dz+ofset)/visee[2];
    printf("## ofset=%d\n",ofset);
    for(i=0;i<4;i++){
     	Ecran[i]= Ecran[i] + visee*tmp;     
    }
     printf("##\n");
    du = maxi[1]-mini[1];
    u=visee;
  }//else visee verticale
   // if(ofset==-1){
//      Ecran[4]=Ecran[3];
//      Ecran[3]= Ecran[0];
//      Ecran[0]= Ecran[4];
//      Ecran[4]=Ecran[2];
//      Ecran[2]= Ecran[1];
//      Ecran[1]= Ecran[4];
//    }
 //validation geom
  FILE *fcan;
  fcan=fopen("proj.can","w");
  // Bounding Box
  fprintf(fcan,"p 1 1 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",vmin[0],vmin[1],vmin[2], vmax[0],vmin[1],vmin[2], vmax[0],vmin[1],vmax[2], vmin[0],vmin[1],vmax[2]);
  fprintf(fcan,"p 1 2 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",vmax[0],vmin[1],vmin[2], vmax[0],vmax[1],vmin[2], vmax[0],vmax[1],vmax[2], vmax[0],vmin[1],vmax[2]);
  fprintf(fcan,"p 1 3 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",vmax[0],vmax[1],vmin[2], vmin[0],vmax[1],vmin[2], vmin[0],vmax[1],vmax[2], vmax[0],vmax[1],vmax[2]);
  fprintf(fcan,"p 1 4 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",vmin[0],vmax[1],vmin[2], vmin[0],vmin[1],vmin[2], vmin[0],vmin[1],vmax[2], vmin[0],vmax[1],vmax[2]);
  fprintf(fcan,"p 1 5 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",vmin[0],vmin[1],vmax[2], vmax[0],vmin[1],vmax[2], vmax[0],vmax[1],vmax[2], vmin[0],vmax[1],vmax[2]);
  //Ecran de projection
  fprintf(fcan,"p 1 0 4 %g %g %g  %g %g %g  %g %g %g  %g %g %g\n",Ecran[0][0],Ecran[0][1],Ecran[0][2], Ecran[1][0],Ecran[1][1],Ecran[1][2], Ecran[2][0],Ecran[2][1],Ecran[2][2], Ecran[3][0],Ecran[3][1],Ecran[3][2]);
  fclose(fcan);


  //Projection parallele a visee sur Ecran
  Point Pp[4];
  Punkt a,b,c;
  Diffuseur *pdiff;
  double distZ,pente;
  Vecteur SvE; // SvE : Scene vers Ecran
  Vecteur delta[3];
  bool up,down,pastoutvu;
  register signed char acv_idx,acv_fin=0;
  //taille ecran stocke en tmp et dz

  imgZ.raz(-10);
  delta[0][0]=delta[2][0]=0.0;
  delta[0][1]=delta[1][1]=0.0;
  delta[0][2]=delta[1][2]=delta[1][2]=0.0;
  delta[1][0]=emax[0]-emin[0];
  delta[2][1]=emax[1]-emin[1];

  if(ofset==1) SvE=Ecran[0]; else SvE=Ecran[1];
  for(Ldiff.debut();! Ldiff.finito();Ldiff.suivant()){
    pdiff=Ldiff.contenu();
    pastoutvu=false;
    //cout <<"Canopy[data3D] primitive = "<<pdiff->primi().name()<<endl;
    //cout <<"Canopy[data3d] P{Re} = ";Ecran[i+1].show();
      /*if( (Ecran[i+1][2]<0.0) || pastoutvu){  // Prim  PARTIELLEMENT pas vue
	if( Ecran[i+1][2]<0.0) {
	pastoutvu=true; 
	//cout <<"Camera[calc_visi] Z neg\n";
	}
	}
	else{
	*/
    //cout <<"Camera[data3d] Pp{Rimage}["<<i<<"]  = ";Pp[i].show();
    //}//else P[i][2]<0.0) || pastoutvu)
    switch(pdiff->acv) {
    case 0: acv_fin=0; break;
    case 1: acv_fin=1; break;
    case 2: 
    case 4: acv_fin=2; break;
    }
    if(!pastoutvu)
      for(acv_idx=0;acv_idx<=acv_fin;acv_idx++) {
	if(acv_idx==1 && pdiff->acv==2)
	  acv_idx++;
	for(i=0;i<3;i++) { // Cas des triangles
	  Ecran[i+1]=pdiff->primi()[i];
	  Ecran[i+1]+=delta[acv_idx];
	  Ecran[i+1]-=SvE;
	  Ecran[i+1]=Ecran[i+1].chgt_base(v,w,u);
	  Pp[i][0]=Ecran[i+1][0];
	  Pp[i][1]=Ecran[i+1][1];
	  Pp[i][2]=Ecran[i+1][2];//distZ;
	  //cout <<"Camera[data3d] Pp{Rimage}["<<i<<"]  = ";Pp[i].show();
	}//for triangle 
	
	// Tri sommets tq Pp[i][1]<<Pp[j][1]<<Pp[k][1] ie A[1] < B[1] < C[1]
	j = (Pp[1][1]>Pp[2][1])? 1: 2; // calc intermed
	k = (Pp[0][1]>Pp[j][1])? 0: j; // indice max pour coord y
	i = (k+1)%3; j= (i+1)%3;
	i = (Pp[i][1]<Pp[j][1])? i : j; // indice min pour coord y
	j = 3- i-k;
	//       cout<<" (i,j,k) = "<<i<<j<<k<<endl;
	if ((i!=k)&& !((A[0]==B[0])&&(B[0]==C[0]))&& !((A[1]==B[1])&&(B[1]==C[1]))){
	  // Pts A,B,Cpas  alignes selon les axes Xou Y
	  // tri Ok
	  up=down=false;
	  if (A[1]==B[1]) // up 
	  { i = (A[0] <B[0])?i:j;
	  j = 3-i-k; 
	  up=(A[0]==B[0])?false: true;
	  }//if up
	  else{
	    if(B[1] == C[1]){ // down 
	      k = (B[0] <C[0])?k:j; 
	      j = 3-i-k;
	      down =(B[0]==C[0])?false: true;
	    }//if down
	    else { 
	      D[1] = B[1];
	      pente=(D[1]-A[1])/(C[1]-A[1]);
	      D[0] = pente*(C[0]-A[0])+A[0];
	      // D[2] = (A[2]*(C[1]-D[1]) + C[2]*(D[1]-A[1]))/(C[1]-A[1]);
	      D[2] = pente*(C[2]-A[2])+A[2];
	      up=down=(B[0]==D[0])?false: true;
	      if(D[0]>B[0]) { l=i; i=j; j=3;}
	      else          { l=i; i=3;     }
	    }//else cas quelconque, ni up , ni down 
	  }
	  //     cout<<"2 (i,j,k) = "<<i<<j<<k<<endl;
	  // rappel syntaxe colorie_triangle()
	  //void colorie_capteur( void * tria,void *Zprim, double **Zbuf,Punkt a,Punkt b,Punkt c, int Tx, int Ty,double tx, double ty,void (*f)(void * Zprim, int i, int j, void* tria)){
	  
	  if(up) {
	    pt2pkt(A,a);
	    pt2pkt(B,b);
	    pt2pkt(C,c);
	    /* printf(" UP : \tA[0]= %g, A[1]=%g, A[2]=%g\n",A[0], A[1],A[2]);
	       printf(" \t\tB[0]= %g, B[1]=%g, B[2]=%g\n",B[0], B[1],B[2]);
	       printf(" \t\tC[0]= %g, C[1]=%g, C[2]=%g\n",C[0], C[1],C[2]);
	       */  
	    colorie_triangle(pdiff,Zprim,Zbuf, c,a,b,imgZ.taille(0),imgZ.taille(1),du,dv,zproj);
	  }
	  if(down) {
	    if (up) { k=j; j=i; i=l; }
	    pt2pkt(A,a);
	    pt2pkt(B,b);
	    pt2pkt(C,c);
	    /* printf(" DOWN : \tA[0]= %g, A[1]=%g, A[2]=%g\n",A[0], A[1],A[2]);
	       printf(" \t\tB[0]= %g, B[1]=%g, B[2]=%g\n",B[0], B[1],B[2]);
	       printf(" \t\tC[0]= %g, C[1]=%g, C[2]=%g\n",C[0], C[1],C[2]);
	       */ 
	    colorie_triangle(pdiff,Zprim,Zbuf, a,b,c, imgZ.taille(0),imgZ.taille(1),du,dv,zproj);
	  }// if down
	}//if pas un triangle plat
      }//if !pastoutvu 
  }// for liste diffuseurs
  //maj de l'image en fonction de Zprim
  double cocnomen,alt;
  //Infinitisation
  if(infty && visee[2]>-1+1e-6) {
    int **roofi;
    
    double cdist;//cste de distance ne depend que de l'inclinaison de la visee
    roofi=new int*[4];
    for(i=0;i<4;i++) {
      roofi[i]=new int[2];
      roof[i]-=SvE;
      roof[i]=roof[i].chgt_base(v,w,u);
      roofi[i][0] = (int)(roof[i][0]*(imgZ.taille(0))/du);
      roofi[i][1] = (int)(roof[i][1]*(imgZ.taille(1))/dv);
      }
    /*
      for(i=0;i<4;i++) {
      for(j=0;j<2;j++) {
	printf(" roof(%d,%d) = %f, roofi(%d,%d) = %d\n, ",i,j,roof[i][j],i,j,roofi[i][j]);
	printf("----------------------------------------\n");
      }
    }
    printf(" Before inifinitoise\n");
    */
    cdist=tan(Macos(-visee[2]))*dv/(double)imgZ.taille(1);
    //cout<<" cdist = " <<cdist;
    infinitise((void ***)Zprim,Zbuf,cdist,roofi,imgZ.taille(0),imgZ.taille(1),true);
    for(i=0;i<4;i++)
      delete [] roofi[i];
    delete [] roofi;
  }//if infty
  //remplissage du No buffer
  for(i=0;i<imgZ.taille(0);i++)
    for(j=0;j<imgZ.taille(1);j++) {
      pdiff=Zprim[i][j];
      alt=(pdiff==NULL)? 0: 100*Zbuf[i][j];
      imgZ.maj(imgZ.taille(0)-1-i,imgZ.taille(1)-1-j,alt);
      if(pdiff!=NULL) {
	pdiff->active(visee);
	Zno[i][j]=pdiff->num();
	pdiff->active(0);
      }
    }
  imgZ.sauve();
  //liberez la memoire!
  imgZ.free();
  // libere les Xbuff
  pZbuf=Zbuf;
  for(i=0;i<imgZ.taille(0);i++,pZbuf++) 
    delete [] (*pZbuf);
  delete []  Zbuf;
  for(i=0,pZprim=Zprim;i<imgZ.taille(0);i++,pZprim++) 
    delete [] (*pZprim);
  delete Zprim;
}//Canopy::data3d()






 /******* Obsolete et bugue???? *************/ 
void color_triangle( void * tria,void *Zprim, double **Zbuf,const Punkt a,const Punkt b,const Punkt c, int Tx, int Ty,double tx, double ty,void (*f)(void *, int, int, void*)){
  double penteL,penteR,xL,xR,zL,yrel,vab,vac,vbc,z,dz; // L Left, R  Right
  int i,j, iR,iL, deby,finy;
  int Ae[2],Be[2],Ce[2];
  register signed char sens;
  double K[2];
  
  //            a           b ______ c 
  //           /\             \    /
  //          /  \      ou     \  /
  //        b/____\c            \/ a
  /*                                    */
  K[0]=(Tx-1)/tx;
  K[1]=(Ty-1)/ty;
  Ae[0]=(int)(a[0]*K[0]);
  Ae[1]=(int)(a[1]*K[1]);
  Be[0]=(int)(b[0]*K[0]);
  Be[1]=(int)(b[1]*K[1]);
  Ce[0]=(int)(c[0]*K[0]);
  Ce[1]=(int)(c[1]*K[1]);
   
  sens   = (Be[1] > Ae[1])? 1 : -1;         
  penteL = (Be[0] - Ae[0]) / (double)(Be[1] - Ae[1]);   
  penteR = (Ce[0] - Ae[0]) / (double) (Ce[1] - Ae[1]);  
    
  if(sens>0){
    deby   = int(Ae[1]+ 0.5);
    finy   = int(Be[1]- 0.5);
  }
  else{
    deby   = int(Be[1]+ 0.5);
    finy   = int(Ae[1]- 0.5);
  }
 
  if( finy<0 || deby>=Ty || deby>finy )
    return;  
  
  deby=max(0,deby);
  finy=min(Ty-1,finy);
  
  //	calcul de grandeurs fixes utilisees pour obtenir l'altitude des pixels
  vab = (b[2] - a[2])/(Be[1] - Ae[1]);
  //vbc = (c[2] - Be[2])/(Be[1] - Ae[1]);
  vbc = (c[2] - b[2])/(Be[1] - Ae[1]);
  
  //	calcul de grandeurs qui seront incrementes dans la boucle sur les lignes
  yrel=(double)deby-0.5-Ae[1];
  xL=penteL*yrel+Ae[0];
  xR=penteR*yrel+Ae[0];
  zL= a[2] + yrel*vab ;//+ ((Ce[2]-Be[2])/(Ce[0]-Be[0])*0.5/K[0]);
  dz=vbc/(penteR-penteL);
  //	boucle sur les lignes 
  for(j=deby;j<=finy;j++){
    iL=max(0,int(xL+0.5));    // centre de gravite 
    iR=min(Tx-1,int(xR-0.5));    // dans le triangle
    if (iL<=iR) {
      z = zL;
      for(i=iL;i<=iR;i++){
	if(z<Zbuf[i][j]){
	  Zbuf[i][j]=z;
	  f(Zprim,i,j, tria);
	} 
	z = z+dz;
      }//for i
    }//if ligne pleine
    xL+=penteL;  
    xR+=penteR;
    zL+=vab;
  }// for j	
}// color_triangle()

#ifdef _ZirBouik
void Canopy::calc_FF_Bfar(SPMAT *FF,VEC **Cfar,char * envname, double denv,int nbsim) {
  int Tx=100,Ty=100
  Tabdyn<double, 2> Zbuf(Tx,Ty);
  Tabdyn<Diffuseur *, 2> Zprim(Tx,Ty);
  Punkt P[3],Pp[4];
  Liste<Pointt> lpt; Punkt Pt1,Pt2; ; Primitive *pprim;
  Param_Inter parag;
  Diffuseur *pdiff;
  register int i,j,k,l;
  double distZ,pente;
  Vecteur &w=prim->normal();
  //cout<<" w=vers? ";w.show();
  Vecteur SvE=E,EvI; // SvE : Scene vers Ecran, EvI : Ecran vers Image
  bool up,down,pastoutvu;
  
  Image imgvisi(Tx,Ty, "visi.ppm");
  imgvisi.raz();
  /* Produit vecteur -matrice :pas implemente cf IF 94
     Matrice<double> SvE;
     
     for(i=0;<3;i++){
     SvE(i,0)=u[i];
     SvE(i,1)=v[i];
     SvE(i,2)=w[i];
     }
     */
  Zbuf.maj(9.9E20);
  Zprim.maj(NULL);
  EvI.formation_vecteur(E,(*prim)[1]);
  EvI=EvI.chgt_base(u,v,w);
  //cout<<"Camera[calc_visi]EvI doit z=0";EvI.show();
  for(Ldiff.debut();! Ldiff.finito();Ldiff.suivant()){
    pdiff=Ldiff.contenu();
    /*     if(pdiff->primi().name()==0) //maj visi sol
	   { (*pdiff)++; (*pdiff)++; (*pdiff)--;
	   break;
	   }      
	   */   if(pdiff->primi().name()==0) pdiff->primi().show("Camera[calc_visi] ");
	   pastoutvu=false;
	   //cout <<"Camera[calc_visi] primitive = "<<pdiff->primi().name()<<endl;
	   for(i=0;i<3;i++){ // Cas des triangles
	     P[i]=pdiff->primi()[i];
	     P[i]-=SvE;
	     P[i]=P[i].chgt_base(u,v,w);
	     cout <<"Camera[calc_visi] P{Re} = ";P[i].show();
	     if( (P[i][2]<0.0) || pastoutvu)  // Prim  PARTIELLEMENT pas vue
	       if( P[i][2]<0.0){ //
		 pastoutvu=true; 
		 //cout <<"Camera[calc_visi] Z neg\n";
	       }
               else {
		 (*pdiff)++; (*pdiff)++; (*pdiff)--;
		 //cout <<"Camera[calc_visi] P vu mais pastoutvu actif!\n";
	       }
             else {
	       Pp[i][0]=P[i][0]*foc/( P[i][2]+foc);
	       Pp[i][1]=P[i][1]*foc/( P[i][2]+foc);
	       Pp[i][2]=0.0;//P[i][2];
	       //cout <<"Camera[calc_visi] Pp{Re}["<<i<<"]  = ";Pp[i].show();
	       // Chgt d'origine E->O (coinBG)
               // distZ=Pp[i].dist2(P[i]);
	       Pp[i]-=EvI;
	       Pp[i][2]=P[i][2];//distZ;
	       cout <<"Camera[calc_visi] Pp{Rimage}["<<i<<"]  = ";Pp[i].show();
	     }//else P[i][2]<0.0) || pastoutvu)
	   }// for points triangle
	   if(!pastoutvu) {
	     // Tri sommets tq Pp[i][1]<<Pp[j][1]<<Pp[k][1] ie A[1] < B[1] < C[1]
	     j = (Pp[1][1]>Pp[2][1])? 1: 2; // calc intermed
	     k = (Pp[0][1]>Pp[j][1])? 0: j; // indice max pour coord y
	     i = (k+1)%3; j= (i+1)%3;
	     i = (Pp[i][1]<Pp[j][1])? i : j; // indice min pour coord y
	     j = 3- i-k;
	     //       cout<<" (i,j,k) = "<<i<<j<<k<<endl;
	     if ((i!=k)&& !((A[0]==B[0])&&(B[0]==C[0]))&& !((A[1]==B[1])&&(B[1]==C[1]))) {
	       // Pts A,B,Cpas  alignes selon les axes Xou Y
	       // tri Ok
	       up=down=false;
	       if (A[1]==B[1]){ // up 
	        i = (A[0] <B[0])?i:j;
	       j = 3-i-k; 
                 up=(A[0]==B[0])?false: true;
	       }//if up
	       else {
		 if(B[1] == C[1]){ // down 
	           k = (B[0] <C[0])?k:j; 
		   j = 3-i-k;
		   down =(B[0]==C[0])?false: true;
		 }//if down
                 else { 
		   D[1] = B[1];
		   pente=(D[1]-A[1])/(C[1]-A[1]);
		   D[0] = pente*(C[0]-A[0])+A[0];
                   // D[2] = (A[2]*(C[1]-D[1]) + C[2]*(D[1]-A[1]))/(C[1]-A[1]);
	           // D[2] = pente*(C[2]-A[2])+A[2];
                   D[2] = 0.0;
		   // calcul du vrai Z de D
                   lpt.ajoute(P[i]);
                   lpt.ajoute(P[j]);
                   lpt.ajoute(P[k]);
                   pprim= new Triangle(lpt,-1);
                   lpt.free_liste();
                   Pt1=oeil;
                   Pt1-=SvE;
                   Pt1=Pt1.chgt_base(u,v,w);
                   Pt2=D+EvI;
                   Vecteur dir(Pt1,Pt2);
                   dir.normalise();
                   parag.change_origine(Pt1);
                   parag.change_direction(dir);
		   cout <<"param k de l'intersection = "<< pprim->intersect(parag,&Pt1) <<endl;
                   delete pprim;
                   cout<<" D = ";Pt1.show();
                   D[2]=Pt1[2];
		   up=down=(B[0]==D[0])?false: true;
		   if(D[0]>B[0]) { l=i; i=j; j=3;}
		   else        { l=i; i=3;     }
		 }//else cas quelconque, ni up , ni down 
	       }
	       //     cout<<"2 (i,j,k) = "<<i<<j<<k<<endl;
	       if(up)
		 colorie_triangle(pdiff,C,A,B, Zbuf,Zprim,&imgvisi);
	       if(down) {
		 if (up) {
		   k=j; j=i; i=l;
		 }
		 colorie_triangle(pdiff,A,B,C, Zbuf,Zprim,&imgvisi);
	       }// if down
	     }//if pas un triangle plat
	   }//if !pastoutvu 
  }// for liste diffuseurs
  imgvisi.sauve(); 
}//Canopy::calc_FF()
#endif



#ifdef _Obsolete
void Canopy::calc_FF_Bfar(SPMAT *FF,VEC **Cfar,char * envname, double denv,int nbsim) {
  register int p,t,b,nb_rec=0,nb_emi=0,nb_test=0,cum_box=0;
  register unsigned char i,nb_box,pr=0;
  int pmax,tmax,inc[3],Gi[3];
  SPROW	*r_sup,*r_inf;
  int n,i_sup,i_inf;
  Diffuseur *diffR, *diffE;
  BSP *box,*Tabox[8];
  Point G,T;
  double mid,S[3],dGS[3],d2env=denv*denv,ff; //attention aux tests entre d et d^2
  //ListeD<Diffuseur *> L2;//permet d'elimnier les doublons
  bool rechauffe;
  
  printf("denv=%lf - taille_vox = %lf\n",denv,mesh.taille());
  init_NFF(envname);
  for(Ldiff.debut();! Ldiff.finito();Ldiff.suivant()){
    diffR=Ldiff.contenu();
    //pmax=diffR->nb_patch();
    //printf("[calc_FF_Bfar] diffR = %d\n",diffR->num());
    //printf("%d ",diffR->num());
    //pr++;
    //if((pr%10)==0) printf("\n");
    //chron.Start();
    //proj_cpt=0;
    //printf("* pmax = %d\n",pmax); 
    //for(p=0;p<pmax;p++) {//cas du diffuseur patche
    //diffR->select_patch(p);
    diffR->active(1);//active la face sup
    init_proj(diffR);
    G=diffR->centre();
    i_sup=diffR->num();
    printf("\n %d -",i_sup);
    r_sup = FF->row+ i_sup;
    if(!diffR->isopaque()) {
      diffR->togle_face();
      i_inf=diffR->num();
      // printf(" %d ",i_inf);
      r_inf = FF->row+ i_inf;
      diffR->togle_face();
    }
    //#define _latotale
    //#ifndef _latotale
    //genere la liste de voxels utiles
     nb_box=0;
    for(i=0;i<3;i++) {
      //Gi[i]=mesh.coord(i,G[i]);
      Gi[i]=(G[i]-mesh.origine()[i])/mesh.taille();
      //printf("i=%d - G[i] = %lf - Gi[i]=%d\n",(int)i,G[i],Gi[i]);
      mid=mesh.milieu_vox(i,Gi[i]);
      if(G[i]<mid && Gi[i]!=0) {
	inc[i]=-1;
	S[i]=mesh.sommet(i,Gi[i]);
      }
      else
	if(G[i]>mid && Gi[i]!=mesh.nb_voxel(i)-1) {
	  inc[i]=+1;
	  S[i]=mesh.sommet(i,Gi[i]+1);
	}
	else {
	  inc[i]=0;
	  S[i]=0.0;
	}
      if(inc[i]!=0 && fabs(G[i]-S[i])>denv) {
	inc[i]=0;
	//cout<<(int)i<<" - "<<nb_rec<<" - "<<pmax<<"bordel acqueux!\n";
      }
    }//loop sur les axes
    nb_box=addbox(mesh(Gi[0],Gi[1],Gi[2]),Tabox,nb_box);
    char axe[3]= {2,2,2},diag=3;
    //axe -> 0 : x-y, 1 : y-z, 2 : z-x
    if(inc[0]!=0) {
      nb_box=addbox(mesh(Gi[0]+inc[0],Gi[1],Gi[2]),Tabox,nb_box);
      axe[0]--; axe[2]--; diag--;
    }
    if(inc[1]!=0) {
      nb_box=addbox(mesh(Gi[0],Gi[1]+inc[1],Gi[2]),Tabox,nb_box);
      axe[0]--; axe[1]--; diag--;
    }
    if(inc[2]!=0) {
      nb_box=addbox(mesh(Gi[0],Gi[1],Gi[2]+inc[2]),Tabox,nb_box);
      axe[1]--; axe[2]--; diag--;
    }
    // precacul distance d(G,S)
    dGS[0]=(G[0]-S[0])*(G[0]-S[0]);
    dGS[1]=(G[1]-S[1])*(G[1]-S[1]);
    dGS[2]=(G[2]-S[2])*(G[2]-S[2]);
    //test des boites jointives par un axe 
    if(axe[0]==0) {
      //test/ dist(G, axe x-y ie z)
      if( (dGS[0]+dGS[1]) < d2env)
	nb_box=addbox(mesh(Gi[0]+inc[0],Gi[1]+inc[1],Gi[2]),Tabox,nb_box);
    }
    if(axe[1]==0) {
      //test/ dist(G, axe y-z ie x)
      if( (dGS[1]+dGS[2]) < d2env)
	nb_box=addbox(mesh(Gi[0],Gi[1]+inc[1],Gi[2]+inc[2]),Tabox,nb_box);
    }
    if(axe[2]==0) {
      //test/ dist(G, axe z-x ie y)
      if( (dGS[2]+dGS[0]) < d2env)
	nb_box=addbox(mesh(Gi[0]+inc[0],Gi[1],Gi[2]+inc[2]),Tabox,nb_box);
    }
    //test de la  boite jointive par le sommet S
    if(diag==0) {
      //test/ dist(G, axe y-z ie x)
      if( (dGS[0]+dGS[1]+dGS[2]) < d2env)
	nb_box=addbox(mesh(Gi[0]+inc[0],Gi[1]+inc[1],Gi[2]+inc[2]),Tabox,nb_box);
    }
    cum_box+=nb_box;
      
     
    //allez zou, on calcule les FF et les Bfar par projection
    //printf("ifndef _latotale\n");
    for(b=0;b<nb_box;b++) {//loop sur les boites dans denv
      box=Tabox[b];
      //loop sur les diffuseurs de la boite
      for(box->Ldiff.debut();! box->Ldiff.finito();box->Ldiff.suivant()){
	diffE=box->Ldiff.contenu();
	/*
	  #ifdef _AntiDoublon
	  //ramene le tps a celui de toute la scene pour le tps de gestion de boucle
	  // a gerer par la matrice en testant si la M(diffR,diffE) est deja prise? tq pas doublon (a voir le temps de parcours de la matrice creuse
	  //printf(" ** card(L2)=%d\n",L2.card());
	  rechauffe=false;
	  if(!L2.est_vide()) {
	  L2.deb_pause();
	  for(L2.debut();! L2.finito();L2.suivant()){
	  //printf(" card(L2)=%d\n",L2.card());
	  if(diffE==L2.contenu()) {
	  rechauffe=true;
	  //printf("anti-doublon marche!\n");
	  break;
	  }
	  }
	  L2.fin_pause();
	  }
	  if(!rechauffe){
	  L2.ajoute(diffE);
	  #endif      //AntiDoublon
	  #else //ifdef _latotale
	  //loop sur les diffuseurs de la scene
	  Ldiff.deb_pause();
	  for(Ldiff.debut();! Ldiff.finito();Ldiff.suivant()) {
	  //printf("ifdef _latotale\n");
	  diffE=Ldiff.contenu();
	    
	  #endif //ifdef _latotale
	  */
	if (diffE!=diffR) {//Diagonale : FF=1
	    
	  //cas du diffuseur patche
	  //tmax=diffE->nb_patch();
	  //for(t=0;t<tmax;t++) {
	  //diffE->select_patch(t);
	  T=diffE->centre();
	  nb_test++;
	  if(G.dist2(T)< d2env) {
	    Vecteur dir(G,T);
	    diffE->active(dir);
	    n=diffE->num();
	    if( sprow_idx(r_sup,n)<0) {
	      //Anti-Doublon System
	      sp_set_val(FF,i_sup,n,-1.0);
	      proj_ortho(diffE);
	      nb_emi++; proj_cpt++;
	    }//if pas deja traite
	  }//if diffE dans l'envt
	  //}//for nb patch diffE (t)	 
	}//if diffE != diffR
      }//for box->Ldiff 
      //#ifndef _latotale
    }//jqa nb_box(b) (ou if rechauffe si AntiDoublon)
    /*#ifdef  _AntiDoublon
      }//jqa nb_box
      if(!L2.est_vide())
      L2.free_liste();
      #endif //	_AntiDoublon
      */
    //#else  //ifdef _latotale
    //Ldiff.fin_pause();
    //#endif //ifdef _latotale
    nb_rec++;
    //Calcul des FF en fct des Buffers
    //chron.Stop();
    //cout<<proj_cpt<<" triangles projete's en "<<chron<<endl;
    //chron.Start();
    NFF(FF,i_sup,i_inf,Cfar,denv);
    //chron.Stop();
    //cout<<proj_cpt<<"NFF en "<<chron<<endl;
    //}//for nb patch diffR
  }//for Ldiff
  stat_NFF();
  printf("[calc_FF_Bfar()] \n  mean(diff/patch) = \t%lf\n  mean(diff_proj/patch) = \t%lf\n mean(box/patch) = \t%lf\n",nb_test/(double)nb_rec,nb_emi/(double)nb_rec,cum_box/(double)nb_rec);
}//Canopy::calc_FF_Bfar()

#endif

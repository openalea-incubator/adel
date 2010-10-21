/************************************************************
*                 radioxity.C - MC 96, 98                   *
*       main() de la radiosite mixte                        *
*                                                           *
* MC Oct05: modif pour Caribu4.4 = Caribu -> Esup et Einf   
* MC june 08: maj Etri.vec pour version Windows
* MC spet 09: cretaion du fichier Etri.vec0 qui compatible avec le .can d'entree
*************************************************************/

#include <iostream> // introduire la notion de namespace
using namespace std ;


#include <ferrlog.h>
#include <system.h>     // raytools::include::bibliotek

#ifdef WIN32
#include <windows.h>	// Mem partagée via CerateFileMapping/MapViewOfFile
char pcClefNum[12] ;	// pour transcription itoa de "clef" signée si y faut
HANDLE	hSharedSeg ;	// Handle du fichier mappé
LPVOID	lpSharedSeg ;	// pointeur LPVOID sur seg. partagé

#include <assert.h>	// associé à detection d'erreurs IO Windows
#endif

#include <outils.h>
#include <chrono.h>
#include <canopy.h>
//#include "lumiere.h"
#include "Mmath.h"

extern "C" {
#include "sparse.h"
#include "iter.h"
}
#ifdef _HD
#include "solver.h"
#include "bzh.h"
#endif

bool bMemoriseMatrix=false;

#include "GetOpt.h"

#define PAUSE(msg)  printf(msg); Ferr <<"- Taper la touche Any" ; getchar() ;

//Prototypes des fct locales
static  void beep(const char *,int);
static  void erreur_syntaxe(char *);
static int options(int argc,char **argv);
static  void genres();

// Variables globales 
extern unsigned int NB;
char verbose=2;

// Variables locales a radioxity.C
static unsigned int i,j;
static int tog, sol;
static FILE * fres;
static Diffuseur **TabDiff,*diff;
static Canopy scene;
static VEC  **B0,**B, **Cenv;
static char opak;
//  Options
static  unsigned int nb_iter,nbsim;
static double denv;
static  bool ffseul, infty, geom, ordre1, 
  ff_print, bio, byseg, byfile, radonly, memsize,bias;
static  double seuil;
static  char *maqname, *envname, *optname, *lightname, *name8; 
static   int clef_shm=-1;
static  char *dirname, *matname;
// Option capteur virtuel - MC0699
static  bool solem; 
static char * nsolem;

ferrlog Ferr((char*)"canestra.log") ;
#ifndef NOMAIN
int main(int argc,char **argv){
#else
  int canestra (int argc,char **argv){
#endif
    reel bornemin[3]={99999999.0,99999999.0,99999999.0};
    reel bornemax[3]={-99999999.0,-99999999.0,-99999999.0};
    Chrono chrono;
    
    // Arret en cas d'erreur sur la ligne de commande
    if(options(argc,argv)) {
      Ferr << "Erreur sur la ligne de commande\n" ;
      return -1; 
    }

    // ? Calcul des flux moyens (eg. sail) par un exec()?
    // ? traitement des surfaces trop grandes (sol, tiges) ?
  	
    //*********** Chargement de la scene et evt des capteurs virtuels ********
    chrono.Start();
    if(byfile){
      scene.parse_can(maqname,optname,name8,bornemin, 
		      bornemax,sol,nsolem,TabDiff);
      Ferr<<__FILE__<<" : byfile"<<'\n';
    }else{
      scene.read_shm(clef_shm,optname,name8,bornemin, 
		     bornemax,sol,nsolem,TabDiff); 
      Ferr<<__FILE__<<" : byshm"<<'\n';
    }
    cout <<"\n Nombre de faces (2*L+S) = "<<scene.radim<<endl;cout.flush();
    if(true ||verbose) {
      for(j=0; j<3; j++)
	//  cout<<" main() : Bornemin["<<j<<"] = "<<bornemin[j]
	//  <<"  Bornemax["<<j<<"] = "<<bornemax[j]<<endl;
	Ferr<<" main() : Bornemin["<<j<<"] = "<<bornemin[j]
	    <<"  Bornemax["<<j<<"] = "<<bornemax[j]<<'\n';
    }
    if(denv<0){// Full-matrix case
      Ferr<<__FILE__<<" : Full Matrix"<<'\n';
      denv=0;
      double x;
      for(j=0; j<3; j++){
	x=bornemax[j]-bornemin[j];
	denv+=x*x;
      }
    
      denv=sqrt(denv)/2.;
      Ferr <<"<*> Full-matrix case: denv<0 ==> denv="  << denv<<"\n" ;
    }
    chrono.Stop();
    Ferr<<"\n>>> Canestra[main] Scene chargee en "<<chrono<< '\n';
    // BUG?F(IPD)
    //fflush(stdout); 
    // fflush(stderr); 
  
  
    //*********** Construction de la double grille ***************
    //pas besoin de grille si reuse ou 1er ordre ou SAIL pur
  
    if(!(radonly || ordre1 || denv==0))  {
      chrono.Start();
      //ex-pb effet de bord avec projplan car modifie les bornes => vmin et vmax
      scene.cstruit_grille(denv);
      chrono.Stop();
      Ferr<<">>> Canestra[main] Grille construite en "<<chrono<< '\n' ;//endl;
    }
    //scene.mesh.visu();
  
    //************ Calcul de l'eclairage direct (soleil & ciel)  ***************
    //    initialisation
    double *Bsource;
    opak=0;
  
    chrono.Start();
    Bsource = new double[scene.radim];
    B= B0 = new VEC*[nbsim]; //B=B0 si pas de calcul des rediffusions
    for(i=0;i<nbsim;i++) {
      B0[i] = v_get(scene.radim);
      for(j=0;j<scene.radim;j++){
	B0[i]->ve[j]=0.0;
      }
    }
  
    //    calcul de visibilite (purely geometric)
    Vecteur dir_source;
    double Esource,rho;
    //     calcul de l'eclairage direct (soleil, ciel)
    ifstream flight(lightname,ios::in);
    do {
      flight>>Esource;
      if(!flight) 
	break;
      flight>>dir_source[0]>>dir_source[1]>>dir_source[2];
    
      for(i=0;i<scene.radim;i++) {
	Bsource[i]=0.0;
      }
      scene.projplan(dir_source,infty,Bsource);
      //recommenter
      Ferr <<"param. projplan : dir = ("  << dir_source[0]<<"," << dir_source[1]
	   <<","  << dir_source[2]<<") - Esun = "  << Esource<<'\n' ;
      for(i=0;i<scene.radim;i++) {   
	if(Bsource[i]!=0.0) {  
	  TabDiff[i]->activ_num(i);
	  if(Bsource[i]>0) 
	    rho=TabDiff[i]->rho()*Bsource[i]/TabDiff[i]->surface();
	  else { 
	    TabDiff[i]->togle_face();
	    rho=-TabDiff[i]->tau()*Bsource[i]/TabDiff[i]->surface();
	  }
	  // Cumule les contrib des differents angles solides
	  B0[0]->ve[i]+=Esource*rho;
	  //Ferr <<"i="  << i<<" : Bsource="  << Bsource[i]<<", B0="  
	  //     << B0[0]->ve[i]<<"\n" ;
	  TabDiff[i]->active(0);// reactive le diff sur la face sup (defaut)
	}
      }   
    }while(flight);
    flight.close();
    chrono.Stop();
    Ferr<<">>> Canestra[main] calcul du direct en "<<chrono<<'\n' ; 
  
    if(byfile){//ecriture du direct dans un fichier E0	
      fres=fopen("E0.dat","w");
      for(j=0;j<scene.radim;j++) {
	//Ferr <<"B0("  << j<<") ="  << B0[0]->ve[j]<<" - B("  << j<<") ="  
	//   << B[0]->ve[j]<<" \n" ;
	fprintf(fres,"%.10lf \n ",B0[0]->ve[j]);
      }
      fclose(fres);
    }//if(byfile)
    /* ******************************************************************** */
  
    //*********** Ecriture de resultats partiels  ***************
    if(geom) {
      Point G;
      Vecteur N;
      int istem;
      double teta,deg=180./M_PI,surfT;
      //Ferr << __FILE__<<" : "<< __LINE__<<'\n' ;
      fres=fopen("geom.dat","w");
      for(i=0;i<scene.radim;i++) {
	diff=TabDiff[i];
	if(TabDiff[i]->isopaque()){
	  istem=1;
	  if(TabDiff[i]->primi().name()==0)
	    istem=-1;
	}
	else{
	  istem=0;
	}
	diff->activ_num(i);
	surfT=diff->surface();
	G=diff->centre();
	N=diff->normal();
	teta=Macos(fabs(N[2]))*deg + ((N[2]<0)?90:0);
	fprintf(fres,"%d   %lf %lf %lf   %lf %d %.10lf\n",
		i,G[0],G[1],G[2],teta,istem,surfT);
	TabDiff[i]->active(0);
      } 
      fclose(fres);
    }
    //if geom
    /**************************************************************************/
	      
    if(!ordre1){// Calcul des rediffusions
      //Calcul des FF et des Bfar
      chrono.Start();
    
      VEC  *r0,*x;
      Cenv = new VEC*[nbsim];
      B = new VEC*[nbsim];
      for(i=0;i<nbsim;i++) {
	Cenv[i] = v_get(scene.radim);
	B[i] = v_get(scene.radim);
	for(j=0;j<scene.radim;j++)
	  Cenv[i]->ve[j]=0.0; }
      if(denv>0){
#ifdef _HD
	if(!radonly){//calcul de la matrice des FF
	  Ferr <<" Version longue : calcul des FF et des Coeff de Bfar"
	       <<'\n' ;
	  //fflush(stderr); // Ferr.flush vient d'etr appele

	  hdmat_init(dirname,matname);
	  scene.calc_FF_Bfar(Cenv,&Esource,envname,bias,denv,nbsim);
	}
	else{
	  //lecture de la mat. : maj des NzName, DgName et BfName
	  //et calcul des Bfar	
	  Ferr <<" Version courte : lecture des FF et des Coeff de Bfar\n" ;
	  hdmat_majname(dirname,matname);
	  if(envname!=NULL)
	    hd_calc_Bfar(Cenv[0],envname,TabDiff,Esource);
	}
#else
	SPMAT *FF;
	r0= v_get(scene.radim);
	FF = sp_get(scene.radim,scene.radim,200);
	//Ferr << __FILE__<< " : "<< __LINE__ << '\n' ;
	for(i=0;i<scene.radim;i++)
#ifndef BCC32
	  r0->ve[i]=drand48();
#else
	r0->ve[i]=rand() / (double) RAND_MAX ;
#endif
	if(verbose>2) 
	  Ferr <<"\n--> initialiastion des vecteurs et matrices faites (alloc.)"
	       <<'\n' ;
	//PAUSE(" ");
	if(ffseul) {
	  //Ferr << __FILE__<< " : "<< __LINE__ << '\n' ;
	  scene.calc_FF(FF);
	  //Ferr << __FILE__<< " : "<< __LINE__ << '\n' ;
	}else {
	  //Ferr << __FILE__<< " : "<< __LINE__ << '\n' ;
	  scene.calc_FF_Bfar(FF,Cenv,&Esource,envname,bias,denv,nbsim);
	  //Ferr << __FILE__<< " : "<< __LINE__ << '\n' ;
	}
#endif
	chrono.Stop();
	Ferr<<">>> Canestra[main] FF et Bfar calcules  en "<<chrono<< '\n';
	//fflush(stderr);
      
	/****************************************************/
      
	// Resolution du systeme lineaire
	int num_steps;
#ifdef _HD
	if(ff_print) {
	  //Ferr << __FILE__<< " : "<< __LINE__ << '\n' ;
	  char carlu;
	  Ferr <<"* ATTENTION - mode HD - Print la matrice a l'ecran\n"
	    " Le voulez vous reellement \?(o/n)\?" ;
	  carlu=getchar();
	  if(carlu=='o' || carlu=='O')
	    print_hd_mat(TabDiff);
	}
	// FF -> syst. lineaire 1-Xi*Fij
	if(verbose>1)  
	  Ferr <<"==>Bi=E+Bfar\n" ;
	for(i=0;i<scene.radim;i++) 
	  B0[0]->ve[i]+=Cenv[0]->ve[i];
      
	// Solve Ax=b; B precondionneur, tol seuil, limit nb_iter_max
	chrono.Start();
	Ferr <<" MGCR-HD : seuil de cvgence = "  << seuil
	     <<" - nb_iter_max = "  << nb_iter<<"\n " ;
      
	//print_hd_mat(TabDiff);
	hd_mgcr(B[0],B0[0],TabDiff,seuil,100, nb_iter, &num_steps);

#else
	if(ff_print) {
	  FILE * fff;
	  fff=fopen("FF.dat","w");
	  for(i=0;i<FF->m;i++) {
	    for(j=0;j<FF->n;j++) {
	      fprintf(fff,"%lf  ",sp_get_val(FF,i,j));
	    }
	    fprintf(fff,"\n");
	  }
	  fclose(fff);
	}//if ff_print  
	// FF -> syst. lineaire 1-Xi*Fij
	SPROW	*r;
	unsigned int  len;
	double trans,refl,*pval;
	opak=0;
	for(i=0;i<FF->m;i++) {
	  B0[0]->ve[i]+=Cenv[0]->ve[i];
	  r = FF->row+i;
	  len=r->len;
	  TabDiff[i]->activ_num(i);
	  refl= -TabDiff[i]->rho();
	  if(!TabDiff[i]->isopaque()) {
	    TabDiff[i]->togle_face();      
	    trans=TabDiff[i]->tau();
	  }
	  for(j=0;j<len;j++) {
	    pval=&(r->elt[j].val);
	  
	    //Ferr <<" FF("  << i<<",ndx "  << j<<") = "  << *pval<<"\n" ;
	    if(*pval!=1) {
	      if(*pval>0)
		*pval*=refl;
	      else
		*pval*=trans;
	    }
	  }
	  TabDiff[i]->active(0);
	}
	if(ff_print) {
	  FILE * fff;
	
	  fff=fopen("M.dat","w");
	  for(i=0;i<FF->m;i++) {
	    for(j=0;j<FF->n;j++) {
	      fprintf(fff,"%lf  ",sp_get_val(FF,i,j));
	    }
	    fprintf(fff,"\n");
	  }
	  fclose(fff);
	}
	//sparse matrix resolution (Conjugate Gradient)
	// Solve Ax=b; B precondionneur, tol seuil, limit nb_iter_max
	//Methode  Leyk's MGCR
	chrono.Start();
	Ferr <<" MGCR : seuil de cvgence = "  << seuil
	     <<" - nb_iter_max = "  << nb_iter<<"\n " ;
      
	B[0]=iter_spmgcr(FF, (SPMAT *)NULL, B0[0],seuil, B[0],
			 20, nb_iter, &num_steps);
#endif
	chrono.Stop();
	if(num_steps<nb_iter)
	  Ferr <<" MGCR CONVERGE en "  << num_steps<<" iteration(s) \n" ;
	else
	  Ferr <<" MGCRN'A PAS CONVERGE' ! \n" ;
	Ferr<<">>> Canestra[main] Resolution du SL par MGCR en "<<chrono<<'\n';
      }//if denv<>0

      /************************************************************************/
      else{//SAIL pur
	//Ferr << __FILE__<< " : "<< __LINE__ << '\n' ;
	scene.sail_pur(Cenv,&Esource,envname);
	for(i=0;i<scene.radim;i++) 
	  B[0]->ve[i]=B0[0]->ve[i]+Cenv[0]->ve[i];
      }//if denv==0 ie SAIL pur
    
      /**********************************************************************/
    
    }//if calcul des rediffusions
  
    //Rendu - Traitement des resultats
    genres();
    // Gestion des fichiers persistants
    if(bMemoriseMatrix==false) {
      EffaceMatrices();
    }
    Ferr <<"This is the end...\n";
    Ferr.close();
    return 0 ;
  } //main()



  /*****************************************************************************
   **********               Fonctions Locales                          *********
   *****************************************************************************/

  //======>  genres(): calcule et genere les fichiers de resultats - MC98 
  void genres(){
    // Impression des resultats : vecteur des  radiosites, if(bio) Eabs.dat et Einc.dat
    int nbf=scene.radim-scene.nbcell;
  
    if(false && !ordre1){// genere les fichiers .dat de debug B0 et Bf generes
      if(envname != NULL){
	fres=fopen("Bf.dat","w");
	for(i=0;i<nbf;i++) 
	  fprintf(fres,"%lf \n",Cenv[0]->ve[i]);
	fclose(fres);
      }
      fres=fopen("B0.dat","w");
      for(j=0;j<nbf;j++) {
	fprintf(fres,"%.10lf \n ",B0[0]->ve[j]);
      }
      fclose(fres);
    }// if fichiers .dat de debug B0 et Bf generes
    // Ecriture des radiosites totales => B.dat
    if(byfile){
      fres=fopen("B.dat","w");
      Ferr <<"==> Impression des résultats radim="  << scene.radim<<", nbcell="  << scene.nbcell<<"\n" ;
      for(j=0;j<nbf;j++) {
	fprintf(fres,"%.10lf \n ",B[0]->ve[j]);
      }
      fclose(fres);
    }
    // Ecriture des ecliarement des capteurs virtues => solem.dat
    if(scene.nbcell>0){
      //id 1er ordre Total en eclairement et surface
      fres=fopen("solem.dat","w");
      for(j=0;j<scene.nbcell;j++) {
	fprintf(fres,"%.0lf\t %.10lf\t %.10lf \t%.6lf\n",
		TabDiff[nbf+j]->primi().name(),
		B0[0]->ve[nbf+j],B[0]->ve[nbf+j],
		TabDiff[nbf+j]->primi().surface());
	if(0) 
	  Ferr <<"SOLEM: " << j<<"/" << scene.nbcell<<" radim="  << scene.radim
	       <<", j+nbf="  << j+nbf<<", B0="  << B0[0]->ve[nbf+j]<<"\n" ;
      }
      fclose(fres);
    }//if nbcell>0
    //Calcul et impression des Eabs et Einc
    if(bio){
      reel *Ei,*Eabs,D,D0,D1,r0=0,r1=0,t0=0,t1=0,E0,E1,Esol,Ssol;
      int ia,shmid2=0;
      FILE *fa=NULL,*fi=NULL,*ft=NULL,*ft0=NULL;
      double *Te=NULL,surf, nom; 
      int Nt; int Nt0=0;
      if(byfile) {//by file
	fa=fopen("Eabs.vec","w");
	fi=fopen("Einc.vec","w");
	ft=fopen("Etri.vec","w");    
	fprintf(ft,"# canestrad: can=%s F8=%s opt=%s light=%s : denv=%.2f direct=%d \n",maqname,name8,optname,lightname,denv,(int)ordre1 );
	fprintf(ft,"# label1 Area Eabs(E/s/m2) Ei(sup) Ei(inf) (Ex=surfacic density of energy <nrj/s/m2>)\n");
	// Version repreannt la liste initiale de triangle du .can pr PyCaribu
	ft0=fopen("Etri.vec0","w");    
	fprintf(ft0,"# canestrad: can=%s F8=%s opt=%s light=%s : denv=%.2f direct=%d \n",maqname,name8,optname,lightname,denv,(int)ordre1 );
	fprintf(ft0,"# No Label1 Area Eabs(E/s/m2) Ei(sup) Ei(inf) (Ex=surfacic density of energy <nrj/s/m2>)\n");
      
      }
      else{//by shared memory
      
	DecodeClefIn(&Nt,&clef_shm,clef_shm); //In caribu
	// TEST: Ne plante plus si ouvre 1 segment =!= de celui de caribu
	//clef_shm+=1 ; 
	cout <<"==> print_Eabs(): Nt="  << Nt<<", clef_shm="  << clef_shm<<"\n" ;
	Ferr <<"==> print_Eabs(): Nt="  << Nt<<", clef_shm="  << clef_shm<<"\n" ;
#ifndef WIN32
	// Mode Unix
	shmid2=shmget((key_t)clef_shm,SEGSIZE*sizeof(double) ,IPC_CREAT|0666);
	if(shmid2==-1){//en cas de pb
	  // stderr2cerr: Parse error here ?
	  // Found _1_ formats but _0_ printable arguments
	  // Normal, y'en n'avait pas 
	  Ferr <<"<!> Ouverture du segment partage no. " <<clef_shm
	       <<" impossible => I terminate now !!\n";
	  exit(16);
	}
	Te=(double *)shmat(shmid2,0,(int)NULL);
#else
	sprintf ( pcClefNum, "%d", clef_shm) ;
	Ferr<<"------------------------oOo-------------------------"<<'\n';
	Ferr<<"Opening shared memory segment."<<'\n';
	Ferr <<"==> print_Eabs(): Nt="  << Nt<<", pcClefNum="  << pcClefNum<<'\n';
	Ferr<<"------------------------oOo-------------------------"<<'\n';
	// Mode Win NT
	assert ( ( hSharedSeg =
		   OpenFileMapping( 
				   FILE_MAP_ALL_ACCESS,
				   FALSE,
				   pcClefNum) ) != INVALID_HANDLE_VALUE ) ;
      
	assert (( lpSharedSeg = 
		  MapViewOfFile (
				 hSharedSeg,
				 FILE_MAP_ALL_ACCESS,
				 0,
				 0,
				 SEGSIZE*sizeof(double ))) != NULL ) ;
      
	Te=(double *)lpSharedSeg ;
#endif
      }
      Ei=new reel[nbf];
      Eabs=new reel[scene.Ldiff.card()-scene.nbcell+1]; // +1:HA

      opak=0;
      ia=0;
      Esol=Ssol=0;
      scene.Ldiff0.debut(); // ! Ldiff.finito();Ldiff.suivant()){
      for(i=0;i<nbf;i++) {
	//Geston de la sortie Etrivec0 identique a liste de triangle en entree - MC09
	while(scene.Ldiff0.contenu()>=0 ){
	  if(scene.Ldiff0.finito()) break;
	  fprintf(ft0,"%d %.0f 0 NaN NaN NaN\n",Nt0,scene.Ldiff0.contenu());
	  Nt0++;
	  // printf("dbg 2, Nt0=%d, Ldiff0()=%d\n", Nt0, scene.Ldiff0.contenu());
	  scene.Ldiff0.suivant();
	}  
	//gestion...
	diff=TabDiff[i];
	surf=diff->surface();
	nom=diff->name();

	if(opak==2)
	  opak=0;
	if(diff->isopaque())
	  opak=0;
	else
	  opak++;
	diff->activ_num(i);
	if(!(sol && diff->primi().name()==0)){
	  //Traitement
       
	  //Ferr <<"opak="  << (int)opak<<"; i="  << i<<", ia="  << ia
	  // <<", rho="  << diff->rho()<<", tau="  << diff->tau()<<"\n" ;
	  switch(opak){
	  case 0: // Opaque
	    if(diff->rho()==0){
	      Ferr <<"<!> Calcul de Einc d'un opaque corps noir impossible : \n"
		" r0*r1 == t0*t1" << '\n' ;
	      Ei[i]=Eabs[ia]=-1;
	    }
	    else{
	      Ei[i]=B[0]->ve[i]/diff->rho();
	      Eabs[ia]=Ei[i]-B[0]->ve[i];
	    }
	    if(byfile){
	      fprintf(fi,"%g\n",Ei[i]);
	      fprintf(fa,"%g\n",Eabs[ia]*surf);
	      fprintf(ft,"%.0f %f  %f  %f %f\n",nom, surf, Eabs[ia], Ei[i],-1.);
	      //liste compatible pycaribu - MC09  
	      fprintf(ft0,"%d %.0f %f  %f  %f %f\n",Nt0,nom, surf, Eabs[ia], Ei[i],-1.);
	      Nt0++;
	      scene.Ldiff0.suivant(); 
	    } else{
	      Te[ia]=Eabs[ia]*surf;
	      //MCoct05: caribu4.4
	      //met dans le SegMem les eclairement des faces sup et inf 
	      // Bug MC nov05: Te[ia+(Nt+1)]=B[0]->ve[i]*surf; //face sup
	      Te[ia+(Nt-1)]=Ei[i]*surf; //face sup
	      Te[ia+2*(Nt-1)]=-surf; //face inf
	      // Ferr <<"Te["  << ia<<"]="  << Te[ia]<<"\n" ;
	    }
	    //MCMarch2006
	    ia++;
	    break;
	  case 1://Transparent[face sup] => preparation
	    r0=diff->rho();
	    t0=diff->tau();
	    break;
	  case 2://Transparent[face inf] => resolution du syst 2eq, 2inc => E0 et E1
	    r1=diff->rho();
	    t1=diff->tau();
	    D=r0*r1 - t0*t1;
	    if(D==0){
	      Ferr <<"<!> Calcul de Einc d'un transparent impossible : \n r0("
		   << r0<<")*r1("<< r1<<") == t0("<< t0<<")*t1("<< t1<<")"<<'\n';
	    
	      Eabs[ia]=-1;
	      Ei[i-1]=Ei[i]=-1;
	      if(r0==t1)
		Eabs[ia]=B[0]->ve[i-1]*(1/r0-1)-B[0]->ve[i];
	    }
	    else{
	      D0= B[0]->ve[i-1]*r1 - B[0]->ve[i]*t1;
	      D1=r0*B[0]->ve[i] - t0*B[0]->ve[i-1];
	      Ei[i-1]=D0/D;
	      Ei[i]=D1/D;
	      Eabs[ia]= Ei[i-1]+Ei[i] - (B[0]->ve[i-1]+B[0]->ve[i]);
	      /* debug
		 Ferr  << r0<<"\t"  <<  t1<<"\t= "  <<  B[0]->ve[i-1]<<"\n" ;
		 Ferr  << t0<<"\t"  <<  r1<<"\t= "  <<  B[0]->ve[i]<<"\n" ;
		 Ferr <<"D0="  << D0<<", D1="  << D1<<", D="  << D<<" => E["  
		 << i-1<<"]="  << Ei[i-1]<<", E["  << i<<"]="  << Ei[i]
		 <<", Ea["  <<  ia<<"]="  <<  Eabs[ia]<<"\n\n" ;
		 /recommenter */
	    }
	    if(byfile){
	      fprintf(fi,"%g\n%g\n",Ei[i-1], Ei[i]);
	      fprintf(fa,"%g\n",Eabs[ia]*surf);
	      fprintf(ft,"%.0f %f  %f  %f %f\n",nom, surf, Eabs[ia], Ei[i-1], Ei[i]);
	      //liste compatible pycaribu - MC09  
	      fprintf(ft0,"%d %.0f %f  %f  %f %f\n",Nt0,nom, surf,  Eabs[ia], Ei[i-1], Ei[i]);
	      Nt0++;
	      scene.Ldiff0.suivant();
	    } else{
	      Te[ia]=Eabs[ia]*surf;
	      //MCoct05: caribu4.4
	      /* Bug 221105 MC
		 Te[ia+(Nt+1)]=B[0]->ve[i-1]*surf;//face sup 
		 Te[ia+2*(Nt+1)]=B[0]->ve[i]*surf; //face inf
	      */
	      Te[ia+(Nt-1)]=Ei[i-1]*surf;//face sup 
	      Te[ia+2*(Nt-1)]=Ei[i]*surf; //face inf
	   
	      if(verbose>2) {
		Ferr <<"Te["  << ia<<"]="  << Te[ia]<<" Ei(sup)="<<Ei[i-1]<<", Ei(inf)="<<Ei[i]<<"\n" ;
		cout <<"Te["  << ia<<"]="  << Te[ia]<<" Ei(sup)["<<ia+(Nt-1)<<"]="<<Ei[i-1]<<", Te(ia+2*(Nt-1)="<<Te[ia+2*(Nt-1)]<<", Ei(inf)["<<ia+2*(Nt-1)<<"]="<<Ei[i]<<", Nt="<<Nt<<"\n" ;
	      }
	    }
	    //MCMarch2006
	    ia++;
	    break;
	  }//switch      
	  //border effect
	  TabDiff[i]->active(0);
	}//if not soil appended
	else{// soil appended and soil primitive
	  /* Old version - Modif MC june08
	     Esol+= B[0]->ve[i]*surf/diff->rho();
	     Ssol+=surf;
	  */
	  Ei[i]=B[0]->ve[i]/diff->rho();
	  Eabs[ia]=Ei[i]-B[0]->ve[i];
	  if(byfile){
	    fprintf(fi,"%g\n", Ei[i]);
	    fprintf(fa,"%g\n",Eabs[ia]*surf);
	    fprintf(ft,"%.0f %f  %f  %f %f\n",nom, surf, Eabs[ia], Ei[i],-2.);
	  }
	  ia++;
	  Esol+=Ei[i]; 
	  Ssol+=surf;
	} 
      }//for nb_faces 
      //vidage de liste au cas ou - MC09
      if(!scene.Ldiff0.finito())
	while(scene.Ldiff0.contenu()>=0 ){
	  fprintf(ft0,"%d %.0f 0 NaN NaN NaN\n",Nt0,scene.Ldiff0.contenu());
	  Nt0++;
	  //printf("dbg 6, Nt0=%d, Ldiff0()=%d\n", Nt0, scene.Ldiff0.contenu());
	  scene.Ldiff0.suivant();
	  if(scene.Ldiff0.finito()) break;
	}
      //printf("dbg 7\n");    
  
      /* Cas du sol mis au placard - MC nov2005
	 if(Ssol>0){
	 // Eclairement du sol
	 Te[ia+(Nt+1)]=Esol/Ssol;
	 }else
	 // MC05 popur avoir toujours un nombre =(num_poly+1)*3
	 Te[ia]=-2;
	 Te[ia]=-2 ;//sol face sup 
	 Te[ia+2*(Nt+1)]=-3; //sol face inf
      */

      //Ferr << "Au max on atteint: Eabs["<<ia<<"]"<<'\n';

      if(byfile){
	fclose(fi); 
	fclose(fa);
	fclose(ft);
	fclose(ft0);
      } else
#ifndef WIN32
	// Unix way
	shmdt((void*)shmid2);
#else
      // Complicated Way
      UnmapViewOfFile(lpSharedSeg) ; // invalidation du ptr sur mem partagee
      CloseHandle(hSharedSeg) ;	   // Fermeture du fichier mappé
#endif
    
      delete [] Ei; delete [] Eabs;
    }//if bio
  
    beep("This is the end...",4);
  }//genres()


  //======>  beep(): fait bip !
  inline void beep(const char *msg="M'enfin ...",int nbeep=1){
    cout<<(char) 7 <<msg<<endl;
    for(register int i=1;i<1;i++) cout<<(char) 7<<endl;
  }//beep()

  //======>  erreur_syntaxe(): imprime les options du prog a l'ecran
  void erreur_syntaxe(char * prog){
    Ferr <<"Syntax Error:  the options of "  << prog<<" are \n" ;
    Ferr <<"  -M filename \t File describing the scene\n"	 
      "  -m shm_key\t Shared memory containing the scene\n"	
      "  -s Ns\t Append a soil to the scene (Ns is a treshold for the number of triangles)\n"
      "  -p filename \t File describingthe optical properties\n"
      "  -l filename \t File describing the light sources\n"	 
      "  -8 filename \t Infinite periodic canopy \n"	 "  -e filename \t Mean fluxes data (computed by SAIL) \n"
      "  -r Rsph \t Radius of the surrounding sphere \n"
      "  -d Dsph \t Diameter of the surrounding sphere \n"
      "  -F \t\t Print the form factors matrix \n"
      "  -R nb \t Resolution of the projection disk [52] \n"	
      "  -S nb \t Number of simulations [1] \n"
      "  -i nb \t Number of iteration of the CG solver [100]\n" 
      "  -a threshold \t Threshold of the CG solver [1e6] \n"
      "  -1 \t\t Compute only the direct lightning \n"
      "  -L nb \t Resolution of the light screen [1536]  \n"
      "  -A \t\t Generate  energy vector (Eabs.dat, Einc.dat)\n"
      "  -g \t\t Generate the geometry file (geom.dat)\n"
      "  -B \t\t Test the effect of the choice of inner triangles (bias?) \n"
      "  -T \t\t Estimate the maximum required memory\n"
      "  -v nb \t The level of verbose\n"
      "  -C filename \t File describing the virtual sensors\n"
#ifdef _HD   
      "  -f filename \t Simulate and store the matrix in filemane \n"
      "  -w filename\t Read the matrix file to simulate an other radiative case, without to compute form factors \n"
      "  -t dirname \t Name of the directory where the FF file is stored\n"
#endif	 
      "  -h \t\t This help message "<<'\n' ; // "%c",7);
  }//erreur_syntaxe()


  //======> options(): traite la ligne de commande argv - MC98
  int options(int argc,char **argv){
    int c;
    GetOpt option(argc,argv,"AC:BFTg1hs:L:M:R:S:8:a:d:e:f:i:l:m:n:p:r:t:v:w:");
  
    // Valeur par defaut des options
    NB=52; nb_iter=1000; nbsim=1;
    denv=0.30; seuil=1e-6; //-1 ie seuil_solver=MACHEPS
    ffseul=infty=geom=ordre1=ff_print=bio=byseg=byfile=radonly=memsize=solem=false;
    bias=true;
    lightname=maqname=envname=optname=name8=dirname=matname=nsolem=NULL;
    sol=0;
    scene.Timg=1536;
    // Traitememnt des options
    if(argc<2){erreur_syntaxe(argv[0]);return 1;}
    while((c=option())!=EOF)
      switch(c) {
      case 'A' : bio =true;                       break;// genere Eabs.dat et Einc.dat
      case 'B' : bias=false;                      break;// pb des a cheval sur la sphere  
      case 'C' : nsolem=option.optarg; solem=true;break;// solem.can     
      case 'F' : ff_print=true;                  break;// FF -> FF.dat
      case 'L' : scene.Timg=atoi(option.optarg); break;//Resolution projplan 
      case 'M' : maqname=option.optarg; byfile=true; break;//maquette .can
      case 'S' : nbsim=atoi(option.optarg);      break;// nombre de simulations  
      case 'R' : NB=atoi(option.optarg);      break;// Resolution FF
      case 'T' : memsize=true;                   break;// Appel maxmem> maxmem.res mem en Ko 
      case '1' : ordre1=true;                    break;//stop apres ordre 1
      case '8' : infty=true;name8=option.optarg; break;//infinity  
      case 'a' : seuil=atof(option.optarg);      break;// seuil de convergence
      case 'd' : denv=atof(option.optarg)/2.;    break;// diam de la sphere
      case 'e' : envname=option.optarg;          break;//donnees de l'envt (eg. sail)
      case 'f' : matname=option.optarg;
	radonly=false;
	bMemoriseMatrix=true;
	break;//calc les FF
      case 'g' : geom=true;                     break;
      case 'h' : erreur_syntaxe(argv[0]); return 1;
      case 'i' : nb_iter=atoi(option.optarg);    break;// nbre d'iterations
      case 'l' : lightname=option.optarg;        break;
      case 'm' : clef_shm=atoi(option.optarg);byseg=true; break;// by segmem clef 
      case 'p' : optname=option.optarg;          break;
      case 'r' : denv=atof(option.optarg);       break;// rayon de la sphere
      case 's' : sol=atoi(option.optarg);;       break;// ajoute un sol
      case 't' : dirname=option.optarg;          break;// specifie le dir des hd mat  ; defaut = /tmp
      case 'v' : verbose=(char) atoi(option.optarg);    break;// verbose
      case 'w' : matname=option.optarg;
	radonly=true;
	bMemoriseMatrix=true;
	break;
	// Without :  lit la matrice pour refaire des calculs en changeant les po ou le sun
      default  : erreur_syntaxe(argv[0]); return 1;
      }// switch
  
    if (clef_shm != -1){
      Ferr <<"-------------o clef_shm = "<<clef_shm<<" o---------------"<<'\n';
    }
    if(byseg&&byfile){
      Ferr <<"<!> Fatal error"  << (char)7<<"\n==> Canestra should be called "
	"with 3 necessary options:\n -M maqname -p optname -l lightname \n " ;
      return 1;
    }
    if(((maqname==NULL)&&byfile) ||((clef_shm==-1)&&byseg) || (lightname==NULL) || (optname==NULL)) {
      Ferr <<"<!> Fatal error"  << (char)7<<"\n==> Canestra should be called "
	"with 3 necessary options:\n"
	"    [-M maqname| -m shmkey]  -p optname -l lightname \n " ;
      return 1;
    }  
  
    if(byfile) cout <<"\n Fichier maquette  :: "<<maqname;
    if(byseg ) cout <<"\n SegMem  maquette  :: "<<clef_shm;
    cout <<"\n Fichier optique   :: "<<optname;
    cout <<"\n Fichier sources   :: "<<lightname; 
    cout <<"\n Seuil convergence :: "<<seuil;
    cout <<"\n Dist envt (rayon) :: "<<denv;
    if(denv==0) cout <<" ==> <!> SAIL pur";
    if(denv<0){
      cout <<" ==> <!> full-matrix Radiosity (not nested)";
      envname=NULL;
    }
    cout<<endl;
    if(matname!=NULL) {
      cout <<" Fichier matrice   :: "<<matname<<endl;
      Ferr <<" Fichier matrice   :: "<<matname<< '\n' ;
    }
    if(sol){
      Ferr <<" Avec Sol          :: oui (nbT < "  <<  sol<<" )\n" ;
      if(ordre1) {
	sol=0;
	Ferr <<"<!> l'option \"1er ordre\" desactive l'option \"sol\" \n" ;
      }
    }
    if(infty){
      cout<<   " Avec Infini       :: oui\n";
      if( denv<0){
	Ferr << "<!> Fatal error "<<(char)7 <<'\n'
	     <<"==> Canestra is called with 2 incompatible options:\n"
	     <<"sphere diameter negative (-d or -r) (that means "
	     <<"full-matrix radiosity) and infinity mode (-8)\n ";
      
	return 1;
      }
    }
    Ferr <<" Projection disk   :: "  << NB<<"x"  << NB<<"\n" ;
    if(!bias)
      Ferr <<" Option test du biais(?) des triangles a-cheval sur la sphere ACTIVEE\n" ;
    if(envname==NULL && !ordre1){
      cout<<" ! Attention :: environnement nul.\n";
      Ferr <<" ! Attention :: environnement is null."<<'\n';
      if( denv==0){
	// stderr2cerr: Parse error here ?
	// Found _0_ formats but _1_ printable arguments
	// Y en avait, deux lignes plus bas: (char)7 ==> sqeezed 
	Ferr << "<!> Fatal error"<< '\n' 
	     << "==> Canestra is called with 2 incompatible options:\n"
	     << "sphere diameter null (-d or -r) and no mean fluxes files (-e)" 
	     << '\n';
	return 1; 
      }
    }
  
    if(dirname==NULL) {
      dirname= new char[100]; strcpy(dirname,".\\");}
    if(memsize){
      char cmd[125];
      sprintf(cmd,"maxmem %s 1 > maxmem.res &",argv[0]);
      if(verbose)
	Ferr <<"Option -T :  Pour avoir le max de memoire occupee = "<<cmd<<"\n";
      system(cmd);
    }
    fflush(stdout); fflush(stderr); 
    //Fin Gestion des Options
    return 0;
  }//options()


/* ***********************************************************************
             Albert OLIOSO   &   Michael Chelle                              
             Mars 1994         Mars 1996 - Avril 97
	     
               Calcul des profils radiatifs                                    
             d'un couvert vegetal en couches
************************************************************************ */
#define _Multicou
#include "multicou.h"
#include <system.h>

// exportées dans multicou.h
double pi, rd;
int N;

// Initialise the Ferr stream 
ferrlog Ferr((char*)"mc-sail.log") ;

int main(int argc, char **argv){
  double clai, ctau, croo, sf, dz, hau, ros;
  int id,npo;
  register int i,j;
  FILE *fpar, *fout, *fenv, *fpvf ;
  ifstream fpar2;
  char line[200];

  //exchange data
  Msailin msailin;
  Limit limit;
  Profout *Tprofout,*Cprofout,profout;
  Mlayout *Tlayout;
  char riri,skyname[200];

  if(argc==1){
    riri=1;
    strcpy(skyname,"incident");
  }else{
    riri=0;
    strcpy(skyname,argv[1]);
  }

  // Facteur de conversion
  pi = atan(1.)*4.;
  rd = pi/180.;

  //mise en place d'une back-door ;-)  pour calculer la BRDF du couvert ss hot spot par mcsail
  // debut MC09
  ifstream fvisee;
  ofstream fbrf;
  
  bool brdf=false;
  int N_visee=1;
  double theta_v, phi_v, BRF[200] ; // en degree
  
  fvisee.open("VISEE",ios::in);
  if(fvisee.is_open()){
    brdf=true;
    fvisee>>N_visee;
    cerr<< "<!> Presence d'un fichier VISEE => activation du mode canopy BRDF - MC09\n Nb de visees = "<<N_visee<<endl;
    fbrf.open("sail_brf.dat", ios::out);
  }

  for (int iv=0; iv<N_visee;iv++){//boucle pour la BRDF du couvert - MC09
    if(brdf) {
      fvisee>>theta_v>>phi_v;
      msailin.tto=theta_v;
      msailin.psi=phi_v;
    } //if brdf     
    //fin MC09


    // Parametres du modele de transferts radiatifs
    // => CROPCHAR
    fpar=fopen("cropchar","r") ; 
    //BUG Scons: MC08 : -DNDEBUG Supprime la fonction assert() donc les fichiers n'étaient pas ouvert
    assert (fpar != NULL ); 
    fgets(line,200,fpar);
    sscanf(line,"%d",&(msailin.nbang));
    fgets(line,200,fpar);
    sscanf(line,"%d %lf",&N,&dz);
    fprintf(stderr,"CROPCHAR : N=%d, nbang=%d, dz=%lf\n",N,msailin.nbang,dz);
    fclose(fpar);
    // N => dynamic allocation of array
    N++;// N couches de vegetation + le sol

    if(iv==0) msailin.alloue(N);
    Tlayout=new Mlayout[N+1];
    Tprofout=new Profout[N+1];
    Cprofout=new Profout[N+1];
    //mise a zero pour pouvoir cumuler les directions de ciel
    for(i=0;i<N+1;i++){
      Cprofout[i].transdir = Cprofout[i].transdif = Cprofout[i].trans = 0;
      Cprofout[i].refdif  =  Cprofout[i].refdir   = Cprofout[i].absc  = 0;
    }
    // => LEAFAREA
    fpar2.open("leafarea",ios::in);  // Hope objets cope with errors
    sf=0;
    double x;
    printf("==> Nb couche N=%d\n",N);
    for(i=N-1;i>0;i--){
      fpar2>>id>>id;
      for(j=0;j<msailin.nbang;j++){
	fpar2>>x;
	msailin.f(i,j)=x*100;
      }
      fpar2>>msailin.bmu[i] >> msailin.bnu[i] >> msailin.l[i]; 
      sf+=msailin.l[i];
      printf("LAI[%d]=%f - sf=%f - f(i,18)=%f\n",i,msailin.l[i],sf,msailin.f(i,msailin.nbang-1));fflush(stdout);
    }//for N couches  
    fpar2.close();
    printf("LAI total = %f\n",sf);fflush(stdout);
    // => SPECTRAL
    fpar2.open("spectral",ios::in);     
    fpar2>> npo >> ros;
    if(npo!=N-1){
      fprintf(stderr, "<!>\tNumber of optical properties layers <> N number of layers\n"
	      "\t only the first line of spectral taken into acount %c\n",7);
      npo=1;
    }
    fpar2>> croo >>  ctau;
    // Bug MC Feb  2006
    // BUG  for(i=1;i<N;i++){
    for(i=N-1;i>0;i--){
      msailin.roo[i]=croo;
      msailin.tau[i]=ctau;
      printf(" roo[%d]=%2.2f\t tau[%d]=%2.2f\n",i,msailin.roo[i],i,msailin.tau[i]);
      if(i>1)
	fpar2>> croo >>  ctau;
    }
    fpar2.close();
    //BUG MCFeb2006: cas ou il n' y a pas de sol dans la maquette =raz de ros
    if(argc>2)
      ros=0.;

    // Cas d'un sol lambertien
    limit.RSdd=limit.RSsd=limit.RSdo=limit.RSso = ros;

    // Parametres du direct 
    fpar2.open(skyname,ios::in);
    if(!fpar2.is_open()){
      fprintf(stderr," File %s does'nt exist => %s aborted...%c\n",skyname,argv[0],7);
      return -1;
    }
    int js,nbs=1;
    double dir[3],Esource;
    for(js=0;js<nbs;js++){
      if(riri){// Version compatible A. Olioso et Riri => INCIDENT
	Esource=1;
	fpar2>> hau >> x >> limit.ed;  
	msailin.tts=90-hau;   
	printf("%lf -%f\n",hau,limit.ed);
      }
      else{// Compatible canestra file.light
	printf("--> Direction ciel no. %d\n",js);
	fpar2>>Esource;
	if(!fpar2.good()) break;
	nbs++;
	fpar2>>dir[0]>>dir[1]>>dir[2];
	msailin.tts=acos(fabs(dir[2])/sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]))/rd;
	printf("    Esource=%.2g, theta_source=%.2g [nbs=%d]\n",Esource,msailin.tts,nbs);
	limit.ed=0; 
      }    
      limit.es=1-limit.ed;
      // appel aux sous programmes de calcul des coefficients de Verhoef
      for(i=0;i<N+1;i++){
	Tprofout[i].transdir = Tprofout[i].transdif = Tprofout[i].trans = 0;
	Tprofout[i].refdif  =  Tprofout[i].refdir   = Tprofout[i].absc  = 0;
      }
      //msail.C
      msailad(msailin);
      //Ferr << "mlayer() : Debut"<<'\n;
      mlayer(msailin, Tlayout);
      //Ferr<< "mlayer() : Fin"<<'\n';
      //  Matrices de reflectance et de transmittance  de chaque couche        
      //Ferr<<"mprofil() : Debut"<<'\n';
      mprofil(limit,Tprofout,Tlayout);
      //Ferr<<"mprofil() : Fin"<<'\n';
      // Cumul des donnees par direction de ciel
      for(i=0;i<N+1;i++){
	// printf("Tprofout[%d].transdir=%.3g - ",i,Tprofout[i].transdir);
	Cprofout[i].transdir += Esource * Tprofout[i].transdir;
	Cprofout[i].transdif += Esource * Tprofout[i].transdif;
	Cprofout[i].trans    += Esource * Tprofout[i].trans;
	Cprofout[i].refdif   += Esource * Tprofout[i].refdif ;
	Cprofout[i].refdir   += Esource * Tprofout[i].refdir;
	Cprofout[i].absc     += Esource * Tprofout[i].absc;
	// printf("Cprofout[%d].transdir=%.3g\n",i,Cprofout[i].transdir);
      }
    }//for directions ciel
    fpar2.close();

    //Ecriture des resultats
    clai = 0;
    msailin.l[0] = 0.;

    //BUG Scons: MC08 : -DNDEBUG Supprime la fonction assert() donc les fichiers n'étaient pas ouvert
    if(!brdf){
      printf("=> Ecriture des resultats\n");
      fout=fopen("profout","w") ; //unit=2
      fenv=fopen("mlsail.env","w") ;//unit=15
      fpvf=fopen("proflux.dat","w"); ;//unit=16
      assert (fout != NULL ) ; //unit=2
      assert (fenv != NULL ) ;//unit=15
      assert (fpvf != NULL ) ;//unit=16
      
      fprintf(fenv,"%d  %lf\n",N-1,dz);
      printf("ic  alti  clai   l(ic)  pene   Fdown  Fup  Rf   Tf\n");fflush(stdout);
      for(i=1;i<N+1;i++){
	profout=Cprofout[i];
	clai += msailin.l[i-1];
	printf("%2d  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f\n",i+1,(i-1)*dz,sf-clai,msailin.l[i],profout.transdir,profout.transdif,Cprofout[i-1].refdif, msailin.roo[i],msailin.tau[i]);
	fprintf(fout,"%f  %f  %f %f %f %f %f %f\n", clai,msailin.l[i],profout.trans,profout.transdir,profout.transdif,profout.absc,profout.refdif,profout.refdir);
	fprintf(fenv,"%f  %f  %f\n",(i-1)*dz,profout.transdif,Cprofout[i-1].refdif);
	fprintf(fpvf,"%8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f\n", (i-1)*dz, sf-clai, msailin.l[i],profout.transdir,profout.transdif,Cprofout[i-1].refdif,Tlayout[i].tss);
      }
      fclose(fout);
      fclose(fenv);
      fclose(fpvf);
    }//if !brdf
    else{
      if(phi_v==0)
	theta_v=-theta_v;
      printf("BRDF: tv=%.1lf, phi=%.1lf %.3lf\n", theta_v, phi_v,  Cprofout[N-1].refdir,  Cprofout[N-1].refdif);
      fbrf<<theta_v << "\t"<< Cprofout[N-1].refdir << "\t" << Cprofout[N-1].refdif << endl;
    }// else if ! brdf
       
    delete [] Cprofout;
    delete [] Tprofout;
    delete [] Tlayout;

    printf("\n=> fin boucle visee: iv=%d; N_visee=%d\n",iv, N_visee);
  }//FIN boucle pour la BRDF du couvert - MC09
  if(brdf){
    fvisee.close();
    fbrf.close();
  }
  Ferr.close() ;
  return 0;
} //main()
#undef _Multicou

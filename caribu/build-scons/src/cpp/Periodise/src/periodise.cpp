/*******************************************************************
 ****                     Periodise                             ****
 ****    genere un .can, qui est un motif pour le periodique    ****
 ****                Michael Chelle - 1996-1997                     ****
 *******************************************************************/
// Compilation Solaris-x86
// g++ -w -O2 -I . -I ../bibliotek -o periodise  periodise.C ../binobj/canopy_io.o ../binobj/diffuseur.o  ../binobj/primitive.o ../binobj/actop.o ../binobj/outils.o ../binobj/anamatide.o ../binobj/T_geometrie.o ../binobj/grille.o ../binobj/boite.o -lm 

// Compilation Radia
//g++ -w -O2 -I . -I ../bibliotek -o periodise  periodise.C $EXE/canopy_io.o $EXE/diffuseur.o  $EXE/primitive.o $EXE/actop.o $EXE/outils.o $EXE/anamatide.o $EXE/T_geometrie.o $EXE/grille.o $EXE/boite.o -lm 

/*
#######################################
05/04/05 : ajout de l'option -p pour permettre de specifier 
le fichier de proprietes optiques
14/09/09 : gestionde primitives non correctes pour conserver le .can original (-> pyCaribu)
#######################################
/*/

#include <cstdio>
using namespace std;


#include "GetOpt.h"


#include "canopyL.h"


inline void beep(const char *msg="M'enfin ..."){
  cout<<(char) 7 <<msg<<endl;
}

int periodise(int argc,char **argv){
  // Gestion des options
  char *outname=NULL,*maqname=NULL,*name8=NULL,*optname=NULL;
  char ok=0;
  GetOpt option(argc,argv,"h8:o:m:p:"); 
  int c;
  while((c=option())!=EOF)
    switch(c){
    case '8' : name8=option.optarg;   ok++; break;//infinity  
    case 'o' : outname=option.optarg;       break;// out.can
    case 'm' : maqname=option.optarg; ok++; break;// canopy.can 
    case 'p' : optname=option.optarg;  break; // *.opt 
    case 'h' : cerr<<"usage : "<<argv[0]<<" [8<infty> o<out> m<in> p<opt>]\n"; return 0 ;
    default  : cerr<<"usage : "<<argv[0]<<" [8<infty> o<out> m<in> p<opt>]\n"; return -1;
    }// switch
  
  if(ok<2) {
    printf(" ==> les options -m  et -8 sont strictement necessaires. \n");
    exit(0);
  }//pas good opt
  FILE * fout,*fz;
  if(outname==NULL) 
    fout=fopen("motif.can","w");
  else
    fout=fopen(outname,"w");
  fz=fopen("Bz.dat","w");
  char i;
  Diffuseur * pdiff;
  Canopy scene(false, false,false,1);
  reel  bmin[3]={99999999.0,99999999.0,99999999.0};
  double minr;
  reel bmax[3]={-99999999.0,-99999999.0,-99999999.0};
  double D[2],delta[3],Gz,Px,Py,Pz;//delta x et y 
   
  int nbs;
  int bps[3];//Bon Pour le Service
  scene.parse_can(maqname,optname,bmin, bmax,false,name8);
  D[0]=bmax[0]-bmin[0];
  D[1]=bmax[1]-bmin[1];
  minr=(0.01*bmax[2]+bmin[2])/1.01;
  printf("==>  minr = %g --> ",minr);
  if(minr<0){
    delta[2]=-minr;
  }
  else
    delta[2]=0.;
  printf(" delta2 = %g\n",delta[2]);
  for(i=0; i<3; i++) 
    cout<<" [1]Bornemin["<<(int)i<<"] = "<<bmin[i]<<"  Bornemax["<<(int)i<<"] = "<<bmax[i]<<endl;
  for(scene.liste_diff_scene.debut();
      !scene.liste_diff_scene.finito();
      scene.liste_diff_scene.suivant()){
    delta[0]=delta[1]=0;
    pdiff=scene.liste_diff_scene.contenu(); 
    if( &(pdiff->primi()) == NULL){// maintient des pas bon triangles - MC09
      fprintf(fout,"p  1 -1 3 0 0 0  0 0 0  0 0 0\n");
    }else{
    Gz=pdiff->primi().centre()[2]+delta[2];
    fprintf(fz,"%g\n",Gz);
    if(!pdiff->isopaque())
      fprintf(fz,"%g\n",Gz);
    nbs=pdiff->primi().nb_sommet();
    fprintf(fout,"p  1 %.0f %d \t",pdiff->primi().name(),nbs);// fflush(fout);
    for (i=0; i<2; i++) {
      bps[i]=pdiff->primi().nb_in(bmin[i],bmax[i], i); 
      if(bps[i]==nbs)
	delta[i]=-D[i];
      if(bps[i]==-nbs)
	delta[i]=D[i];
      //printf("delta(%d) = %g\n",i,delta[i]);
    }
    //printf("nbs = %d\n",nbs);
    for(i=0;i<nbs;i++){
      Px=pdiff->primi().sommets(i)[0]+delta[0];
      Py=pdiff->primi().sommets(i)[1]+delta[1];
      Pz=pdiff->primi().sommets(i)[2]+delta[2];
      //printf(" %d => %g %g %g \n",i,Px,Py,Pz);
      // la ligne suivante generait un Bus error => scinde en 3 lignes et ca marche ??
      //fprintf(fout," %.6g %.6g %.6g  ",Px,Py,Pz);
     fprintf(fout," %.6g  ",Px);
     fprintf(fout," %.6g",Py);
     fprintf(fout," %.6g",Pz);
     
    }
   fprintf(fout,"\n");
    }//if pas bon triangle
  }//for Ldiff
  fclose(fout);
  fclose(fz);
 return 0;
}//main::periodise


/* varaible globale des nom de fichier generer par tempnam(dir, pref) */
#include <iostream>
using namespace std ;
#include <ferrlog.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "diffuseur.h"
extern "C" {
#include "sparse.h"
}
#include "outils.h"
#include <bzh.h>
#include <system.h>

char pcNzName[128];
char pcDgName[128];
char pcBfName[128];


// Test de la fonction ANSI remove
//MC2005
#include <errno.h>
void Tremove(char * buffer){
  Ferr<<">>> bzh.cpp: Tremove()- destruction par remove() du fichier "<<buffer<<'\n';
  
  if (remove(buffer)) {
    switch(errno) {
      case EACCES:
        Ferr<<">>> bzh.cpp: Tremove()- Fichier protégé contre l'écriture !"<<'\n';
        cout<<">>> bzh.cpp: Tremove()- Fichier protégé contre l'écriture !"<<'\n';
        break;
      case ENOENT:
        Ferr<<">>> bzh.cpp: Tremove()- Fichier non trouvé !"<<'\n';
        cout<<">>> bzh.cpp: Tremove()- Fichier non trouvé !"<<'\n';
        break;
      case EINVAL:
        Ferr<<">>> bzh.cpp: Tremove()- Caractères invalides pour un nom !"<<'\n';
        Ferr<<">>> bzh.cpp: Tremove()- Caractères invalides pour un nom !"<<'\n';
	break;
      default:
        Ferr<<">>> bzh.cpp: Tremove()- Erreur indéterminée !"<<'\n';
        cout<<">>> bzh.cpp: Tremove()- Erreur indéterminée !"<<'\n';
    }
  }
}//Tremove()

//MC2005 : Ajout de la destruction des fichiers matrice si option -f
void placenette(char *ficname){
  /*  char cmd[255];
      int k; 
      k=access(ficname,F_OK);
      fprintf(stderr,"=> placenette(%s) k=%d\n",ficname,k);
      if(k==0){
      Ferr<< " <!> Danger: fichier "<<ficname<<" existant => destruction!"<<'\n'; 
      //Version hard codé pour Win$
      sprintf(cmd,"del %s",ficname);
      system(cmd);
      }
  */
  // Comment effacer un fichier (ANSI)?  
  //   #include <stdio.h>
  //   int remove(const char *pathname);
  
  
  remove(ficname);
}

/** Efface la matrice diagonale, les non_zeros et les Bfar */
void PlaceNette(void) {
  placenette(pcNzName);
  placenette(pcDgName);
  placenette(pcBfName);
  // Tremove(".\\toto");
}

/** Gestion des donnees persistantes de Canestra
* proto declare dans bzh.h */
void EffaceMatrices(void) {
  PlaceNette();
}

void hdmat_init(char *dir,char *rootname){
  char tempo[128];
    
  if(rootname==NULL)
    strcpy(tempo,"xxx");
  else
    strcpy(tempo,rootname);
  
  //

  sprintf(pcDgName,"%s%s%s",dir, DG_NAME, tempo);
  sprintf(pcNzName,"%s%s%s",dir, NZ_NAME, tempo);
  sprintf(pcBfName,"%s%s%s",dir, BF_NAME, tempo);
  
  if(true || verbose) 
    Ferr <<"HDMatrices : NZ="  << pcNzName<<", DG="  << pcDgName
	 <<", BF=" << pcBfName<<'\n' ;
 // cout <<"HDMatrices : NZ="  << pcNzName<<", DG="  << pcDgName <<", BF=" << pcBfName<<'\n' ;

  // Destruction des fichiers si existants - MC05
  PlaceNette();
  Ferr << " hdmat_init() : fin\n";
  Ferr << '\n' ;//fflush(stderr);
}

void hdmat_majname(char *dir,char*suff){

  char tempo[128];
    
  if(suff==NULL)
    strcpy(tempo,"xxx");
  else
    strcpy(tempo,suff);
  
  sprintf(pcDgName,"%s%s%s",dir, DG_NAME, tempo);
  sprintf(pcNzName,"%s%s%s",dir, NZ_NAME, tempo);
  sprintf(pcBfName,"%s%s%s",dir, BF_NAME, tempo);

  if(true || verbose) {
    Ferr <<"HDMatrices : NZ="  << pcNzName<<", DG="  << pcDgName
	 <<", BF="  << pcBfName<<'\n' ;
    Ferr << '\n' ; //fflush(stderr);
    if(verbose>2)  
      printf(" hdmat_majname() : fin\n");
  }
}
//  hdmat_majname

void hd_calc_Bfar(VEC *Cenv,char *pcEnvName,Diffuseur ** TabDiff,double Eclt){
  int	i,is,j, Nc,iff,nbp,*diag=NULL;//i : indice prim, is indice face
  double rho[2],tau[2],po;
  float cl[2];
  FILE *fic;
  Tabdyn<double,2> Tenv;
  char transp;

  if(verbose>2) 
    Ferr<<"*  hd_calc_Bfar() : Debut"<<'\n';
  fic=fopen(pcEnvName,"r");
  fscanf(fic,"%d %lf",&nbp,&po);
  //printf("Nc = %d - dz = %lf\n",nbp,po);
  Tenv.alloue(nbp+1,2);
  for(i=0;i<=nbp;i++) {
    fscanf(fic,"%lf %lf %lf ",&po,&(Tenv(i,0)),&(Tenv(i,1)));
    if(verbose>1)printf("z=%lf : trans[%d] = %lf - ref[%d] = %lf\n ",po,i,Tenv(i,0),i,Tenv(i,1));
  }
  fclose(fic);
  
  fic=fopen(pcBfName,"rb");
  if(verbose>1)printf("\t-> Lecture de %s\n",pcBfName);
  fread(&Nc,sizeof(int),1,fic);
  if(verbose>1) 
    Ferr<<"Nc = "<<Nc;
  if(Nc!=nbp){
    Ferr <<" Erreur nombre de couche de "  << pcEnvName<<" et de " 
	 << pcBfName<<" differents\n" ;
    exit(3);
  }
  fread(&nbp,sizeof(int),1,fic);
  if(verbose>1) 
    Ferr<<"nb_prim = "<<nbp;
  //Lecture de BF_name et remplissage de Cenv
  is=0;
  for (i=0; i<nbp;i++){
    //init des variables diffuseur
    transp=!TabDiff[is]->isopaque();
    TabDiff[is]->activ_num(is);
    rho[0]= TabDiff[is]->rho();
    if(transp) {// A verifier ordre de Tau
      tau[0]=TabDiff[is]->tau();
      TabDiff[is]->togle_face();
      rho[1]= TabDiff[is]->rho();
      tau[1]=TabDiff[is]->tau();
      TabDiff[is]->togle_face();
    }
    TabDiff[is]->activ_num(is);
    for (j =0; j<=Nc; j++) {
      fread(cl,sizeof(int),2,fic);
      Cenv->ve[is]+=rho[0]*(cl[0]*Tenv(j,0)+cl[1]*Tenv(j,1));
      if(transp)
	Cenv->ve[is+1]+=tau[0]*(cl[0]*Tenv(j,0)+cl[1]*Tenv(j,1));
    }

    if(transp) {
      for (j =0; j<=Nc; j++) {
	fread(cl,sizeof(int),2,fic);
	Cenv->ve[is]+=tau[1]*(cl[0]*Tenv(j,0)+cl[1]*Tenv(j,1));
	if(transp)
	  Cenv->ve[is+1]+=rho[1]*(cl[0]*Tenv(j,0)+cl[1]*Tenv(j,1));
      }
      is++;
      Cenv->ve[is]*=Eclt;
    }
    Cenv->ve[is-1]*=Eclt;
    is++;
  }

  fclose(fic);

  if(verbose>2)  
    Ferr<<"*  hd_calc_Bfar() : Fin\n";
  Tenv.free();
}//hd_calc_Bfar()


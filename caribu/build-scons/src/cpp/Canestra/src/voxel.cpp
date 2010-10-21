#include <iostream>
using namespace std ;

#include "voxel.h"
// reglage de la subdivision adaptative en run-time

#include <fstream> //.h>

//#define SEUIL_NB_DIFF 2
int SEUIL_NB_DIFF= 20; // 10fort.big.maq 184 lines

template <class Type>
inline char T_niveau_egal(Type *val,Type *val_max) {
  return ( (    (val[0]==val_max[0])
		&& (val[1]==val_max[1])
		&& (val[2]==val_max[2])
		) ? 1 : 0
	   );
} 

/* Eleve 2 a la puissance n */
int deux_puiss_n(int n){
  int res=1;
  res= res<<n;
  return res;
}

//Voxel
Voxel::~Voxel(){
  char str[250];
  BSP **Tptr;
  int N=0;
  bool newpt;
  Tptr=new BSP*[SI.maxi()[0]*SI.maxi()[1]*SI.maxi()[2]];
  Ferr <<"~Voxel: debut"<< '\n' ;
  if(1) {// Debug MCHA 151203
    // Pb de desallocation si 1boite pour +sieurs voxels
    // (free de zone deja free -MC
    // Si ptr non null, desallouer PUIS affecter ptr <- NULL // HA
    // MC161203: ne marche pas, car le delete ne doit pas mettre a null 
    // MC161203 : gestion a la mano des desallocation par un tableau de pointeur
    // NB une facon plus propre serait de gerer dans BSP un nbre d'instance

    register int i,j,k,n; 
    for(i=0;i<SI.maxi()[0];i++)
      for(j=0;j<SI.maxi()[1];j++)
	for(k=0;k<SI.maxi()[2];k++){
	  newpt=true;
	  for (n=0;n<N;n++){
	    if(SI(i,j,k)==Tptr[n]){
	      newpt=false;    
	      break;
	    }
	  }
	  if(newpt==true){
	    Tptr[N]=SI(i,j,k);
	    N++;
	  }
	  SI(i,j,k)=NULL;
	}
    for (n=0;n<N;n++)
      delete Tptr[n];
    delete [] Tptr;
  }// if debug   
  Ferr <<"~Voxel: fin"<< '\n' ;
  if(verbose>2)  
    cerr <<"~Voxel() FIN\n" ;
}// ~Voxel()

BSP* Voxel::operator() (int i, int j, int k){
  // Ferr <<"Voxel::operator() : DEBUT\n" ;
  return SI(i,j,k);
  }

Voxel& Voxel::operator = (Voxel& g){
  SI=g.SI;
  taille_vox=g.size();
  O=g.origine();
  return *this;
}

void  Voxel::division(char  axe, unsigned char * niveau_sub_max,BSP* B){
  BSP *B1;
  BSP *B2;
  register int i;
  int indi;
  bool taille_min[3];
  double frontier;
  
  
  niveau_sub_max[axe]++;
  if(verbose>3) 
    Ferr <<"[division] DEBUT axe "  << (int)axe<<" : niveau subdiv = " 
	 << (int) niveau_sub_max[axe]<<"\n " ;
  fflush(stderr);
  for(i=0;i<3;i++) {
    taille_min [i]= fabs((B->maxi(i)-B->mini(i))- taille_vox)<1e-4;
    if(verbose>3)
      Ferr <<"--> delta["  <<  i<<"] = "  << B->maxi(i)-B->mini(i)<<"\n" ;
  }
  
  //Ferr <<"Voxel[division] DEBUT\n" ;
  if(verbose>3)
    Ferr <<"nb_diffuseurs de B= "  <<  B->nb_diffuseur()<<"\n" ;
  
  if (taille_min[0] && taille_min[1] && taille_min[2] ){//Arret si BSPfeuille==voxel
    if(verbose>3) 
      Ferr <<"*** Voxel[division] Rgt boite : Niv_Max - nb_diff = "  
	   << B->nb_diffuseur()<<"\n" ;
    numerotation(B);
  }
  else{
    if (B->nb_diffuseur() <= SEUIL_NB_DIFF){
      if(verbose>3)
	Ferr <<"Voxel[division] Rgt boite : Nb_Max\n" ;
      numerotation(B); 
    }
    else{
      B1=new BSP;
      B2=new BSP;
      double tmp;
      
      //Ferr <<"decoupage de la boite selon l'axe "  <<  (int)axe<<" " ;

      // colle la limite des 2 boites-fille a la grille de voxels
      frontier=(B->maxi(axe)+B->mini(axe))/2.0;
      indi=coord(axe,frontier);
      if(verbose>3) 
	Ferr <<" [division] coupure selon l'axe "  << axe
	     <<" au debut de la boite "  << indi<<"\n" ;
      tmp=indi*taille_vox+O[axe];
      if( (frontier-tmp)>taille_vox/2.0)
	frontier=tmp+taille_vox;
      else
	frontier = tmp;
      
      //Ferr <<" frontier = "  << frontier<<"\n" ;
      B->decoupe_boite(axe, B1, B2,frontier);
      delete B;
      axe=(axe+1)%3;
      while( taille_min[axe])
	axe=(axe+1)%3;
      if(verbose>3) 
	Ferr <<"division des nouvelles boites\n" ;
      division(axe, niveau_sub_max, B1);
      division(axe,niveau_sub_max, B2);
    }
  }
  
  //Ferr <<"Voxel[division] FIN\n" ;
}//Voxel::division()

void  Voxel::creation(reel* bornemin, reel* bornemax,int* nb_vox, ListeD<Diffuseur *>& Ldiff){
  // Ferr << "Voxel[creation] DEBUT\n"; cerr.flush();
  BSP *B;
  int i, j;
  unsigned char subdiv_max[3]= {0,0,0};
  enum {axeX,axeY,axeZ};
  B=new BSP;
  if(verbose>2) 
    Ferr <<"Voxel[creation] volume englobant\n" ;
  B->volume_englobant_scene(bornemin, bornemax, Ldiff);
  if(verbose>2) 
    Ferr <<"Voxel[creation] division grille\n" ;
  division(axeZ, subdiv_max, B);
  /* Lbox.debut();
     B=Lbox.contenu();
     Ferr <<"Voxel[creation]boite.nb_diffuseur= "  <<  B->nb_diffuseur()<<"\n" ;
  */
  if(verbose){
    Ferr <<"Voxel[creation]declaration de SI\n" ;Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
    Ferr <<"\t max subdiv[0] = "  << (int)subdiv_max[0]<<"\n" ;
    Ferr <<"\t max subdiv[1] = "  <<  (int)subdiv_max[1]<<"\n" ;
    Ferr <<"\t max subdiv[2] = "  <<  (int)subdiv_max[2]<<"\n" ;
  }
  //Ferr <<"Voxel[creation] FIN\n" ;

}//Voxel::creation


int num_box=0;
void Voxel::numerotation(BSP* B){
  int i, j, k;
  int lim[3][2],ii;
  
  //Ferr <<"Voxel[numerotation] DEBUT box.num = "  << num_box<<"\n" ;

  B->baptise(num_box++);
  // Calcul des voxels contenus dans la boite B
  for(ii=0;ii<3;ii++) {
    lim[ii][0]=(int)((B->mini(ii)-O[ii])/taille_vox);
    if((B->mini(ii)- (lim[ii][0]*taille_vox+O[ii]))>taille_vox/2.0)
      lim[ii][0]++;
    lim[ii][1]=(int)((B->maxi(ii)-O[ii])/taille_vox);
    if((B->maxi(ii)- (lim[ii][1]*taille_vox+O[ii]))>taille_vox/2.0)
      lim[ii][1]++;
    /*
      Ferr <<" bmin["  << ii<<"]="  << B->mini(ii) <<" - bmax["  
      << ii<<"]="  << B->maxi(ii)<<"\n" ;
      Ferr <<" lim_inf["  << ii<<"]="  <<  lim[ii][0]<<" - lim_sup["
      << ii<<"]="  <<  lim[ii][1]<<"\n" ;
      Ferr <<" lim_min["  << ii<<"]="  <<  lim[ii][0]*taille_vox+O[ii]
      <<" - lim_max["  << ii<<"]="  <<  lim[ii][1]*taille_vox+O[ii]<<"\n" ;
    */

  }//for ii
  //Remplissage de SI avec B
  for(i=lim[0][0];i<lim[0][1];i++)
    for(j=lim[1][0];j<lim[1][1];j++)
      for(k=lim[2][0];k<lim[2][1];k++) {
	SI(i,j,k)=B;
      }//for k
  //Ferr <<"Voxel[numerotation] FIN\n" ;

}//Voxel::numerotation()

//************* Voxel:: visu()  ***********
void Voxel::visu(){
  register int i,j,k;
  BSP * B;
  
  Ferr <<"Voxel[visu()] DEBUT\n" ;
  for(i=0;i<SI.maxi()[0];i++)
    for(j=0;j<SI.maxi()[1];j++)
      for(k=0;k<SI.maxi()[2];k++){ 
	Ferr <<"SI("<<i<<","<<j<<","<<k<<") :"; 
        B=SI(i,j,k);
	//Ferr <<" Boite "  <<  B->nom()<<" contenanFerr << __FILE__ " : "<< __LINE__ << '\n' ;t "  << B->nb_diffuseur()<<" diffuseurs\n" ;
      }
  Ferr <<"Voxel[visu()] FIN" << "\n" ;
  //fflush(stdout);
}// Voxel::visu()

//-****************** Voxel::construction()  *****************
void Voxel::construction(reel* bornemin, reel* bornemax,
			 double Renv, ListeD<Diffuseur*>& Ldiff) {
  double delta[3],reste,ratio=1.0 ;
  register int i;
  FILE *grille_par;
  char grille_nom[]="grille_par";
  
  if(verbose>2)  
    Ferr <<"Voxel[construction] DEBUT\n" ;
  //calcul des dimensions de la grille en fct de Renv et ajustement des bornes
  grille_par=fopen(grille_nom,"r");
  if(grille_par==NULL){
    Ferr <<"<!> Voxel::construction() -> fichier \""  << grille_nom 
	 <<"\" pas trouve \n   => valeurs par defaut des parametres de subivision adaptative\n" ;
    SEUIL_NB_DIFF=50; //arret de la subdivision
    ratio=1; //ration entre arete de la grille et diametre de la sphere englobante
    strcpy(grille_nom,"Valeurs par defaut ");
  }
  else
    fscanf(grille_par,"%d %g",SEUIL_NB_DIFF,ratio);
  Ferr <<"Seuil dans "  << grille_nom<<" = "  << SEUIL_NB_DIFF
       <<" - ratio = "  << ratio<<"\n" ;
  O=bornemin;
  taille_vox=2.0*ratio*Renv;
  if(verbose>1) 
    Ferr <<"Taille voxel = "  << taille_vox<<"\n" ;
  for(i=0;i<3;i++){
    delta[i] = bornemax[i]-bornemin[i];
    //Ferr <<"delta["  << i<<"]= "  << delta[i]<<" \n" ;

    reste=delta[i]/taille_vox;
    nb_vox[i]=(int)reste;
    if(nb_vox[i]<reste)
      (nb_vox[i])++;
    delta[i]=(taille_vox*nb_vox[i]-delta[i])/2.0;
    Ferr <<"delta2["  << i<<"]= "  << delta[i]<<" \n" ;
    bornemax[i]+=delta[i];
    bornemin[i]-=delta[i];
    Ferr <<"xxx  bmin["  <<  i<<"]= "  << bornemin[i]<<" - bmax[" << i<<"]= "
	 << bornemax[i]<<" - nb_vox["  << i<<"]= "  << nb_vox[i]<<" \n" ;
  }
  SI.alloue(nb_vox[0],nb_vox [1],nb_vox [2]);
  SI.maj((BSP *)0 );
  // printf"Voxel[construction] creation de la grille\n)"; 
  creation(bornemin, bornemax, nb_vox, Ldiff);
  if(verbose) 
    Ferr <<"Voxel[construction] nbre de boites cre'e'es = "  << num_box<<"\n" ;
  // visu();
  Ferr <<"Voxel[construction] FIN\n" ;
}//Voxel::construction()




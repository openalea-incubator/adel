#include "bsp.h"
#include "outils.h"
// BSP

BSP::~BSP(){
  // 0cout << "destructeur de BSP\n"; cout.flush();
  Ldiff.free_liste();
}

void BSP::destruction_boite(){
  delete this;
}
//modif pour traiter un volume fixe par l'user  (cas infini) - MC 95
void BSP::volume_englobant_scene(reel *bornemin, reel *bornemax,ListeD<Diffuseur*>& Ldiff_scene){
  //modif : pour l'infini : inclure les a-cheval
  register int i,nbs;
  Diffuseur *pdiff;
  Point pt;
  // char bps[3];//Bon Pour le Service 
  for (i=0; i<3; i++) {
    position_max[i]=bornemax[i];
    position_min[i]=bornemin[i];
  }
  if (Ldiff_scene.est_vide() == VRAI){
    cout << "ERREUR - Pas de diffuseur dans la scene\n"; cout.flush();
    exit (2);
  }
  else{
    //visualise le tri pour l'infini
    Primitive* prim;
    FILE* ftri;
    ftri=fopen("trinf.can","w");
    
    nb_diffuseurs=0;
    for( Ldiff_scene.debut(); ! Ldiff_scene.finito();
	 Ldiff_scene.suivant() ){
      pdiff=Ldiff_scene.contenu();
      /* test BPS deporte dans parse_can()
	 //teste si dans le volume
	 for (i=0; i<2; i++) {
	 bps[i]=pdiff->primi().nb_in(vmin[i],vmax[i], i);
	 printf("bps[%d]=%d\n",i,(int)bps[i]);
	 }
	 nbs=pdiff->primi().nb_sommet();
	 if(fabs(bps[0])<nbs && fabs(bps[1])<nbs && fabs(bps[2])<nbs) {//Ote les exterieurs
	 if( bps[0]<=0 && bps[1]<=0 ) {//ote les a cheval sur bornes max
	 */
      Ldiff.ajoute(pdiff);
      nb_diffuseurs++;
      prim=&(pdiff->primi());

/*       fprintf(ftri," p  1 %ld %d \t",prim->name(),prim->nb_sommet()); HA */
      /* d'apres canopy_io.cpp */
      fprintf(ftri," p  1 %.0f %d \t",prim->name(),prim->nb_sommet());

      for(i=0;i<prim->nb_sommet();i++)
	fprintf(ftri," %lf %lf %lf  ",prim->sommets(i)[0],prim->sommets(i)[1],prim->sommets(i)[2]);
      fprintf(ftri," \n");
  
    }
    fclose(ftri);
  }
  if(verbose>1 ) printf("BSP[Vol englob] : nb_diff_scene= %d\n", nb_diffuseurs);
}//volume_englobant_scene(Liste)

int BSP::partage_diffuseur(Diffuseur* no_diff, BSP* B1, int axe){
  if (no_diff->primi().tout_point_inf(B1->position_max[axe], axe))
    return 1;
  else
    if (no_diff->primi().tout_point_sup(B1->position_max[axe], axe))
      return 2;
    else
      return 0;
}

inline void BSP::ajouter(Diffuseur* no_diff){
  Ldiff.ajoute(no_diff);	
  nb_diffuseurs++;
}

void BSP::copie_boite(BSP* B1){
  register int i;
  
  for (i=0; i<3; i++){
    B1->position_max[i]=position_max[i];
    B1->position_min[i]=position_min[i];
  }
}

void BSP::what_in() {

  cout<<"[] BSP : "<<this;
  if(nb_diffuseurs==0) printf("\t BSP vide\n");
  else
    for(Ldiff.debut();!Ldiff.finito();Ldiff.suivant()){
      printf("\t");
      Ldiff.contenu()->primi().qui();
      printf("\n");;
    }
  Ldiff.debut();
}//BSP::what_in()

void BSP::remplir_diffuseur(BSP* B1, BSP* B2, int axe){
  int partage;
  int i=0;
  
  for(Ldiff.debut();!Ldiff.finito();Ldiff.suivant()){
    /*  printf("BSP::remplir_diffuseur: diff %.0lf(%d) - axe=%d\n",
	Ldiff.contenu()->primi().name(), Ldiff.contenu()->num(), (int) axe); */
    partage=partage_diffuseur(Ldiff.contenu(), B1, axe);
    // cout << "partage=" << partage << endl;
    switch (partage) {
    case 0 :
      B1->ajouter(Ldiff.contenu());
      B2->ajouter(Ldiff.contenu());
      break;
    case 1 : B1->ajouter(Ldiff.contenu()); break;                 
    case 2 :  B2->ajouter(Ldiff.contenu()); break;         
    }//switch
  }//for Ldiff
}//remplir_diffuseur()

void BSP::decoupe_boite(int axe, BSP* B1, BSP* B2, double frontier){
  // adapte les bornes des filles a la grille de voxels
  
  copie_boite(B1);
  copie_boite(B2);
  B1->nb_diffuseurs=B2->nb_diffuseurs=0;
  B1->position_max[axe]=B2->position_min[axe]=frontier;
  remplir_diffuseur(B1, B2, axe);
}

bool BSP::contient(Point &P) {
  register int ii;
  bool dedans=true;
  for(ii=0;ii<3;ii++) {
    dedans &= (P[ii]<position_max[ii]) && (P[ii]>position_min[ii]);
  }
  return dedans;
}//BSP::contient()







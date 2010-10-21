/*
************************************************************
* 05/04/05 : dans la fontion parse_can, modification de la facon d'afficher 
* les erreurs d'ouverture de fichier quand la variable est nulle
* 05/04/08 : dans la fontion parse_can, nouvelle modification 
* de la facon d'afficher again 
************************************************************
*/

#include <iostream>
using namespace std;
#include <sstream>

#include "canopyL.h"
#include <cmath>
#include "outils.h"
#include <unistd.h>
// CANOPY

/************************************************************
*                      Prototypes                           *
*************************************************************/

void Canopy::gencan(char *fname){
  //genere a partir de la liste le fichier can avec indice neg pour les diff 8 et
  // en otant les non-actifs
  FILE *fic;
  Diffuseur *diff;
  Primitive *prim;
  register int i,nbs,nbd=0,nbp=0;
  int label;

  fic=fopen(fname,"w");
  for( liste_diff_scene.debut(); ! liste_diff_scene.finito();
       liste_diff_scene.suivant() ){
    diff=liste_diff_scene.contenu();
    if(diff->actif()){
      prim=&(diff->primi());
      if(diff->isreal())
	label=0;
      else
	label=8;
      nbs= prim->nb_sommet();
      fprintf(fic,"p 2 %.0f %d   %d  ",(double)prim->name()*1000,label,nbs);
      for (i=0; i<nbs; i++){
	fprintf(fic," %f %f %f ",(double)prim->sommets(i)[0],(double)prim->sommets(i)[1],(double)prim->sommets(i)[2]);
      }
       fprintf(fic,"\n");
       nbd++;
    }
    else
      nbp++;
  }
  printf("Canopy::gencan(\"%s\"): %d primitives + %d captflux = %d\n",fname,nbd,nbp,nbd+nbp);
  fclose(fic);
}//gencan()



void syntax_error(char *);

// FONCTION QUI CREE LA LISTE DES DIFFUSEURS DE LA SCENE
// FORMAT FICHIER : 1entier : libelle ; 3 coordonees (double) pour 3 sommets d'un triangle

const bool opak=true;
//lectopt() : lit les proprietes optiques (fonction locale)
Actop* lectop(ifstream &fopti, bool opac=false){
  double popt[9]={0, 0, 0, 0, 0, 0, 0, 0, 0};
  Actop *actop=NULL;
  char c;
  
  fopti>>c; //cerr<<c;
  switch(c) {
  case 'd':// diffuse (Lambert's law)
    fopti>>popt[0]; // reflectance
      if(!opac)
	fopti>>popt[1];  // transmitance
      actop = new Lambert(popt[0], popt[1]);    
    break; 
  case 's':// specular (Fresnel's law + attenation(van der Bild)
    if(!opac){
      cerr<<"speculaire pas prevu today pour non-opaques!\n";
      exit(-1);
    }
    fopti>>popt[0]; // indice de refraction n
    fopti>>popt[1]; // "hair index" , default 0.1
    actop = new Fresnel(popt[0], popt[1]);    
    break;
  case 'm':  //mixte (Lambert + specular)
    fopti>>popt[0]; // indice de refraction n
    fopti>>popt[1]; // "hair index"  
    fopti>>popt[2]; // reflectance lambert
    if(!opac)
      fopti>>popt[3]; //transmittance lambert
    actop = new Specdifu(popt[0],popt[1],popt[2],popt[3]); 
    break;
  case 'g':  //speculaire a lobe Gaussien
    if(!opac){
      cerr<<"speculaire gaussien pas prevu today pour non-opaques!\n";
      exit(-1);
    }
    fopti>>popt[0]; // indice de refraction n
    fopti>>popt[1]; // "hair index"  
    fopti>>popt[2]; // zeta : lobe spread
    actop = new Gauss(popt[0],popt[1],popt[2]); 
    break;
  case 'r':  //Ross model (Lambert + Gaussian specular)
    fopti>>popt[0]; // indice de refraction n
    fopti>>popt[1]; // "hair index"  
    fopti>>popt[2]; // zeta : lobe spread
    fopti>>popt[3]; // reflectance lambert
    if(!opac)
      fopti>>popt[4]; //transmittance lambert
    actop = new Ross(popt[0],popt[1],popt[2],popt[3],popt[4]); 
    break;
  case 'S':  //SOIL model (Lambert + Gaussian hot spot (~SOILSPEC))
    fopti>>popt[0]; // albedo de simple diffusion w
    fopti>>popt[1]; // rugosite h
    fopti>>popt[2]; // param de la fct de phase : b
    fopti>>popt[3]; //                            c
    fopti>>popt[4]; //                            bb
    fopti>>popt[5]; //                            cc
    fopti>>popt[6]; // zeta : lobe spread
    fopti>>popt[7]; // reflectance lambert
    if(!opac)
      fopti>>popt[8]; //transmittance lambert
    actop = new Sol(popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8]); 
    break;
  case 'n': //neural network's regression
    cerr<<"@!$ description des prop. opt. par reseau de neurones pas implementee\n";
    exit(-1);	 
    break;
  default  :
    syntax_error((char*)"file.opt");  
  }//switch c
  return actop;
}//lectopt()



void syntax_error(char * nomfic){
  cerr<<" Canopy_io.C: Syntax Error during parcing "<<nomfic<<"\n";
  exit(-1); 
}  

//-********************   Canopy::parse_can()    ***********************
//-**** format des libelles : espece, no plante, no tige/feuille, no triangle

//-**************** not_yet() *************************************
inline void not_yet(char * type){
  cerr<<"desole mais le type\""<<type<<"\" n'est pas encore implemente : Ligne ignoree... \n";
}// not_yet()


// char linie[500];
// inline 
char* endline(ifstream & fin){
  long int k=0;
  char car;
  char *linie=new char[500];
  //  ostrstream ligne;
  do {
    fin.get(car);
    linie[k]=car;
    k++;
  }while(car!='\n');
  linie[k]=0;
  //cerr<<"\tCanopy[parse_can] finligne : chaine = "<<linie<<"de longueur "<<k<<endl;
  return linie;
}//endline()

void Canopy::parse_can(char *ngeom,char *nopti,reel *bornemin,reel*bornemax,bool sol,char *name8){
  bool rejet=false,infty=false;
  register int i=0,nbp=0;
  Diffuseur* diff;
  ifstream fopti;
  char c, line[256];
  int nbopt=0,ii=0, No=20;
  Tabdyn<Actop*,1> tabopaque;
  Tabdyn<Actop*,2> tabtransp;
  
  ifstream fgeom; 

  if (ngeom == NULL) {
    cerr << "ERREUR - Pas de nom de maquette : ";
    cerr << "Je quitte" <<endl; 
    cerr.flush();
    exit (1);
  } else {
    if(access(ngeom,R_OK )){
      cerr << "ERREUR - Impossible d'ouvrir la maquette ";
      cerr <<ngeom<<endl;
      cerr << "Je quitte" <<endl; 
      cerr.flush();
      exit (1);
    }
  }
  // Ok
  fgeom.open(ngeom,ios::in);
    
  //   if (nopti == NULL){
  if (access(nopti,R_OK )){
    if (nopti != NULL){
      cerr << "ERREUR  - Impossible d'ouvrir le nir : ";
      cerr <<nopti<<endl ;
    } else {
      cerr << "Attention - Pas de proprietes optiques"<<endl;
    }
    cerr << " ==> Utilisation des valeur par defaut du NIR" <<endl ;
    cerr.flush();
    //exit(-1);
    // Chargement des proprietes par default pour 5 especes en NIR
    tabtransp.alloue(No,2);
    tabopaque.alloue(No+1);
    //sol
    tabopaque(0) = new Lambert(.35, 0);
    for(short ie=0;ie<No;ie++){
      tabopaque(ie+1) = new Lambert(.40, 0);
      tabtransp(ie,0) = new Lambert(.40, 0.45);
      tabtransp(ie,1) = new Lambert(.40, 0.45);
    }
  }
  else{
    // lecture des proprietes optiques (fichier '.opt')
    fopti.open(nopti,ios::in);
    do{
      fopti>>c; //cerr<<c;
      if(!fopti) break;
      switch(c) {
      case '#':
	fopti.getline(line,256);
	break;
      case 'n':
	fopti>>ii;
	cerr<<" nb especes (optiques) = "<<ii<<endl;
	tabtransp.alloue(ii,2);
	tabopaque.alloue(ii+1); 
	
	fopti.getline(line,256);		        
	break; 
      case 's':
	if(ii==0) syntax_error(nopti);
	cout<<"p.o. sol lues\n";
	//cerr<<"Canopy[parse_can]ficoptik : sol lu\n";
	tabopaque(nbopt) = lectop(fopti, opak); 
	raus(tabopaque(nbopt)==NULL,"Canopy[parse_can] allocation tabopaque impossible!");  
	nbopt++;
	fopti.getline(line,256);
	break;
      case 'e':
	if(ii==0) syntax_error(nopti);
	//cerr<<"Canopy[parse_can]ficoptik : espece no "<<nbopt<<endl;
	tabopaque(nbopt) =  lectop(fopti, opak); 
	raus(tabopaque(nbopt)==NULL,"Canopy[parse_can] allocation tabopaque impossible!");
	tabtransp(nbopt-1,0)= lectop(fopti); 
	raus(tabtransp(nbopt-1,0)==NULL,"Canopy[parse_can] allocation tabtransp_sup impossible!");
	tabtransp(nbopt-1,1)= lectop(fopti); //face inf
	raus(tabtransp(nbopt-1,1)==NULL,"Canopy[parse_can] allocation tabtransp_inf impossible!");      
	nbopt++;
	fopti.getline(line,256);
	break; 
      default  :
	syntax_error(nopti);  
      }//switch c
      // cout <<"fopti="<<!fopti<<endl; 
    } while(fopti && (nbopt<=ii));
    if(nbopt<ii)  syntax_error(nopti);
  }  
  //cerr<<"-_-_-_-_-_  Proprietes optiques chargee\n";
  
  // lecture des primitives geometriques (fichier '.can')
  //-*****************  Parser de fichier .can *************  
  Primitive* prim = NULL;
  char T,*pch;
  int nbid;
  //Tabdyn<long int, 1> tabid; :ne marche que sur l'alpha : 64bits
  
  Tabdyn<double, 1> tabid;
  double nom,espid;
  bool valid;
  reel min[3],max[3];
  register int id;
  unsigned char specie; 
  bool opak;
  espid=1000000;
  espid*=100000;
  printf("espid=%g\n",espid);
  do {
    fgeom>> T; 
    //cerr<<T;
    if(fgeom.eof()) {
      cerr <<" -_-_-_-_-_  Primitive chargé chargees\n"<<(char)7<<endl;
      break;
    }
    valid=false;
    switch(T) {
    case '#': break;
    case 'p': valid=true; break;
    case 'n': not_yet((char*)"poly avec normales"); break;
    case 'd': not_yet((char*)"disque"); break;
    case 'y': not_yet((char*)"cylindre"); break;
    case 's': not_yet((char*)"sphere"); break;
    case 'c': not_yet((char*)"cone"); break;
    default : syntax_error(ngeom);  
    }//switch T

    if(!valid) delete  endline(fgeom);
    else{
      //-** saisie des identifiants
      fgeom>>nbid;
      //cerr<<"nbid = "<<nbid<<endl;
      if(nbid>0){
	tabid.alloue(nbid);
	//cerr<<"just apres tabid.alloue(nbid)\n"; 
	for(id=0;id<nbid;id++){
	  fgeom>>tabid(id);
	  //cerr<<tabid(id)<<":";
	}
      }//if nbid >0
      if(nbid<1){
	cerr<<"Attention : nbid<1 ==> nom = 1\n";
	tabid.alloue(1);	
	// tabid(0)=100000001000;
	tabid(0)=espid+1000.;
      }//if erreur de syntaxe nbid<1
      //format label : esp*1E11 + plante*1e6+ feuille*1e3+ triangle
      // MC Avril 98 : Pb sur sun de division de tabid(0)/1e11 !!!!
      //               => Creation de la var. double espid=1e11
      // specie=(short)(tabid(0)/100000000000);
      specie=(unsigned char)(tabid(0)/espid);
      opak=((long)(tabid(0)/1000)%1000 ==0)? true : false;
      nom=tabid(0);
      //printf(" tabdid(0)=%g,nom=%.0f, specie=%d, opak=%d \n ",tabid(0),nom,(int)specie,opak?0:1);
      //-** saisie de la geometrie
      pch=endline(fgeom);   
      switch(T) {
      case 'p': 
	prim=new Polygone(pch,nom,min,max);
	break;
      default : syntax_error(ngeom);  
      }//switch T
  
      assert (prim != 0);
      if(min[0]>max[0]){ //primitive rejete
	cerr <<" *****  Primitive rejete : libelle = "<<nom<<endl;
	//rejet=true; - MC09
	delete prim;
	prim= NULL;
      }
      else{
	for (i=0; i<3; i++){
	  bornemin[i]=T_min(bornemin[i],min[i]);
	  bornemax[i]=T_max(bornemax[i],max[i]);
	}
      }//mc09
	/* ajout d'un diffuseur a la liste */
      if(opak)  
	diff=new DiffO(prim, tabopaque(specie));
      else
	diff=new DiffT(prim, tabtransp(specie-1,0),tabtransp(specie-1,1) );   
      assert (diff != 0);
      liste_diff_scene.ajoute(diff);
      nbp++;
      // MC09 }//else rejected primi
      tabid.free();
    }//else  !valid
  }while (fgeom);
  fgeom.close();
  if(rejet)
    cout <<"Canopy[parse_can] *************  Segment(s) rejete(s) *******\n";
  cout << "Canopy [parse_can] nbre de primitives = "<<nbp<<endl;
  
  if(name8==NULL)
  infty=false;
  else{
    infty=true;
    ifstream fdim(name8,ios::in);
    fdim>>bornemin[0]>>bornemin[1];
    fdim>>bornemax[0]>>bornemax[1];
    fdim.close();
    for(i=0; i<3; i++)
      cout<<"Infty: Bornemax["<<i<<"] = "<<bornemax[i]<<"  Bornemin["<<i<<"] = "<<bornemin[i]<<endl;
  }//if !infty
  
 
   bornemin[2]-=(bornemax[2]-bornemin[2])/100.0;
  tabopaque.free();
  tabtransp.free();
}//parse_can()















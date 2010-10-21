
/*******************************************************************
*                 canopy_io.C - mike 95                            *
*            classe scene de la radioxity                          *
*          Fonctions-membres qui  gerent les E/S                   *
********************************************************************/

#include <iostream>	// pour user namespace (compile mieux)
using namespace std ;

#ifdef WIN32
#include <windows.h>	// le + près possible de namespace
#ifndef MINGW
#undef min
#undef max
#endif
#endif
#include <sstream>

#include <cmath>

#include "canopy.h"
#include "outils.h"

/*
char clef_seg_in[12] ;	//  version char* de la clef numerique
HANDLE	hSharedSegIn ;	// Handle du fichier mappé
LPVOID	lpSharedSegIn ;	// pointeur LPVOID sur seg. partagé
*/
/**************************************************************************
 *       Chargement de la scene : parse_can, read_shm                     *
 **************************************************************************/
void syntax_error(char * nomfic)
 { Ferr<<" Class Canopy Erreur de syntaxe dans le fichier "<<nomfic<<"\n";
   exit(6); 
 }  
void  err_syntax(char*msg){
  Ferr <<" \n\n<!>dans canopy_io.C, pb Capteur virtuel =>" ;
  Ferr  << msg<<"" ;
  Ferr << '\n' ;
  exit(7);
  }

void  err_syntax( string msg){
  Ferr <<" \n\n<!>dans canopy_io.C, pb Capteur virtuel =>" ;
  Ferr  << msg<<"" ;
  Ferr << '\n' ;
  exit(8);
  }

bool opak=true;

//lectopt() : lit les proprietes optiques (fonction locale)
Actop* lectop(ifstream &fopti, bool opac=false){
  double popt[4]={0.0,0.0,0.0,0.0};
  Actop *actop;
  char c;
  
  fopti>>c; //Ferr <<c;
  switch(c) {
  case 'd':// diffuse (Lambert's law)
    fopti>>popt[0]; // reflectance
      if(!opac)
	fopti>>popt[1];  // transmitance
      actop = new Lambert(popt[0], popt[1]);    
    break; 
  case 's':// specular (Fresnel's law + attenation(van der Bild)
    syntax_error((char*)" .opt : Speculaire pas implemente\n");
    break;
  case 'm':  //mixte (Lambert + specular)
    if(verbose) 
	Ferr<<"*** Attention : prise en compte uniquement de la"
	<<" partie diffuse"<<"\n";
    fopti>>popt[0]; // indice de refraction n
    fopti>>popt[1]; // "hair index"  
    fopti>>popt[2]; // reflectance lambert
    if(!opac)
      fopti>>popt[3]; //transmittance lambert
    actop = new Lambert(popt[2], popt[3]); 
    break;
  case 'n': //neural network's regression
    syntax_error((char*)" .opt :  description des prop. opt. par reseau de neurones pas implementee\n");
    break;
  default  :
    syntax_error((char*)"file.opt");  
  }//switch c
  return actop;
}//lectopt()





//-********************   Canopy::parse_can()    ***********************
//-**************** not_yet() *************************************
inline void not_yet(char * type){
  Ferr<<"desole mais le type\""<<type<<"\" n'est pas encore implemente : Ligne ignoree... \n";
}// not_yet()

inline char * endline(ifstream & fin){
  register long int iKompteur=0;
  char car;
  ostringstream ligne; // was ostrstream
  do {
    fin.get(car);
    ligne.put(car);
    iKompteur++;
  }while(car!='\n');
  
  // Apprends à utiliser les strings, Herve !!! // HA
  istringstream isTmp (ligne.str()) ;  
  char *tmp1 = new char[LONG_LIGNE_CAN],
    *pline = new char[LONG_LIGNE_CAN] ;
  strcpy (pline,"");
  if ((tmp1 != NULL) && (pline != NULL)) {
    while (isTmp.good() ){
      isTmp >> tmp1 ;
      strcat (pline, " ");
      strcat (pline, tmp1);
    }
  }
  return pline;
}//endline()

long int Canopy::parse_can(char *ngeom,char *nopti,char * name8,reel *bornemin,reel*bornemax,int sol,char *nsolem,Diffuseur **&TabDiff){
  bool rejet=false;
  int i=0,j;
  long nbp=0;
  Diffuseur* diff;
  ifstream fopti(nopti,ios::in);
  char c, line[256];
  double popt[4];
  int nbopt=0,ii=0;
  Tabdyn<Actop*,1> tabopaque;
  Tabdyn<Actop*,2> tabtransp;
  Actop *testopt; 
   
  if (!fopti){
    Ferr << "ERREUR - Impossible d'ouvrir :"<<nopti<<'\n' ;//endl;
    //Ferr.flush();
    exit(9);
  }
  
  ifstream fgeom(ngeom,ios::in);
  if (!fgeom){
    Ferr << "ERREUR - Impossible d'ouvrir :"<<ngeom<<'\n' ;//endl;;
    //Ferr->flush();
    exit(10);
  }
  // lecture des proprietes optiques (fichier '.opt')
  do{
    fopti>>c;
    if(!fopti) break;
    switch(c) {
    case '#':
      fopti.getline(line,256);
      break;
    case 'n':
      fopti>>ii;
      if(verbose>1) 
	Ferr<<" nb especes (optiques) = "<<ii<<'\n' ;//endl;
      tabtransp.alloue(ii,2);
      tabopaque.alloue(ii+1); 
      
      fopti.getline(line,256);		        
      break; 
    case 's':
      if(ii==0) syntax_error(nopti);
      //Ferr <<"Canopy[parse_can]ficoptik : sol lu\n";
      tabopaque(nbopt) = lectop(fopti, opak); 
      raus(tabopaque(nbopt)==NULL,"Canopy[parse_can] allocation tabopaque impossible!");  
      nbopt++;
      fopti.getline(line,256);
      break;
    case 'e':
      if(ii==0) syntax_error(nopti);
      //Ferr <<"Canopy[parse_can]ficoptik : espece no "<<nbopt<<'\n' ;//endl;
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
    // cout <<"fopti="<<!fopti<<'\n' ;//endl; 
  } while(fopti && (nbopt<=ii));
  if(nbopt<ii)  
	syntax_error(nopti);  
  if(verbose>1) 
	Ferr<<"-_-_-_-_-_  Proprietes optiques chargees\n";
  
  //cas infini
  if(name8!=NULL) {
    ifstream finf(name8,ios::in);
    finf>>bornemin[0]>>bornemin[1];
    finf>>bornemax[0]>>bornemax[1];
    finf.close();
    delta[0]=bornemax[0]-bornemin[0];
    delta[1]=bornemax[1]-bornemin[1];
    if(verbose)  printf("parse_can() : %s - %f - %f\n",name8,bornemin[0],bornemax[1]);
    infty=true;
  }//if infty
  else
    infty=false;
  // lecture des primitives geometriques (fichier '.can')
  //-*****************  Parser de fichier .can *************  
  Primitive* prim;
  char T,*pch=NULL;
  int nbid,nl=0;
  Tabdyn<double, 1> tabid;//long ==> double
  long int idb=0;
  bool valid;
  reel min[3],max[3];
  double smax=-1,espid,nom;
  int id;
  short specie;
  char acv;
  espid=1000000;
  espid*=100000;
  //printf("espid=%g\n",espid);

  do {
    nl++;fgeom>> T; //printf("T(%d)=%c\n",nl,T);
    if(fgeom.eof()) {
      //Ferr << "fin du fichier\n";
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

    if(!valid) delete endline(fgeom);
    else{
      //idb++; printf("+ ligne %ld lue: ",idb);fflush(stdout);
      //-** saisie des identifiants
      fgeom>>nbid;
      //Ferr <<"nbid; = "<<nbid<<'\n' ;//endl; fflush(stderr);
      if(nbid>0){
	tabid.alloue(nbid);	
	for(id=0;id<nbid;id++){
	  fgeom>>tabid(id);
	  //Ferr <<tabid(id)<<"*\n";
	}
      }//if nbid >0
      if(nbid<1){
	Ferr<<"Attention : nbid<1 ==> nom = 1\n";
	tabid.alloue(1);	
	tabid(0)=espid+1000.;
      }//if erreur de syntaxe nbid<3
      //format label : esp*1E11 + plante*1e6 + feuille*1E3 + triangle
      // MC Avril 98 : Pb sur sun de division de tabid(0)/1e11 !!!!
      //               => Creation de la var. double espid=1e11
      // specie=(short)(tabid(0)/100000000000);
      specie=(short)(tabid(0)/espid);
      
      opak=((long)(tabid(0)/1000)%1000 ==0)? true : false;
      // MCjune08: nom=tabid(0)/1000;
      nom=tabid(0);
/* debug ?
      if (false){
	if(opak)
	  printf("===>  espece %d opaque %c\n",specie,7);
	else
	  printf("===>  espece %d transparente %c\n",specie,7);
      }
*/

      //-** saisie de la geometrie
      pch=endline(fgeom);   
      switch(T) {
      case 'p': prim=new Polygone(pch,nom,min,max); break;
      default : syntax_error(ngeom);  
      }//switch T
      delete pch;
      acv=0;
      assert (prim != 0);
      rejet=false;
      if(min[0]>max[0]){ //primitive rejete
	//Ferr <<" *****  Primitive rejete : libelle = "<<nom<<'\n' ;//endl;
	rejet=true;
      }
      /* Traitement des a-cheval ici et non dans BSP::volume_englobant,
	 co parcinopy a cause de visu3d.C*/
      else 
	if(infty) {
	  //printf("infty!!!!!");
	  Point G;
	  int bps[3],cpts;//Bon Pour le Service
	  rejet=true;
	  G=prim->centre();
	  if(G[0]>bornemin[0] &&G[1]>bornemin[1]&&G[0]<=bornemax[0]&&G[1]<=bornemax[1]){
	    //centroide (G) dans le cube
	    rejet=false;
	    cpts=prim->nb_sommet();
	    for (i=0; i<2; i++) {
	      bps[i]=prim->nb_in(bornemin[i],bornemax[i], i);
	      //printf("bps[%d]=%d\n",i,(int)bps[i]);
	    }
	    if(bps[0]!=0) acv++;
	    if(bps[1]!=0) acv+=2;
	  }
	  else
	    rejet=true;
	}//if infty
      if(rejet) {
	//printf( "Rejetee !\n");
	Ldiff0.ajoute(nom);
	delete prim;
      }
      else {
	//printf( " Ok\n");
	//cout<<"prim="<<prim<<'\n' ;//endl;
	//printf("prim(%ld) = %ld\n",nbp,(long) prim);
	//printf(":::> smax=%g\n",smax);
	//printf("prim(%ld) = %ld\n",nbp,(long) prim);
	smax=(prim->surface()>smax)? prim->surface() : smax;
	for (i=0; i<3; i++){
	  if(infty)
	    i=2;
	  bornemin[i]=T_min(bornemin[i],min[i]);
	  bornemax[i]=T_max(bornemax[i],max[i]);
	}
	/* ajout d'un diffuseur a la liste */
	if(opak)  
	  diff=new DiffO(prim, tabopaque(specie));
	else
	  diff=new DiffT(prim, tabtransp(specie-1,0),tabtransp(specie-1,1) );   
	assert (diff != 0);
	//cout<<"numero = "<<diff->num()<<'\n' ;//endl;
	diff->acv=acv;
	Ldiff.ajoute(diff);
	Ldiff0.ajoute(-1); //bon triangle : code label <0 - MC10
	nbp++;
      }//else rejected primi
      tabid.free();
    }//else  !valid
  }while (fgeom);
  fgeom.close();
  //  if(rejet) cout <<"Canopy[parse_can] *************  Segment(s) rejete(s) *******\n";
  if(verbose)  cout << "Canopy [parse_can] nbre de primitives  ss sol = "<<nbp<<'\n' ;//endl;
  if(verbose>1)  cout << "Canopy [parse_can] surface max primitive      = "<<smax<<'\n' ;//endl;
  
  // sol
  if(sol!=0){
    reel dx,dy,p[2];
    int nbs;
    
    nbs=(int)(sqrt((bornemax[0]-bornemin[0])*(bornemax[1]-bornemin[1])/4/smax))+1;
    printf("* surface-based number of soil triangles = %d\n", nbs*nbs*2);
    if(sol>0 && sol< nbs*nbs*2){
      nbs=(int)sqrt(sol/2.);
      printf("<!>The mesh of soil use the real threshold value (%d T <= %d)\n", 
	     sol, nbs*nbs*2 );	
    }
    dx=(bornemax[0]-bornemin[0])/(reel)nbs;
    dy=(bornemax[1]-bornemin[1])/(reel)nbs;
    p[0]=bornemin[0];
    p[1]=bornemin[1];
    for(i=0;i<nbs;i++) {
      for(j=0;j<nbs;j++) {
	{//triangle du bas
	  ostringstream ligne ;
	  ligne<<3<<" ";
	  ligne<<p[0] <<" ";
	  ligne<<p[1] <<" ";
	  ligne<<bornemin[2] <<"  ";
	  ligne<< (p[0]+dx) <<" ";
	  ligne<<p[1] <<" ";
	  ligne<<bornemin[2] <<"  ";
	  ligne<<p[0] <<" ";
	  ligne<< (p[1]+dy) <<" ";
	  ligne<<bornemin[2] <<"  ";

	  string sTmp = ligne.str() ;

	  //  ajout du sol 1 a la liste 
	  diff=new DiffO(new Polygone(sTmp,0,min,max),tabopaque(0));
	  assert (diff != 0);
	  Ldiff.ajoute(diff);
	  nbp++;
	 
	}//bloc necessaire pour use de ligne!
	{//triangle du haut
	  ostringstream ligne;
	  int iKompteur ;
	  ligne<<3<<" ";
	  ligne<< (p[0]+dx) <<" ";
	  ligne<<p[1] <<" ";
	  ligne<<bornemin[2] <<"  ";
	  ligne<< (p[0]+dx) <<" ";
	  ligne<< (p[1]+dy) <<" ";
	  ligne<<bornemin[2] <<"  ";
	  ligne<<p[0] <<" ";
	  ligne<< (p[1]+dy) <<" ";
	  ligne<<bornemin[2] <<"  ";
	  // gerer des #define si istringstream passe pas avec Bcc32
	  // memes remarques que dans le bloc precedent
	  string sTmp = ligne.str() ;
	  // iKompteur = sTmp.size();
	  istringstream isTmp (sTmp) ;
 	  //pch = new char[ iKompteur +1 ] ;
	  //isTmp >> pch ;
	  //  ajout du sol 1 a la liste 
	  // MC09 - BUG !!!! diff=new DiffO(new Polygone(pch,0,min,max),tabopaque(0));
	  //	  diff=new DiffO(new Polygone(sTmp,0,min,max),tabopaque(0)); //MC09 debug
	  //MC09 
	  // string  -> char * 	  const char *c_str() const; 
	  // ou const char *data() const;
	  pch= (char *) sTmp.c_str();
	  diff=new DiffO(new Polygone(pch,0,min,max),tabopaque(0));
	  //90CM
	  assert (diff != 0);
	  Ldiff.ajoute(diff);
	  nbp++;
	  //delete pch;
	}//bloc necessaire pour use de ligne!
	p[1]+=dy;
      }//for j
      p[0]+=dx;
      p[1]=bornemin[1];
    }//for i
    if(verbose>1) cout << "Canopy [parse_can] nbre de primitives  avec sol = "<<nbp<<'\n' ;//endl; 
  }// if sol
  for(i=0;i<3;i++) {
    vmin[i]=bmin[i]=bornemin[i];
    vmax[i]=bmax[i]=bornemax[i];
    if(verbose>1) printf("parse_can() : vmin[%d] = %g - vmax[%d] = %g \n",i,vmin[i],i,vmax[i]);
    //bornemin[i]-=(bmax[i]-bmin[i])/100.0; 
    //bornemax[i]+=(bmax[i]-bmin[i])/100.0;
  }
  tabopaque.free();
  tabtransp.free();
  nb_face=diff->idx;
  nbprim=Ldiff.card();

  //Ajout des capteurs virtuels
  if(nsolem!=NULL){
    ifstream fgeom(nsolem,ios::in);
    //if (!fgeom) 
    if (!fgeom.good()) {
      ostringstream ErrMsg;;
      ErrMsg << "ERREUR ligne "<< __LINE__;
      ErrMsg <<" - Impossible d'ouvrir le fichier \""<<nsolem<<"\"." ;
	err_syntax(ErrMsg.str());
    }
    // lecture de solem.can
    //  lecture d la premire ligne "#nb_cell"
    fgeom>> T; //Ferr <<T;
    if(T !='#')
      err_syntax((char*)"in the file .can, the first char IS NOT  '#'\n==> program aborted\n");
    fgeom>>nbcell;
    Ferr <<" "  << nbcell<<" cellules Solem positionees dans le couvert\n" ;
    if(nbcell==0) err_syntax((char*)"in the file .can, the second char should be the number of cell > 0\n==> program aborted\n");
    delete endline(fgeom);
    // Lecture des cellules define au format can
    for(i=0;i<nbcell;i++){
      valid=false;
      fgeom>> T; //Ferr <<T;
      switch(T) {
      case '#': break;
      case 'p': valid=true; break;
      case 'n': not_yet((char*)"poly avec normales"); break;
      case 'd': not_yet((char*)"disque"); break;
      case 'y': not_yet((char*)"cylindre"); break;
      case 's': not_yet((char*)"sphere"); break;
      case 'c': not_yet((char*)"cone"); break;
      default : err_syntax((char*)"Erreur de syntaxe dans le fichier can");  
      }//switch T
      if(!valid) {
	delete endline(fgeom);
	i--;
      }
      else{
	//-** saisie des identifiants
	int lab;
	fgeom>>nbid;
	for(j=0;j<nbid;j++){
	  fgeom>>id;
	  if(j==0) lab=id;
	}
	pch=endline(fgeom);
	switch(T) {
	case 'p': prim=new Polygone(pch,lab,min,max); break;
	default : err_syntax((char*)"erreur de syntaxe dans le .can");  
	}//switch T
	//  Ajout du sol 1 a la liste 
	diff=new DiffP(prim);
	assert (diff != 0);
	Ldiff.ajoute(diff);
	delete pch;
      }
    }//for i	 
    fgeom.close();
    if(verbose>1) cout<<"Canopy [parse_can] nbre de capteurs virtuels = "<<nbcell<<'\n' ;//endl;
  }// if capteur virtuel
  else {
    nbcell=0;
  }
  radim=diff->idx;
  if(verbose>1) cout<<"Canopy [parse_can] nbre de faces visibles = "<<radim<<'\n' ;//endl;
  if(verbose)   cout<<"Canopy [parse_can] nbre de primitives = "<<Ldiff.card()<<'\n' ;//endl;
 
  //mise en tableau
  FILE* fcan;
  fcan=fopen("scene.can","w");
  TabDiff= new Diffuseur *[radim];
  for(Ldiff.debut(),i=0;! Ldiff.finito();Ldiff.suivant()){
    diff=Ldiff.contenu();
    // Ecriture du .can (apres nettoyage de parse_can et ajout du sol)
    prim=&(diff->primi());
    fprintf(fcan,"p  1 %.0f %d \t",prim->name()*1000.,prim->nb_sommet());
    for(char iii=0;iii<prim->nb_sommet();iii++)
      fprintf(fcan," %lf %lf %lf  ",prim->sommets(iii)[0],prim->sommets(iii)[1],prim->sommets(iii)[2]);
    fprintf(fcan," \n");
    //Remplissage du tableau TabDiff
    TabDiff[i++]=diff;
    if(!diff->isopaque())
      TabDiff[i++]=diff;
  }// for Ldiff
  return radim;
}//parse_can()

//-********************   Canopy::read_shm()    ***********************
long int Canopy::read_shm(
	int clef,
	char *nopti,
	char *name8,
	reel *bornemin,
	reel *bornemax,
	int sol,
	char *nsolem,
	Diffuseur **&TabDiff)
{
  bool rejet=false;
  int i=0,j,it,nbp=0,Nt;
  Diffuseur* diff;
  ifstream fopti(nopti,ios::in);
  char c, line[256];
  double popt[4];
  int nbopt=0,ii=0;
  Tabdyn<Actop*,1> tabopaque;
  Tabdyn<Actop*,2> tabtransp;
  Actop *testopt; 
  Patch *Ts;

 
  if (!fopti){
    Ferr << "ERREUR - Impossible d'ouvrir :"<<nopti<<'\n' ;//endl;;
    //Ferr->flush();
    exit(11);
  }

  // initialisation du segpar 
  //BUG MC feb 2006 !!!:  Nt=clef/100;
  Nt=clef/100-1;
  clef%=100;

#ifndef WIN32
  // Unix
  int shmid;
  shmid=shmget((key_t)clef,SEGSIZE*sizeof(Patch) ,IPC_CREAT|0666);
  if(shmid==-1){//en cas de pb
    Ferr << "<!> ouverture du segment partage no. "<< clef
	 <<" impossible =>exit" << "\n" ;
    exit(12);
  }
  Ts=(Patch *)shmat(shmid,0,NULL);
#else
  // Win NT
     char clef_seg_in[12] ;	//  version char* de la clef numerique
     HANDLE	hSharedSegIn ;	// Handle du fichier mappé
     LPVOID	lpSharedSegIn ;	// pointeur LPVOID sur seg. partagé

      sprintf ( clef_seg_in, "%d", clef) ;
      assert ((hSharedSegIn = OpenFileMapping (
	FILE_MAP_ALL_ACCESS,
	FALSE,
	clef_seg_in)) != NULL ) ;

      lpSharedSegIn = MapViewOfFile (
	hSharedSegIn,
	FILE_MAP_ALL_ACCESS,
	0,
	0,
	SEGSIZE*sizeof(Patch) ) ;
      assert ( lpSharedSegIn != NULL ) ;
      Ts = (Patch *)lpSharedSegIn ;
#endif
      if(true||verbose) 
	Ferr <<"Canestra{read_shm] Clef="  << clef<<", Nt="  << Nt<<"\n" ;
  // lecture des proprietes optiques (fichier '.opt')
  do{
    fopti>>c; 
    if(!fopti) break;
    switch(c) {
    case '#':
      fopti.getline(line,256);
      break;
    case 'n':
      fopti>>ii;
      if(verbose>1) 
	Ferr<<" nb especes (optiques) = "<<ii<<'\n' ;//endl;
      tabtransp.alloue(ii,2);
      tabopaque.alloue(ii+1); 
      
      fopti.getline(line,256);		        
      break; 
    case 's':
      if(ii==0) syntax_error(nopti);
      //Ferr <<"Canopy[read_shm]ficoptik : sol lu\n";
      tabopaque(nbopt) = lectop(fopti, opak); 
      raus(tabopaque(nbopt)==NULL,"Canopy[read_shm] allocation tabopaque impossible!");  
      nbopt++;
      fopti.getline(line,256);
      break;
    case 'e':
      if(ii==0) 
	syntax_error(nopti);
      //Ferr <<"Canopy[read_shm]ficoptik : espece no "<<nbopt<<'\n' ;//endl;
      tabopaque(nbopt) =  lectop(fopti, opak); 
      raus(tabopaque(nbopt)==NULL,"Canopy[read_shm] allocation tabopaque impossible!");
      tabtransp(nbopt-1,0)= lectop(fopti); 
      raus(tabtransp(nbopt-1,0)==NULL,"Canopy[read_shm] allocation tabtransp_sup impossible!");
      tabtransp(nbopt-1,1)= lectop(fopti); //face inf
      raus(tabtransp(nbopt-1,1)==NULL,"Canopy[read_shm] allocation tabtransp_inf impossible!");      
      nbopt++;
      fopti.getline(line,256);
      break; 
    default  :
      syntax_error(nopti);  
    }//switch c
    // cout <<"fopti="<<!fopti<<'\n' ;//endl; 
  } while(fopti && (nbopt<=ii));
  if(nbopt<ii)  
	syntax_error(nopti);  
  if(verbose>2) 
	Ferr<<"-_-_-_-_-_  Proprietes optiques chargees\n";
  
  //cas infini
  if(name8!=NULL) {
    ifstream finf(name8,ios::in);
    finf>>bornemin[0]>>bornemin[1];
    finf>>bornemax[0]>>bornemax[1];
    finf.close();
    delta[0]=bornemax[0]-bornemin[0];
    delta[1]=bornemax[1]-bornemin[1];
    if(true||verbose)
      printf("INFTY read_shm() : %s =[%.3f,%.3f]x[%.3f,%.3f]\n",name8,bornemin[0],bornemax[0],bornemin[1],bornemax[1]);
    infty=true;
  }//if infty
  else
    infty=false;
  // lecture des primitives geometriques (shm clef)
  char *pch=NULL;  
  Primitive* prim;
  double nom;
  bool valid;
  reel min[3],max[3];
  double smax=-1;
  short specie;
  char acv;
  int Nrejet=0;
  Ferr <<" read_shm() nombre de triangles="  << Nt<<"\n" ;
  for(it=0;it<Nt;it++){
    //identifiants
    //format label : -i opaque, 0 sol, i transparent de l'esp i
    specie=(short)(fabs(double(Ts[it].t))); // HA cast double 11 2003
    opak=(Ts[it].t<=0)? true : false;
    nom=(double)(Ts[it].t);
    //-** saisie de la geometrie
    prim=new Polygone(Ts[it].P,nom,min,max);

    //Ferr <<"=>  it="  << it<<", prim="  << (long)prim<<"\n" ;

    acv=0;
    assert (prim != 0);
    rejet=false;
    if(min[0]>max[0]){ //primitive rejete
      Ferr <<" *****  Primitive rejete car min>max!?: libelle = "<<nom<<" - No="<<it<<'\n' ;
      rejet=true;
    }
    /* Traitement des a-cheval ici et non dans BSP::volume_englobant,
       co parcinopy a cause de visu3d.C*/
    else {
      if(infty) {
	Point G;
	char bps[3],cpts;//Bon Pour le Service
	rejet=true;
	G=prim->centre();


	if(   G[0]>bornemin[0]  && G[1]>bornemin[1]
	   && G[0]<=bornemax[0] && G[1]<=bornemax[1]){
	  //centroide (G) dans le cube
	  rejet=false;
	  cpts=prim->nb_sommet();
	  for (i=0; i<2; i++) {
	    bps[i]=prim->nb_in(bornemin[i],bornemax[i], i);
	    //printf("bps[%d]=%d\n",i,(int)bps[i]);
	  }
	  if(bps[0]!=0) acv++;
	  if(bps[1]!=0) acv+=2;
	}
	else{
	  rejet=true;
	  Ferr<<"* Prim no "<<it<<": Centroide G hors cube\n";
	  Ferr<<" bmin="<<bornemin[0]<<","<<bornemin[1]<<" - bmax="<<bornemax[0]<<","<<bornemax[1]<<"\n Centroid G =";
	  G.show();
	  Ferr <<" Primitive =";
	  (*prim)[0].show();
	  (*prim)[1].show();
	  (*prim)[2].show();
	}
      }//if infty
      if(rejet) {
	Ferr<<"* Prim no "<<it<<" ==> rejettee"<<"\n";
	delete prim;
	Nrejet++;
      }
      else {
	smax=(prim->surface()>smax)? prim->surface() : smax;
	for (i=0; i<3; i++){
	  if(infty)
	    i=2;
	  bornemin[i]=T_min(bornemin[i],min[i]);
	  bornemax[i]=T_max(bornemax[i],max[i]);
	}
	/* ajout d'un diffuseur a la liste */
	if(opak)  
	  diff=new DiffO(prim, tabopaque(specie));
	else
	  diff=new DiffT(prim, tabtransp(specie-1,0),tabtransp(specie-1,1) );   
	assert (diff != 0);
	//cout<<"numero = "<<diff->num()<<'\n' ;//endl;
	diff->acv=acv;
	Ldiff.ajoute(diff);
	nbp++;
      }//else rejected primi
    }//if valid
  }//for it

  //liberation du shm
#ifndef WIN32
  // Unix
  shmdt((void*)shmid);
#else
  // Win NT
    UnmapViewOfFile(lpSharedSegIn) ; // invalidation du ptr sur mem partagee
    // Ts = NULL ;		   // tester avant ...
    CloseHandle(hSharedSegIn) ;	   // Fermeture du fichier mappé
#endif

   if(rejet){
     char Tmsg[100];
     sprintf(Tmsg,">>>  Canopy[read_shm] ****** %d  rejected triangles ****",Nrejet);
     cout <<Tmsg<<"\n";
     Ferr <<Tmsg<<"\n";
   }
  if(verbose)  cout << "Canopy [read_shm] nbre de primitives  ss sol = "<<nbp<<'\n' ;//endl;
  if(verbose>1)  cout << "Canopy [read_shm] surface max primitive      = "<<smax<<'\n' ;//endl;
  // sol
  printf("**** sol=%d\n", sol);
  if(sol){
    reel dx,dy,p[2];
    int nbs;
    
   nbs=(int)(sqrt((bornemax[0]-bornemin[0])*(bornemax[1]-bornemin[1])/4/smax))+1;
   printf("* surface-based number of soil triangles = %d\n", nbs*nbs*2);
   if(sol>0 && sol< nbs*nbs*2){
		nbs=(int)sqrt(sol/2.);
		printf("<!>The mesh of soil use the real threshold value (%d T <= %d)\n", 
					  	 sol, nbs*nbs*2 );	
	}
	dx=(bornemax[0]-bornemin[0])/(reel)nbs;
    dy=(bornemax[1]-bornemin[1])/(reel)nbs;
	p[0]=bornemin[0];
    p[1]=bornemin[1];
    for(i=0;i<nbs;i++) {
      for(j=0;j<nbs;j++) {
	{//triangle du bas
	  ostringstream ligne;
	  ligne<<3<<" ";
	  ligne<< (p[0]) <<" ";
	  ligne<< (p[1]) <<" ";
	  ligne<< bornemin[2] <<"  ";
	  ligne<< (p[0]+dx) <<" ";
	  ligne<< (p[1]) <<" ";
	  ligne<< bornemin[2] <<"  ";
	  ligne<< (p[0]) <<" ";
	  ligne<< (p[1]+dy) <<" ";
	  ligne<< bornemin[2] <<"  ";

	  // pch=ligne.str();
	  string sTmp = ligne.str() ;
	  //  ajout du sol 1 a la liste 
	  diff=new DiffO(new Polygone(sTmp,0,min,max),tabopaque(0));
	  assert (diff != 0);
	  Ldiff.ajoute(diff);
	  nbp++;

	}//bloc necessaire pour use de ligne!
	{//triangle du haut
	  ostringstream ligne;
	  ligne<<3<<" ";
	  ligne<< (p[0]+dx) <<" ";
	  ligne<< p[1] <<" ";
	  ligne<< bornemin[2] <<"  ";
	  ligne<< (p[0]+dx) <<" ";
	  ligne<< (p[1]+dy) <<" ";
	  ligne<< bornemin[2] <<"  ";
	  ligne<< p[0] <<" ";
	  ligne<< (p[1]+dy) <<" ";
	  ligne<<bornemin[2] <<"  ";

	  string sTmp = ligne.str() ;
	  //  ajout du sol 1 a la liste 
	  diff=new DiffO(new Polygone(sTmp,0,min,max),tabopaque(0));

	  assert (diff != 0);
	  Ldiff.ajoute(diff);
	  nbp++;

	}//bloc necessaire pour use de ligne!
	p[1]+=dy;
      }//for j
      p[0]+=dx;
      p[1]=bornemin[1];
    }//for i
    if(true||verbose) cout << "Canopy [read_shm] nbre de primitives  avec sol = "<<nbp<<'\n' ;//endl; 
  }// if sol
  for(i=0;i<3;i++) {
    vmin[i]=bmin[i]=bornemin[i];
    vmax[i]=bmax[i]=bornemax[i];
    if(verbose>1) printf("read_shm() : vmin[%d] = %g - vmax[%d] = %g \n",i,vmin[i],i,vmax[i]);
    //bornemin[i]-=(bmax[i]-bmin[i])/100.0; 
    //bornemax[i]+=(bmax[i]-bmin[i])/100.0;
  }
  tabopaque.free();
  tabtransp.free();
  nb_face=diff->idx;
  nbprim=Ldiff.card();

  //Ajout des capteirs virtuels
  if(nsolem!=NULL){
    ifstream fgeom(nsolem,ios::in);
    char T;
    int nbid; 
    double id;
    if (!fgeom) {
      ostringstream ErrMsg;;
      ErrMsg <<" ERREUR ligne "<< __LINE__;
      ErrMsg <<" - Impossible d'ouvrir le fichier \""<<nsolem<<"\"." ;
      err_syntax(ErrMsg.str());
    }
    // lecture de solem.can
    //  lecture d la premire ligne "#nb_cell"
    fgeom>> T; //Ferr <<T;
    if(T !='#')
      err_syntax((char*)"in the file .can, the first char IS NOT  '#'\n==> program aborted\n");
    fgeom>>nbcell;
    Ferr <<" "  << nbcell<<" cellules Solem positionees dans le couvert\n" ;
    if(nbcell==0) 
	err_syntax((char*)"in the file .can, the second char should be the number of cell > 0\n==> program aborted\n");
    delete endline(fgeom);
    // Lecture des cellules define au format can
    for(int ic=0;ic<nbcell;ic++){  
      valid=false;
      fgeom>> T; //Ferr <<T;
      switch(T) {
      case '#': break;
      case 'p': valid=true; break;
      case 'n': not_yet((char*)"poly avec normales"); break;
      case 'd': not_yet((char*)"disque"); break;
      case 'y': not_yet((char*)"cylindre"); break;
      case 's': not_yet((char*)"sphere"); break;
      case 'c': not_yet((char*)"cone"); break;
      default : err_syntax((char*)"Erreur de syntaxe dans le fichier can");  
      }//switch T
      if(!valid) {
	delete endline(fgeom);
	i--;
      }
      else{
	//-** saisie des identifiants
	fgeom>>nbid;
	for(j=0;j<nbid;j++)fgeom>>id;
	pch=endline(fgeom);
	switch(T) {
	case 'p': prim=new Polygone(pch,-(i+1),min,max); break;
	default : err_syntax((char*)"erreur de syntaxe dans le .can");  
	}//switch T
	//  Ajout du sol 1 a la liste 
	diff=new DiffP(prim);
	assert (diff != 0);
	Ldiff.ajoute(diff);
	delete pch;
      }
    }//for i	 
    fgeom.close();
    if(verbose>1) cout<<"Canopy [read_shm] nbre de capteurs virtuels = "<<nbcell<<'\n' ;//endl;
  }// if capteur virtuel
  else {
    nbcell=0;
  }

  radim=diff->idx;
  if(verbose>1) cout<<"nbre de faces visibles = "<<radim;
  if(verbose)  {
    cout<<"\nCanopy [read_shm] nbre de primitives = "
        <<Ldiff.card()<<'\n' ;//endl;
  }

  //mise en tableau
  FILE* fcan;
  fcan=fopen("scene.can","w");
  TabDiff= new Diffuseur *[radim];
  for(Ldiff.debut(),i=0;! Ldiff.finito();Ldiff.suivant()){
    diff=Ldiff.contenu();
    // Ecriture du .can (apres nettoyage de parse_can et ajout du sol)
    prim=&(diff->primi());
    fprintf(fcan,"p  1 %.0f %d \t",prim->name()*1000.,prim->nb_sommet());
    for(char iii=0;iii<prim->nb_sommet();iii++)
      fprintf(fcan," %lf %lf %lf  ",prim->sommets(iii)[0],prim->sommets(iii)[1],prim->sommets(iii)[2]);
    fprintf(fcan," \n");
    TabDiff[i++]=diff;
    if(!diff->isopaque())
      TabDiff[i++]=diff;
  }
  fclose(fcan);
  return radim;
}//read_shm()

///////////////////////// FIN READ_SHM

#ifdef _NRJ
// OBSOLESCENT- JUIN 97
//-********************   Canopy:: xabs()    *******************
//%%%%%%%%%%%%%%%%%%%       isobary()       %%%%%%%%%%%%%%%%%%
inline Point isobary(Point&A, Point&B, Point&C){
  Point G;
  
  G[0]=(A[0]+B[0]+C[0])/3;
  G[1]=(A[1]+B[1]+C[1])/3;
  G[2]=(A[2]+B[2]+C[2])/3;
  G[3]=1;
  return G;
}//isobary()

void Canopy::xabs(char* nx3d,double *bornemin,double *bornemax,bool normee)
{ register int i,j,nb_prim=0,nb_seg=0,nb_face=0,color;
   double M,Emax=-1.0,Emin=9.9e12;
   int nbcol=203;
   Point G;
   Diffuseur *pdiff; 
   ofstream fout(nx3d,ios::out);
   if (!fout)
     { Ferr << "ERREUR - Impossible d'ouvrir :"<<nx3d<<'\n' ;//endl;
       exit(13);
     }

// Ecriture du fichier nx3d au format X3D
   for(M=-99999999.9,j=0; j<3; j++)
     M=T_max(M,bornemax[j]-bornemin[j]);
   M=3000.0/M;

      // Entete
      fout<<"#     "<<nx3d<<'\n' ;//endl;
      fout<<"# genere par raddir - visu des Eabs d'un couvert - MC94\n";
      nb_prim=Ldiff.card();
      // Couleurs
      fout<<"# number of colors used in object\n";
      fout<<nbcol<<'\n' ;//endl;
      fout<<"# colors used in object (color number red green blue)\n";
      for(i=0;i<nbcol-3;)               // degrade de verts
        fout<<i++<<"  0 "<<i<<" 0\n";
      fout<<nbcol-3<<"  0   0   255\n";      // Bleu
      fout<<nbcol-2<<"  255 255 255\n";  // Blanc
      fout<<nbcol-1<<"  255 0   0\n";      // Rouge
      // Points
      fout<<"# number of points used in object\n";
      fout<<nb_prim*(normee?5:3) +4<<'\n' ;//endl;

      for(i=0,Ldiff.debut();!Ldiff.finito();Ldiff.suivant())
        { pdiff=Ldiff.contenu();
          Emin=(Emin<pdiff->nrj())? Emin : pdiff->nrj();
          Emax=(Emax>pdiff->nrj())? Emax : pdiff->nrj();
//          cout<<"Canopy[xabs]  Eabs = "<<  pdiff->nrj()<<'\n' ;//endl;
          Primitive &prim=pdiff->primi();
//          prim->qui();
          for(j=0;j<3;j++,i++)
           { fout<<i<<"\t"<<M*(prim.sommets(j)[0]-bornemin[0])-1000.0;
             fout<<"  "<<M*(prim.sommets(j)[1]-bornemin[1])-1000.0;
             fout<<"  "<<M*(prim.sommets(j)[2]-bornemin[2])-1000.0<<'\n' ;//endl;
           }
          if(normee)
	   { G=isobary(prim.sommets(0),prim.sommets(1),prim.sommets(2));
             fout<<i++<<"\t"<<M*(G[0]-bornemin[0])-1000.0;
             fout<<"  "<<M*(G[1]-bornemin[1])-1000.0;
             fout<<"  "<<M*(G[2]-bornemin[2])-1000.0<<'\n' ;//endl;
             G+=prim.normal()*(50.0/M);
             fout<<i++<<"\t"<<M*(G[0]-bornemin[0])-1000.0;
             fout<<"  "<<M*(G[1]-bornemin[1])-1000.0;
             fout<<"  "<<M*(G[2]-bornemin[2])-1000.0<<'\n' ;//endl;
	   }//if normee 
        }
      fout<<i++<<" "<<-1001<<"  "<<-1001<<"  "<<-1001<<"\n";
      fout<<i++<<" "<<1001<<"  "<<-1001<<"  "<<-1001 <<"\n";
      fout<<i++<<" "<<-1001<<"  "<<1001<<"  "<<-1001 <<"\n";
      fout<<i++<<" "<<-1001<<"  "<<-1001<<"  "<<1001 <<"\n";

      // Segments
      int inc_seg=(normee)?4:3;
      nb_seg=nb_prim*inc_seg;
      fout<<"# number segments used in object\n";
      fout<<nb_seg+3<<'\n' ;//endl;
      for(i=0;i<nb_seg;i+=inc_seg) 
        if(normee)
	  { for(j=0;j<3;j++)
              fout<<i+j<<" "<<  nbcol-2<<"  "<< i+((i/4))+j<<"   "<<i+((i/4))+((j+1)%3)<<'\n' ;//endl;
	      fout<<i+3<<" "<< nbcol-1<<"  "<< i+((i/4))+3<<"   "<<i+((i/4))+4<<'\n' ;//endl;
	  }
         else
          for(j=0;j<3;j++)
            fout<<i+j<<" "<<  nbcol-2<<"  "<< i+j<<"   "<<i+((j+1)%3)<<'\n' ;//endl;
      j=nb_prim*(normee?5:3);             
      fout<<i  <<"  "<< nbcol-3<<"  "<<j <<"  "<<j+1<<'\n' ;//endl;
      fout<<i+1<<"  "<< nbcol-3<<"  "<<j <<"  "<<j+2<<'\n' ;//endl;
      fout<<i+2<<"  "<< nbcol-3<<"  "<<j <<"  "<<j+3<<'\n' ;//endl;  

      // Faces
      cout<<"Canopy[xabs] Emin = "<<Emin<<" - Emax = "<<Emax<<'\n' ;//endl;
      fout<<"# number of polygons used in object\n";
      fout<<nb_prim<<'\n' ;//endl;
      for(i=0,j=0,Ldiff.debut();!Ldiff.finito();Ldiff.suivant(),i++)
       { pdiff=Ldiff.contenu();
         color=(int) ((double)(nbcol-4)*(pdiff->nrj()-Emin)/(Emax-Emin));
         fout<<"# "<<pdiff->primi().name()<<'\n' ;//endl;
         fout<<i<<"  "<<color<<" 3  "<<j++<<"\t"<<j++<<"\t"<<j++<<'\n' ;//endl;
         if(normee) j++;
       }
       fout<<"#*********************************************\n";
       fout<<"#$$$$$ Emin = "<<Emin<<" - Emax = "<<Emax<<"$$$$$$\n";
       fout<<"#*********************************************\n";
     fout.close();
 
 }//Canopy::xabs()

//-********************   Canopy:: xrad()    *******************

void Canopy::xrad(char* nvar,double *bornemin,double *bornemax,bool normee)
 { register int i,j,nb_prim=0,nb_seg=0,nb_face=0,color;
   double M,Emax=-1.0,Emin=9.9e12;
   int nbcol=203;
   Point G;
   Diffuseur *pdiff; 
   ofstream fout(nvar,ios::out);
   if (!fout)
     { Ferr << "ERREUR - Impossible d'ouvrir :"<<nvar<<"\n";//endl
       exit(14);
     }

// Ecriture du fichier nvar au format X3D
   for(M=-99999999.9,j=0; j<3; j++)
     M=T_max(M,bornemax[j]-bornemin[j]);
   M=3000.0/M;

      // Entete
      fout<<"#     "<<nvar<<'\n' ;//endl;
      fout<<"# genere par raddir - visu des radiosite d'un couvert - MC95\n";
      nb_prim=Ldiff.card();
      // Couleurs
      fout<<"# number of colors used in object\n";
      fout<<nbcol<<'\n' ;//endl;
      fout<<"# colors used in object (color number red green blue)\n";
      for(i=0;i<nbcol-3;)               // degrade de verts
        fout<<i++<<"  0 "<<i+5<<" 0\n";
      fout<<nbcol-3<<"  0   0   255\n";      // Bleu
      fout<<nbcol-2<<"  255 255 255\n";  // Blanc
      fout<<nbcol-1<<"  255 0   0\n";      // Rouge
      // Points
      fout<<"# number of points used in object\n";
      fout<<nb_prim*(normee?5:3) +4<<'\n' ;//endl;

      for(i=0,Ldiff.debut();!Ldiff.finito();Ldiff.suivant())
        { pdiff=Ldiff.contenu();
          Emin=(Emin<pdiff->Radmax())? Emin : pdiff->Radmax();
          Emax=(Emax>pdiff->Radmax())? Emax : pdiff->Radmax();
//          cout<<"Canopy[xrad] V(Eabs) = "<<  pdiff->Radmax()<<'\n' ;//endl;
          Primitive &prim=pdiff->primi();
//          prim->qui();
          for(j=0;j<3;j++,i++)
           { fout<<i<<"\t"<<M*(prim.sommets(j)[0]-bornemin[0])-1000.0;
             fout<<"  "<<M*(prim.sommets(j)[1]-bornemin[1])-1000.0;
             fout<<"  "<<M*(prim.sommets(j)[2]-bornemin[2])-1000.0<<'\n' ;//endl;
           }
          if(normee)
	   { G=isobary(prim.sommets(0),prim.sommets(1),prim.sommets(2));
             fout<<i++<<"\t"<<M*(G[0]-bornemin[0])-1000.0;
             fout<<"  "<<M*(G[1]-bornemin[1])-1000.0;
             fout<<"  "<<M*(G[2]-bornemin[2])-1000.0<<'\n' ;//endl;
             G+=prim.normal()*(50.0/M);
             fout<<i++<<"\t"<<M*(G[0]-bornemin[0])-1000.0;
             fout<<"  "<<M*(G[1]-bornemin[1])-1000.0;
             fout<<"  "<<M*(G[2]-bornemin[2])-1000.0<<'\n' ;//endl;
	   }//if normee 
        }
      fout<<i++<<" "<<-1001<<"  "<<-1001<<"  "<<-1001<<"\n";
      fout<<i++<<" "<<1001<<"  "<<-1001<<"  "<<-1001 <<"\n";
      fout<<i++<<" "<<-1001<<"  "<<1001<<"  "<<-1001 <<"\n";
      fout<<i++<<" "<<-1001<<"  "<<-1001<<"  "<<1001 <<"\n";

 Emin=0.0;
      // Segments
      int inc_seg=(normee)?4:3;
      nb_seg=nb_prim*inc_seg;
      fout<<"# number segments used in object\n";
      fout<<nb_seg+3<<'\n' ;//endl;
      for(i=0;i<nb_seg;i+=inc_seg) 
        if(normee)
	  { for(j=0;j<3;j++)
              fout<<i+j<<" "<<  nbcol-2<<"  "<< i+((i/4))+j<<"   "<<i+((i/4))+((j+1)%3)<<'\n' ;//endl;
	      fout<<i+3<<" "<< nbcol-1<<"  "<< i+((i/4))+3<<"   "<<i+((i/4))+4<<'\n' ;//endl;
	  }
         else
          for(j=0;j<3;j++)
            fout<<i+j<<" "<<  nbcol-2<<"  "<< i+j<<"   "<<i+((j+1)%3)<<'\n' ;//endl;
      j=nb_prim*(normee?5:3);             
      fout<<i  <<"  "<< nbcol-3<<"  "<<j <<"  "<<j+1<<'\n' ;//endl;
      fout<<i+1<<"  "<< nbcol-3<<"  "<<j <<"  "<<j+2<<'\n' ;//endl;
      fout<<i+2<<"  "<< nbcol-3<<"  "<<j <<"  "<<j+3<<'\n' ;//endl;  

      // Faces
      cout<<"Canopy[xrad] Emin = "<<Emin<<" - Emax = "<<Emax<<'\n' ;//endl;
      fout<<"# number of polygons used in object\n";
      fout<<nb_prim<<'\n' ;//endl;
      for(i=0,j=0,Ldiff.debut();!Ldiff.finito();Ldiff.suivant(),i++)
       { pdiff=Ldiff.contenu();
         color=(int) ((double)(nbcol-4)*(pdiff->Radmax()-Emin)/(Emax-Emin));
         fout<<"# "<<pdiff->primi().name()<<'\n' ;//endl;
         fout<<i<<"  "<<color<<" 3  "<<j++<<"\t"<<j++<<"\t"<<j++<<'\n' ;//endl;
         if(normee) j++;
       }
       fout<<"#*********************************************\n";
       fout<<"#$$$$$ Emin = "<<Emin<<" - Emax = "<<Emax<<"$$$$$$\n";
       fout<<"#*********************************************\n";
     fout.close();
 
 }//Canopy::xrad()

#endif














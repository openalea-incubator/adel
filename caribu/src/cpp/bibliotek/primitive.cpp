
#include <primitive.h>
#include <sstream>

#include <Mmath.h>

// PRIMITIVE

Primitive :: Primitive(Point p)
{
	pt_intersection=p;
}

// POLYGONE

void Polygone::init(Liste<Point>& liste_sommet,double name=0){
  register int i=1;
  
  nom=name;
  for( liste_sommet.debut();
      !liste_sommet.est_fin();
       liste_sommet.suivant() )
    i++;        
  
  nb_sommets=i;
  sommet=new Point[i];
  assert(sommet);
  liste_sommet.debut();
//        qui();
  for (i=0; i<nb_sommets; i++){
    sommet[i]=liste_sommet.contenu();
//                cout<<" sommet " <<i<<" "; sommet[i].show();
    liste_sommet.suivant();
  }
  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
}//init

Polygone::Polygone( Polygone& frere):Primitive(){
  register int i;
  nb_sommets=i=frere.nb_sommet();
  printf("Polygone::Polygone(Polygone&) nb_sommets=%d\n",i);fflush(stdout);
  sommet=new Point[i];
  assert(sommet);
  for (i=0; i<nb_sommets; i++)
    sommet[i]=frere.sommets(i);
  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
}//Cstructeur par copie
Polygone::Polygone(Polygone* frere){
  register int i;
  
  nb_sommets=i=frere->nb_sommet();
  sommet=new Point[i];
  assert(sommet);
  for (i=0; i<nb_sommets; i++)
    sommet[i]=frere->sommets(i);
  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
}//Cstructeur par adresse

///////////////////////////////// constructeur par nom (char*)
Polygone::Polygone(char* line,double name,
		reel*mini,reel*maxi,
		bool valid){
  init(line,name,mini,maxi,valid);
}
// constructeur par nom (char*)

///////////////////////////////// constructeur par nom (string)
Polygone::Polygone(string ligne,double name,
		   reel*mini,reel*maxi,
		   bool valid){
  //Ferr <<__FILE__<<" : "<<__LINE__<<"\n" ;
  istringstream isTmp (ligne);

  // Note: apprendre à vider TOUTE une istringstream d'1 coup

  char *tmp1 = new char[LONG_LIGNE_CAN],
    *tmp2 = new char[LONG_LIGNE_CAN] ;
  strcpy (tmp2,"");
  if ((tmp1 != NULL) && (tmp2 != NULL)) {
    while (isTmp.good() ){
      isTmp >> tmp1 ;
      strcat (tmp2, " ");
      strcat (tmp2, tmp1);
    }
    init(tmp2,name,mini,maxi,valid);
  } else {
    Ferr <<"Plus de memoire" << "\n" ;
    exit (23);
  }
  if (tmp1 != NULL) delete (tmp1);
  if (tmp2 != NULL) delete (tmp2);
}
// constructeur par nom string 

//-*******************   Polygone::init(char *,...)   ***********************
void Polygone::init(char* line,double name,
		reel*mini,reel*maxi,
		bool valid){
  //istrstream parag(line);
  //Ferr <<__FILE__<<" : "<<__LINE__<<"\n" ;
  int i,j,ii;
  Point P;
  double x;
  char *str,*str1,str2[5],str3[125],c;
  bool bb=!((maxi==NULL)||(mini==NULL));

  //Ferr << "line = " << &line<< " ; len(line) = " <<strlen(line) << "\n" ;
  //printf(" => Py.init(\"%s\",%.0g,min,max,false);\n",line,name);
  //Ferr <<line << "\n" ;

  if(bb)
    for(j=0;j<3;j++){
      mini[j]=99999999.0;
      maxi[j]=-99999999.0;
    }    
  nom=name;
  strcpy(str3,line);
  sscanf(str3,"%d",&ii);
  sprintf(str2,"%d",ii);
  str=strstr(str3,str2);//strtok(str3,str2);
  str++;
  nb_sommets=i=ii;

  // Ferr <<"init: " <<ii<<" sommets"<<"\n";
  // Ferr << __FILE__<<" : "<< __LINE__ <<"\n" ;

  sommet=new Point[i];
  assert(sommet);

  // Ferr << __FILE__<<" : "<< __LINE__ <<"\n" ;
  // qui();
 
  for (i=0; i<nb_sommets; i++){
    // Ferr << __FILE__<<" : "<< __LINE__ <<"\n" ;
    for(j=0;j<3;j++){
      // Ferr << __FILE__<<" : "<< __LINE__ <<"\n" ;
      for(str1=str;*str1== ' '|| *str1=='\t';str1++);
      str=str1;
      for(str1=str;*str1!= ' '&& *str1!='\t'&& *str1!='\n';str1++);
      *str1=0;
      sscanf(str,"%lf",&x);
      str=str1+1;
      
      //printf("=> sommet(%d,%d)=%.3g <= str=%s, str1=%s\n",
      //     i,j,x,str,str1);
      //fflush(stderr);
      sommet[i][j]=x;
      // fprintf(stderr,"<= sommet(%d,%d)=%g,   x=%f,str=%s\n",
      //      i,j,sommet[i][j],x,str);
      // fflush(stderr);
      
      if(bb){
	mini[j]=T_min((double)x,(double)mini[j]);
	maxi[j]=T_max((double)x,(double)maxi[j]);
      }
    }
  }
  // rejet test des 3 1st point
  //Ferr << "valid = "<<valid << "\n";Ferr.flush();
  if(valid){
    //printf("Polygone::init() Test de validite \n");
    if(nb_sommets>=3) {
      Vecteur u(sommet[0],sommet[1]),v(sommet[0],sommet[2]);
 
      // Ferr << __FILE__<<" : "<< __LINE__ <<"\n" ;
      x=u.norme();

      if(x==0.) {
	Ferr <<"\tPolygone::init(char*) =>"
	  " Comme les sommets 0 et 1 sont egaux, " ;
	  Ferr <<'\n' ;
	// Ferr <<"S0: "; sommet[0].show() ;
	// Ferr <<"S1: "; sommet[1].show() ;
	Ferr <<'\n' ;
	if(bb){
	  mini[0] =1;
	  maxi[0]=0;
	}
	return;
      }
      u/=x;
      x=v.norme();
      if(x==0.) {
	Ferr <<"\tPolygone::init(char*) =>"
	  " Comme les sommets 0 et 2 sont egaux, \n" ;
	if(bb){
	  mini[0] =1;
	  maxi[0]=0;
	}
	return;
      }
      v/=x;
      if(sommet[1]==sommet[2]){
	Ferr <<"\tPolygone::init(char*) =>"
	  " Comme les sommets 1 et 2 sont egaux, \n" ;
	if(bb){
	  mini[0] =1;
	  maxi[0]=0;
	}
	return;
      }
      if(fabs(u.prod_scalaire(v))>0.9999999999 ){
	Ferr <<"\tPolygone::init(char*) => "
	  "Comme les sommets 0, 1 et  2  sont aligne's, \n" ;
	if(bb){
	  mini[0] =1;
	  maxi[0]=0;
	}
	return; 
      }
      v=u.prod_vectoriel(v);
      x=v.norme();
      if(x<=0){
Ferr <<"\tPolygone::init(char*) => Norme de (POP1)^(POP2 Nulle \n" ;
	if(bb){
	  mini[0] =1;
	  maxi[0]=0;
	}
	return;
      }
    }
  }// if valid
  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
}//Polygone::init(char*) 

Polygone::Polygone(float(*T)[3],double name,reel*mini,reel*maxi){
  register int i,j;
  int ii;
  Point P;
  double x;

  for(j=0;j<3;j++){
    mini[j]=99999999.0;
    maxi[j]=-99999999.0;
  }    
  nom=name;  
  nb_sommets=i=3;
  sommet=new Point[i];
  assert(sommet);
  // qui();
  for (i=0; i<nb_sommets; i++){
    for(j=0;j<3;j++){
      P[j]=T[i][j];
      mini[j]=T_min(P[j],mini[j]);
      maxi[j]=T_max(P[j],maxi[j]);
    }    
    sommet[i]=P;
    //Ferr<<" sommet " <<i<<" "; sommet[i].show();
  }
  // rejet test des 3 1st point
  if(nb_sommets>=3) {
    Vecteur u(sommet[0],sommet[1]),v(sommet[0],sommet[2]);
    x=u.norme();
    if(x==0.) {
Ferr <<"\tPolygone(reel **) => Comme les sommets 0 et 1 sont egaux, \n" ;
      mini[0] =1;
      maxi[0]=0;
      return;
    }
    u/=x;
    x=v.norme();
    if(x==0.) {
Ferr <<"\tPolygone(reel **) => Comme les sommets 0 et 2 sont egaux, \n" ;
      mini[0] =1;
      maxi[0]=0;
      return;
    }
    v/=x;
    if(sommet[1]==sommet[2]){
Ferr <<"\tPolygone(reel **) => Comme les sommets 1 et 2 sont egaux, \n" ;
      mini[0] =1;
      maxi[0]=0;
      return;
    }
    if(fabs(u.prod_scalaire(v))>0.9999999999 ){
Ferr <<"\tPolygone(reel **) => Comme les sommets 0, 1 et  2  sont aligne's, \n" ;
      mini[0] =1;
      maxi[0]=0;
      return; 
    }
    v=u.prod_vectoriel(v);
    x=v.norme();
    if(x<=0){
Ferr <<"\tPolygone(reel **) => Norme de (POP1)^(POP2 Nulle \n" ;
      mini[0] =1;
      maxi[0]=0;
      return;
    }
  }
  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
}//Polygone(reel **)

void Polygone::show(const char* msg,ostream &out){
  out << msg<<"-Polygone :"; qui();
  for (register int i=0; i<nb_sommets; i++){
    out << "\t" <<sommet[i][0]<<" "<<sommet[i][1]<<" "<<sommet[i][2]<<endl;
  }
}//show()

Polygone::~Polygone(){
	free();
}

inline void Polygone::free(){
  delete [] sommet;
}
/*
Polygone& Polygone::operator = (Polygone& p)
{       nom=p.nom;
	nb_sommets=p.nb_sommets;
	sommet=p.sommet;
	normale=p.normale;
	cst_equ_plan=p.cst_equ_plan;

	return *this;
}
*/

Point& Polygone::operator [] (int i)
{
	assert (i>=0 && i<nb_sommets);
	return sommet[i];
}

int Polygone::tout_point_inf(reel& position, const int& axe){
  int i, compris=1;

  // cout << "position=" << position << endl;
  // cout << "axe=" << axe << endl;
  // cout << "nb_sommets=" << nb_sommets << endl;
  for (i=0; i<nb_sommets; i++){
    // cout << "sommet[i][axe]=" << sommet[i][axe] << endl;
    if (sommet[i][axe]>position)
      compris=0;
  }
  
  // cout << "compris_inf=" << compris << endl;
  return compris;
}

int Polygone::tout_point_sup(reel& position, const int& axe){
  int i, compris=1;
  
  // cout << "position=" << position << endl;
  // cout << "axe=" << axe << endl;
  for (i=0; i<nb_sommets; i++){
    // cout << "sommet[i][axe]=" << sommet[i][axe] << endl;
    if (sommet[i][axe]<position)
      compris=0;
  }
  // cout << "compris_sup=" << compris << endl;
  return compris;
}

int Polygone::nb_in(reel& amin,reel& amax, const int& axe) {
  register int nbp=0,i;
  bool info=false;
  for (i=0; i<nb_sommets; i++) {
    if(sommet[i][axe]>amax) {
      if(nbp<0) {
	Ferr<<" (!) Pb Polygone::nb_in : primitive no. "<<nom<<" plus grande que les bornes \n";
	info=true;
      }	  
      nbp++;
    }
    else
      if(sommet[i][axe]<amin) {
	if(nbp>0) {
	  Ferr<<" (!) Pb Polygone::nb_in : primitive No. "<<nom<<" plus grande que les bornes\n!";
	info=true;
	}	  
	nbp--;
      }
  }//for sommets
  if(info) {
      Ferr<<"--> Axe "<<axe<<" - bornes = ("<<amin<<", "<<amax<<")\n";
      for (i=0; i<nb_sommets; i++) {
	Ferr<<"\tP["<<i<<"] = "<<sommet[i][axe]<< '\n' ;
      }
  }
  return nbp;  
}//Polygone::nb_in()

void Polygone::calcul_normale_cst_equ(Point& p1, Point& p2, Point& p3){
  Vecteur v1, v2,v0(0,0,0);
  
  /* calcul de la normale au plan */
  v1.formation_vecteur(p1, p2);
  v2.formation_vecteur(p1, p3);
  //v1.show();
  //v2.show();
  normale=v1.prod_vectoriel(v2);
  if(normale==v0){
    Ferr<<"Polygone::calcul_normale_cst_equ : normale nulle ->";
    qui();
    exit(24);
  }
  normale.normalise();
  //        normale.show();
  /* calcul de la constant de l'equation du plan P*N-D=0 */
  Vecteur OA;
  Point O(0.0, 0.0, 0.0);
  OA.formation_vecteur(O,p1);
  cst_equ_plan=OA.prod_scalaire(normale);
  //Calcul de l'isobarycentre et de la + grde distance de isob aux sommets (aprt_sphere)
  register int i,j;
  float d2;
  isob[0]=isob[1]=isob[2]=0.;
  for(i=0;i<nb_sommets;i++)
    for(j=0;j<3;j++){
      isob[j]+=sommet[i][j];
    }
  isob=isob/(double) nb_sommets;
  d2=-1;
  for(i=0;i<nb_sommets;i++){
    dGS=isob.dist2(sommet[i]);
    if(dGS>d2){
      d2=dGS;
    }
  }
  dGS=sqrt(d2);
}

Vecteur Polygone::azi(){
// renvoie l'azimuth zero du polygone , qui correspond au 1er segment
  Vecteur u(sommet[0],sommet[1]);
  
  return u;
}//Polygone::azi()

// TRIANGLE
Triangle::Triangle(Liste<Point>& liste_sommet,double name=0) : Polygone(liste_sommet,name)
{
  if (nb_sommets != 3)
    { Ferr << "ERREUR - nombre de sommets incoherent pour un triangle\n"; 
      Ferr << "nombre de sommets=" << nb_sommets << '\n';
      exit(25);
    }
//  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
}
Triangle::Triangle(Point &A,Point &B,Point &C,double name=0)
{ sommet = new Point[3];
  sommet[0]=A;
  sommet[1]=B;
  sommet[2]=C;
  nom=name;
  nb_sommets=3;
  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
}


inline bool signe(const double& x) 
 { return(x>=0.0);
 } 
bool Polygone::dedans(const Point &I){
// teste la position de  I / au polygone poly
   int i,j,k,inc;
   Vecteur u,v;
   bool sign;

   i=0;
   if (fabs(normale[1]) > fabs(normale[0])) i=1;
   if (fabs(normale[2]) > fabs(normale[1])) i=2;   
   j=(i+1)%3;
   k=(j+1)%3;
 
   u.formation_vecteur(I,sommet[0]);
   v.formation_vecteur(I,sommet[1]);
   sign=signe(u[j]*v[k]-(u[k]*v[j]));
   for(inc=2;inc<=nb_sommets;inc++) {
     u=v;
     v.formation_vecteur(I,sommet[inc%nb_sommets]);
     if(signe(u[j]*v[k]-(u[k]*v[j]))!=sign) return false;
    }//for
   //   cout<<"Polygone|[dedans] vrai: ";qui();
   return true;   
}//dedans

double Polygone::intersect(Param_Inter& parag,Point *I) {
  // algo de SNYDER
  /* Calcul de l'intersection entre la droite definie par le segment et
     le plan defini par le polygone.

     Le plan est defini par un point A et un vecteur normal N. Un point P
     appartient alors au plan si (P-A)*n = 0. ou * = prod scal
     D'ou l'equation du plan : P*n - D = 0 ou D=A*n

     La droite est definie par: P = O + t*u ou 0 est un point de la droite
     et u, un vecteur directeur

     Intersection :  ti= (D-O*n)/(n*u)
     I = O + ti*u
     */

  double scal=normale.prod_scalaire(parag.direct());
  
  //  cout<<" Polygone[intersect] DEBUT - scal = "<<scal ; qui();
  if (fabs(scal) < 1e-6)     /* scal== 0 : droite // plan */
    return 2e9;
  double t;
  Vecteur O;
  Point OO;
  O.formation_vecteur(OO,parag.origin());
  t = ( cst_equ_plan - normale.prod_scalaire(O))/ scal;
  //cout<<"Polygone [intersect] t="<<t;qui();

  if (t <= 0.) // si I avant O
    return 3e9;
  *I=parag.origin()+ parag.direct()*t ;
  //cout<<"Polygone[intersect] inter I : "<<(*I)[0]<<"  "<<(*I)[1]<<"  "<<(*I)[2]; qui();
  // Le point appartient-il au polygone ?
  if (dedans(*I)){
    //cout<<"Polygone[intersect] dedans I : "<<(*I)[0]<<"  "<<(*I)[1]<<"  "<<(*I)[2];qui();
    return t;
  }
  else
    return 4e9;

}//intersect()

//-***************************   Polygone::surface()   ************************
double Polygone::surface(){
  double x;  
  Vecteur u,v, CB,CD;
  // Ferr<<"Polygone::surface()";qui(); //show();
  switch(nb_sommets){
  case 1: //polygone de faussaire
  case 2:
    Ferr<<"Polygone[surface] Irrtum nb_sommet = "<<nb_sommets<< '\n';
    exit(26);
    break;
  case 3: //triangle
    u.formation_vecteur(sommet[0],sommet[1]);
    v.formation_vecteur(sommet[0],sommet[2]);
    u.normalise();
    v.normalise();
    x=Macos(u.prod_scalaire(v));
    x=sommet[0].dist(sommet[1])*sommet[0].dist(sommet[2])*sin(x)/2.0;
    raus(x<0.0,"Polygone[surface] surface < 0.0 !!??"); // warning trigraphe (??.)
    break;
  case 4:
    //rectangle
    u.formation_vecteur(sommet[0],sommet[1]);
    v.formation_vecteur(sommet[0],sommet[3]);
    CB.formation_vecteur(sommet[2],sommet[1]);
    CD.formation_vecteur(sommet[2],sommet[3]);
    if(u.prod_scalaire(v)==CB.prod_scalaire(CD)){
      x=(u.prod_vectoriel(v)).norme();
      //Ferr<<"Polygone[surface]rectangle\n";
      break;
    }
  default:
    Point G=centre();
    x=0.0;
    for(register int i=0;i<nb_sommets;){
      u.formation_vecteur(G,sommet[i++]);
      v.formation_vecteur(G,sommet[(i<nb_sommets)?i:0]);
      x+=(u.prod_vectoriel(v)).norme();
      //Ferr<<"Polygone[surface]poly qcque\n";
    }
    break;
  }//switch
  return(x);   
 }// Polygone::surface() 


//-***************************   Triangle::surface()   ************************
double Triangle::surface(){
  double x;
  Vecteur AB(sommet[0],sommet[1]), AC(sommet[0],sommet[2]);   

  AB.normalise();
  AC.normalise();
  x=Macos(AB.prod_scalaire(AC));
  //   cout <<"@@@@@@@@@@ M_PI_2 = "<< M_PI_2<<" - alpha = "<<x<<endl;
  x=sommet[0].dist(sommet[1])*sommet[0].dist(sommet[2])*sin(x)/2.0;
  raus(x<0.0,"Triangle[surface] surface < 0.0 !!??");
  return(x);   
}// Triangle::surface()

//-***************************  Polygone::distance_point() ************************
#ifdef _TEST
extern char verbeux;
#else
static char verbeux=0;
#endif
//static   long int tpsD,tpsF,tps1=0,tps2=0; 
static   time_t tpsD,tpsF,tps1=0,tps2=0; 


double Polygone::distance2_point(Point &C) {
  /* Calcule le carre de la distance d'un point a un polygone convexe T en 2 etapes :
     - I, intesection du plan de T avec la droite (C,n), où n est la normale de T
     (cf. intersect())
     - Determination de I', point minimisant la distance (T,I) dans le plan  
     */
  double t,scal,alpha,pvec,d2=-1;
  Point I,G=centre(),OO;
  Vecteur O,u,v;
  bool sign,signG[3];
  int i,j,k,inc,ii=0,idx[2];
  //char verbeux=0;

  if(verbeux<0)  time(&tpsD);
  O.formation_vecteur(OO,C);
  t = ( cst_equ_plan - normale.prod_scalaire(O));
  if(verbeux>1) printf(" cst_eq_plan = %g -- t = %g\n",cst_equ_plan,t);
  I=C+ normale*t ;
  if(verbeux>1) cout<<"Polygone[distance_point] inter I : "<<I[0]<<" "<<I[1]<<" "<<I[2]<<endl;
  //Position du point I / au polygone ? 
  i=0;
  if (fabs(normale[1]) > fabs(normale[0])) i=1;
  if (fabs(normale[2]) > fabs(normale[1])) i=2;   
  j=(i+1)%3;
  k=(j+1)%3;
  for(inc=0;inc<nb_sommets;inc++){
    u.formation_vecteur(G,sommet[inc]);
    v.formation_vecteur(G,sommet[(inc+1)%nb_sommets]);
    signG[inc] =signe(u[j]*v[k]-(u[k]*v[j]));
    u.formation_vecteur(I,sommet[inc]);
    v.formation_vecteur(I,sommet[(inc+1)%nb_sommets]);
    pvec=u[j]*v[k]-(u[k]*v[j]);
    signG[inc]^=signe(pvec);
    if(signG[inc])
      idx[ii++]=inc;
    if(verbeux) printf("sG(%d)=%d - ",inc,(int)signG[inc]);
  }
  if(verbeux<0)  {time(&tpsF); tps1+=tpsF-tpsD;}
  if(verbeux) printf("\n");
  switch(ii){
  case 0:
    d2 = C.dist2(I); break;
  case 1: // printf("ii=1\n");
    u.formation_vecteur(sommet[idx[0]],I);
    v.formation_vecteur(sommet[idx[0]],sommet[(idx[0]+1)%nb_sommets]);
    scal=u.prod_scalaire(v);
    if(scal<=0){
      if(verbeux) printf("1.1 dist(C,S(%d))\n",idx[0]);
      d2 = C.dist2(sommet[idx[0]]);
    }
    else{
      u.formation_vecteur(sommet[(idx[0]+1)%nb_sommets],I);
      v=-v;
      scal=u.prod_scalaire(v);
      if(scal<=0){
	if(verbeux) printf("1.2 dist(C,S(%d))\n",(idx[0]+1)%nb_sommets);
	d2 =  C.dist2(sommet[(idx[0]+1)%nb_sommets]);
      }
      else{
	alpha=scal/v.prod_scalaire(v);
	// Cast oblige le 310798 avec egcs 2.91.52....
	I= Point(sommet[idx[0]]*alpha) + sommet[(idx[0]+1)%nb_sommets]*(1-alpha);
	if(verbeux) {printf("1.3 dist(C,S%dS%d))\n",idx[0],(idx[0]+1)%nb_sommets);I.show();}
	d2 = C.dist2(I);
      }
    }
    break;
  case 2: if(verbeux) printf("ii=2\n");
    u.formation_vecteur(sommet[idx[1]],I);
    v.formation_vecteur(sommet[idx[1]],sommet[idx[0]]);
    scal=u.prod_scalaire(v);
    if(scal>=0){
      alpha=scal/v.prod_scalaire(v);
      I=Point(sommet[idx[1]]*alpha) + sommet[idx[0]]*(1-alpha);
      if(verbeux){ printf("2.1 dist(C,S%dS%d))\n",idx[1],idx[0]);I.show();}
      d2 = C.dist2(I); 
    }
    else{
      v.formation_vecteur(sommet[idx[1]],sommet[(idx[1]+1)%nb_sommets]);
      scal=u.prod_scalaire(v);
      if(scal>=0){
	alpha=scal/v.prod_scalaire(v);
	I=Point(sommet[idx[1]]*alpha) + sommet[(idx[1]+1)%nb_sommets]*(1-alpha);
	if(verbeux){ printf("2.2 dist(C,S%dS%d))\n",idx[1],(idx[1]+1)%nb_sommets);I.show();}
	d2 = C.dist2(I); 
      }
      else{// I=sommet
	if(verbeux) printf("2.3 dist(C,S(%d))\n",idx[1]);
	d2 =  C.dist2(sommet[idx[1]]);
      }
    }
    break;
  }
  if(verbeux<0)  {
    time(&tpsD);
    tps2+=tpsD-tpsF;
    printf("=> tps1=%ld, tps2=%ld\n",tps1,tps2);
  }
  return d2;
}//distance_point()

bool Polygone::appart_sphere(Point & O, double R){
  // sous-estimateur de la distance de O au poly par efficacite
  double d;
  
  d = isob.dist(O) - dGS;
  // printf("appart_sphere()  d=%g (dOG=%g, dGS=%g)\n",d,isob.dist(O),dGS);//I.show();
 return (d<R)? true: false;
}//Polygone::appart_sphere()




#include "primitive.h"
#include "Mmath.h"

#include <ctype.h>




// PRIMITIVE

Primitive :: Primitive(Point p){
	pt_intersection=p;
}
Primitive :: Primitive(Primitive &frere){
  // cout<<"*** Primitive(Primitive &frere)"<<endl;
  pt_intersection=frere.pt_intersection;
  nom=frere.name();
  psurf=frere.surface();
}

// POLYGONE

void Polygone::init(Liste<Point>& liste_sommet,double name=0){
   int i=1;
  
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
  calc_surface();
}//init


Polygone::Polygone(Polygone& frere):Primitive(frere){
  int i;
  // cout<<"*** Polygone(Polygone& frere)"<<endl;
  nb_sommets=i=frere.nb_sommet();
  sommet=new Point[i];
  assert(sommet);
  for (i=0; i<nb_sommets; i++)
    sommet[i]=frere.sommets(i);
  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
}//Cstructeur par copie

char * next_number(char *ch){
  char c;
  int i; 
  
  for (i=0;  (c=ch[i])>0  && isspace(c)!=0; i++) 
    ; //printf(">> ch[%d]=%c\n",i,c); 
  return &(ch[i]); 

}
char * next_nonumber(char *ch){
  char c;
  int i; 
  
  for (i=0; (c=ch[i])>0 && isspace(c)==0; i++)
    ; // printf(">> ch[%d]=%c\n",i,c); 
  return &(ch[i]); 

}


Polygone::Polygone(char* line,double name,reel*mini,reel*maxi){
  //  istrstream Lparag(line);
  register int i,j;
  int ii;
  Point P;
  char *ch;
  double x;
  float Pi;
  
  //  printf("==> line=%s\n", line);
  ch=line;
  for(j=0;j<3;j++){
    mini[j]=99999999.0;
    maxi[j]=-99999999.0;
  }    
  nom=name;  
  // Lparag>>ii;
  ch=next_number(ch);
  sscanf(ch,"%d",&ii);
  ch=next_nonumber(ch);
  nb_sommets=i=ii;
  sommet=new Point[i];
  assert(sommet);
  // qui();
  for (i=0; i<nb_sommets; i++){
    for(j=0;j<3;j++){
      // Lparag>>P[j];
      ch=next_number(ch);
      sscanf(ch,"%f",&Pi);
      ch=next_nonumber(ch);
      P[j]=Pi;
      mini[j]=T_min(P[j],mini[j]);
      maxi[j]=T_max(P[j],maxi[j]);
    }    
    sommet[i]=P;
    // cerr<<" sommet " <<i<<" "; sommet[i].show();
  }
  // rejet test des 3 1st point
  if(nb_sommets>=3) {
    Vecteur u(sommet[0],sommet[1]),v(sommet[0],sommet[2]);
    x=u.norme();
    if(x==0.) {
      fprintf(stderr,"\t=> Comme les sommets 0 et 1 sont egaux, \n");
      mini[0] =1;
      maxi[0]=0;
      correct=false;
      return;
    }
    u/=x;
    x=v.norme();
    if(x==0.) {
      fprintf(stderr,"\t=> Comme les sommets 0 et 2 sont egaux, \n");
      mini[0] =1;
      maxi[0]=0;
      correct=false;
      return;
    }
    v/=x;
    if(sommet[1]==sommet[2]){
      fprintf(stderr,"\t=> Comme les sommets 1 et 2 sont egaux, \n");
      mini[0] =1;
      maxi[0]=0;
      correct=false;
      return;
    }
    if(fabs(u.prod_scalaire(v))>0.9999999999 ){
      fprintf(stderr,"\t=> Comme les sommets 0, 1 et  2  sont aligne's, \n");
      mini[0] =1;
      maxi[0]=0;
      correct=false;
      return; 
    }
    v=u.prod_vectoriel(v);
    x=v.norme();
    if(x<=0){
      fprintf(stderr,"\t=> Norme de (POP1)^(POP2 Nulle \n");
      mini[0] =1;
      maxi[0]=0;
      correct=false;
      return;
    }
  }
  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
  calc_surface();
  correct=true;
}//Polygone

Polygone::Polygone(float(*T)[3],double name,reel*mini,reel*maxi){
  register int i,j;
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
    //cerr<<" sommet " <<i<<" "; sommet[i].show();
  }
  // rejet test des 3 1st point
  if(nb_sommets>=3) {
    Vecteur u(sommet[0],sommet[1]),v(sommet[0],sommet[2]);
    x=u.norme();
    if(x==0.) {
      fprintf(stderr,"\t=> Comme les sommets 0 et 1 sont egaux, \n");
      mini[0] =1;
      maxi[0]=0;
      correct=false;
      return;
    }
    u/=x;
    x=v.norme();
    if(x==0.) {
      fprintf(stderr,"\t=> Comme les sommets 0 et 2 sont egaux, \n");
      mini[0] =1;
      maxi[0]=0;
      correct=false;
      return;
    }
    v/=x;
    if(sommet[1]==sommet[2]){
      fprintf(stderr,"\t=> Comme les sommets 1 et 2 sont egaux, \n");
      mini[0] =1;
      maxi[0]=0;
      correct=false;
      return;
    }
    if(fabs(u.prod_scalaire(v))>0.9999999999 ){
      fprintf(stderr,"\t=> Comme les sommets 0, 1 et  2  sont aligne's, \n");
      mini[0] =1;
      maxi[0]=0;
      correct=false;
      return; 
    }
    v=u.prod_vectoriel(v);
    x=v.norme();
    if(x<=0){
      fprintf(stderr,"\t=> Norme de (POP1)^(POP2 Nulle \n");
      mini[0] =1;
      maxi[0]=0;
      correct=false;
      return;
    }
  }
  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
  calc_surface();
  correct=true;
}//Polygone(reel **)

void Polygone::show(const char* msg,ostream &out)
{   out << msg<<"-Polygone :"; qui();
for (register int i=0; i<nb_sommets; i++)
	{ out << "\t" <<sommet[i][0]<<" "<<sommet[i][1]<<" "<<sommet[i][2]<<endl;
		
	}
}//show
Polygone::~Polygone()
{
	delete sommet;
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

int Polygone::tout_point_inf(reel& position, char& axe){
  short i, compris=1;

  //	cout << "position=" << position << endl;
  //	cout << "axe=" << axe << endl;
  //	cout << "nb_sommets=" << nb_sommets << endl;
  for (i=0; i<nb_sommets; i++) {
    
    //		cout << "sommet[i][axe]=" << sommet[i][axe] << endl;
    if (sommet[i][axe]>position)
      compris=0;
  }
  //	cout << "compris_inf=" << compris << endl;
  return compris;
}

int Polygone::tout_point_sup(reel& position, char& axe){
  short i, compris=1;
 
  //	cout << "position=" << position << endl;
  //	cout << "axe=" << axe << endl;
  for ( i=0; i<nb_sommets; i++){
    //		cout << "sommet[i][axe]=" << sommet[i][axe] << endl;
    if (sommet[i][axe]<position)
      compris=0;
  }
  //	cout << "compris_sup=" << compris << endl;
  return compris;
}

int Polygone::nb_in(reel& amin,reel& amax, char& axe) {
  register int nbp=0,i;
  bool info=false;
  for (i=0; i<nb_sommets; i++) {
    if(sommet[i][axe]>amax) {
      if(nbp<0) {
	cerr<<" (!) Pb Polygone::nb_in : primitive no. "<<nom<<" plus grande que les bornes \n";
	info=true;
      }	  
      nbp++;
    }
    else
      if(sommet[i][axe]<amin) {
	if(nbp>0) {
	  cerr<<" (!) Pb Polygone::nb_in : primitive No. "<<nom<<" plus grande que les bornes\n!";
	info=true;
	}	  
	nbp--;
      }
  }//for sommets
  if(info) {
      cerr<<"--> Axe "<<axe<<" - bornes = ("<<amin<<", "<<amax<<")\n";
      for (i=0; i<nb_sommets; i++) {
	cerr<<"\tP["<<i<<"] = "<<sommet[i][axe]<<endl;
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
    cerr<<"Polygone::calcul_normale_cst_equ : normale nulle ->";
    qui();
    exit(-1);
  }
  normale.normalise();
  //        normale.show();
  /* calcul de la constant de l'equation du plan P*N-D=0 */
  Vecteur OA;
  Point O(0.0, 0.0, 0.0);
  OA.formation_vecteur(O,p1);
  cst_equ_plan=OA.prod_scalaire(normale);
}

Point  Polygone::centre(){
  register int i,j;
  Point G(0,0,0);

  for(i=0;i<nb_sommets;i++)
    for(j=0;j<3;j++){
      //cout<<"Polygone::centre() : "<<(sommet[i])[j]<<endl;
      G[j]+=(sommet[i])[j];
    }
  G=G/(double) nb_sommets;
return G;
}//Polygone::centre()

Vecteur Polygone::azi(){
// renvoie l'azimuth zero du polygone , qui correspond au 1er segment
  Vecteur u(sommet[0],sommet[1]);
  
  return u;
}//Polygone::azi()

// TRIANGLE
Triangle::Triangle(Liste<Point>& liste_sommet,double name=0) : Polygone(liste_sommet,name)
{
  if (nb_sommets != 3)
    { cout << "ERREUR - nombre de sommets incoherent pour un triangle\n"; 
      cout << "nombre de sommets=" << nb_sommets << "\n"; cout.flush();
      exit(1);
    }
//  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
}
Triangle::Triangle(Point &A,Point &B,Point &C,double name=0){
  sommet = new Point[3];
  sommet[0]=A;
  sommet[1]=B;
  sommet[2]=C;
  nom=name;
  nb_sommets=3;
  calcul_normale_cst_equ(sommet[0],sommet [1],sommet[2]);
  calc_surface();
}


inline bool signe(const double& x) 
 { return(x>=0.0);
 } 
bool Polygone::dedans(Point I)
// teste si le point I est "dans" le polugone poly
 { 

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
   for(inc=2;inc<=nb_sommets;inc++)
    { u=v;
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
  Point OO,X(parag.origin());
  O.formation_vecteur(OO,X);
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

//-***************************   Polygone::calc_surface()   *************************
void Polygone::calc_surface(){
  Vecteur u,v, CB,CD;
  //qui();
  switch(nb_sommets){
  case 1: //polygone de faussaire
  case 2:
    cerr<<"Polygone[calc_surface] Irrtum nb_sommet = "<<nb_sommets<<endl;
    exit(-1);
    break;
  case 3: //triangle
    u.formation_vecteur(sommet[0],sommet[1]);
    v.formation_vecteur(sommet[0],sommet[2]);
    u.normalise();
    v.normalise();
    psurf=Macos(u.prod_scalaire(v));
    psurf=sommet[0].dist(sommet[1])*sommet[0].dist(sommet[2])*sin(psurf)/2.0;
    //  raus(psurf<0.0,"Polygone[calc_surface] surface < 0.0 !!??!!?!!!");
    break;
  case 4:
    //rectangle
    u.formation_vecteur(sommet[0],sommet[1]);
    v.formation_vecteur(sommet[0],sommet[3]);
    CB.formation_vecteur(sommet[2],sommet[1]);
    CD.formation_vecteur(sommet[2],sommet[3]);
    if(u.prod_scalaire(v)==CB.prod_scalaire(CD)){
      psurf=(u.prod_vectoriel(v)).norme();
      //cerr<<"Polygone[calc_surface]rectangle\n";
      break;
    }
  default:
    Point G=centre();
    psurf=0.0;
    for(register int i=0;i<nb_sommets;){
      u.formation_vecteur(G,sommet[i++]);
      v.formation_vecteur(G,sommet[(i<nb_sommets)?i:0]);
      psurf+=(u.prod_vectoriel(v)).norme();
      //cerr<<"Polygone[calc_surface]poly qcque\n";
    }
    break;
  }//switch   
 }// Polygone::calc_surface() 


//***************************   Triangle::calc_surface()   *************************
void Triangle::calc_surface(){
  Vecteur AB(sommet[0],sommet[1]), AC(sommet[0],sommet[2]);   

  AB.normalise();
  AC.normalise();
  psurf=Macos(AB.prod_scalaire(AC));
  //   cout <<"@@@@@@@@@@ M_PI_2 = "<< M_PI_2<<" - alpha = "<<psurf<<endl;
  psurf=sommet[0].dist(sommet[1])*sommet[0].dist(sommet[2])*sin(psurf)/2.0;
  //raus(psurf<0.0,"Triangle[calc_surface] surface < 0.0 !!??!!?!!!");
}// Triangle::calc_surface() 


#ifdef _SEGURA
#define BIG -999999
  double L8=BIG;
  Point O,Q1;


double signo3D ( Point p1 , Point p2 , Point p3 , Point p4 ){
  double suma;
  
  suma = p1[0]*( (p2[1]*p3[2]) + (p2[2]*p4[1]) + (p3[1]*p4[2])
                -(p3[2]*p4[1]) - (p4[2]*p2[1]) - (p2[2]*p3[1]) 
	       ) 
	-p2[0]*(  (p1[1]*p3[2]) + (p1[2]*p4[1]) + (p3[1]*p4[2])
                 -(p3[2]*p4[1]) - (p1[2]*p3[1]) - (p4[2]*p1[1])
	       ) 
        +p3[0]*(  (p1[1]*p2[2]) + (p2[1]*p4[2]) + (p1[2]*p4[1])
                 -(p2[2]*p4[1]) - (p4[2]*p1[1]) - (p1[2]*p2[1])
	       ) 
       -p4[0]*(  (p1[1]*p2[2]) + (p2[1]*p3[2]) + (p1[2]*p3[1])
                -(p2[2]*p3[1]) - (p3[2]*p1[1]) - (p1[2]*p2[1])
	      ) ;

        return suma;
}
void signo3D2 ( Point p1a, Point p1b , Point p2 , Point p3 , Point p4, double&s1a,  double&s1b){
  double tmp[4];

   tmp[0] =  (p2[1]*p3[2]) + (p2[2]*p4[1]) + (p3[1]*p4[2])
            -(p3[2]*p4[1]) - (p4[2]*p2[1]) - (p2[2]*p3[1]); 
   tmp[1] =  (p2[0]*p3[2]) + (p2[2]*p4[0]) + (p3[0]*p4[2])
            -(p3[2]*p4[0]) - (p4[2]*p2[0]) - (p2[2]*p3[0]); 
   tmp[2] =  (p2[0]*p3[1]) + (p2[1]*p4[0]) + (p3[0]*p4[1])
            -(p2[1]*p3[0]) - (p4[1]*p2[0]) - (p3[1]*p4[0]); 
   tmp[3] =  p4[2]*(p2[0]*p3[1]) +  p3[2]*(p2[1]*p4[0]) + p2[2]*(p3[0]*p4[1])
            -p2[2]*(p3[1]*p4[0]) -  p3[2]*(p4[1]*p2[0]) - p4[2]*(p2[1]*p3[0]); 

   // printf("~~> tmp=%g, %g,%g,%g\n",tmp[0],tmp[1],tmp[2],tmp[3] );
  s1a = (p1a[0]*tmp[0] - p1a[1]*tmp[1] + p1a[2]*tmp[2] - tmp[3]);
  s1b = (p1b[0]*tmp[0] - p1b[1]*tmp[1] + p1b[2]*tmp[2] - tmp[3]);
}

double Triangle::intersect(Param_Inter& parag,Point *I) {
  // algo de Segura&Leito, C&G98
  double t,a,b,c,h1,h2;
 
  if(parag.prems){
    O=parag.origin();
    Q1= O + parag.direct()*L8 ;
    parag.prems=false;
  }
  a = signo3D (O ,Q1 , sommet[2] , sommet[0]); 
  b = signo3D (O ,Q1 , sommet[1] , sommet[2]); 
  c = signo3D (O ,Q1 , sommet[0] , sommet[1]); 
 
  if( (a>0 && b>0 && c>0) || (a<0 && b<0 && c<0) ){// Il y a intersection
    // printf("*** h1=%g, h2=%g --> t=%g\n",h1,h2,fabs(h1/(h1-h2))*L8);
    signo3D2(O,Q1, sommet[0] , sommet[1], sommet[2] , h1,h2);
    t=fabs(h1/(h1-h2))*L8;
    *I=parag.origin()+ parag.direct()*t ;
    //printf("*** h1=%g, h2=%g  => t=%g\n",h1,h2,t);
    return t;
  }
  else
    return 4e9;

}//Triangle::intersect()

Triangle::Triangle (char*str,double nome,reel* mini,reel* maxi):Polygone (str,nome,mini,maxi)
       {L8=(L8<(maxi[1]-mini[1]))?maxi[1]-mini[1]:L8;}
#else

Triangle::Triangle(char*str,double nome,reel* mini,reel* maxi):Polygone (str,nome,mini,maxi)
{}
#endif

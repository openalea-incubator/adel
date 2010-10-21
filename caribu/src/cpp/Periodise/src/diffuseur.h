#ifndef DIFFUSEUR
#define DIFFUSEUR

#include "primitive.h"
#include "actop.h"
#define ORDRE_MAX 6
//#include "outils.h"
//modif pour le cas infini : 050496 : 
//   creation d'ne nouvelle hierachie

enum Visiblite {pasvu, ptet, vu};

// DIFFUSEUR: classe mere virtuelle
class Diffuseur{
public:
  // Actop * optback;
  Diffuseur() {}
  //{cerr<<"Diffuseur [constructeur] Warning : virtual class !\n";}
  virtual  Primitive& primi()=0;
  virtual  Param_Inter interact(Param_Inter&,bool,int)=0;
  virtual  double lux(Param_Inter&)=0;
  virtual bool isreal()=0;
  virtual bool actif()=0;
  virtual bool isopaque()=0;
  virtual void maj_stat(int)=0;
  virtual void calc_stat(int,double Stoit)=0;
  virtual void norm_surf(int,double)=0;
  virtual double nrj(int)=0;
  virtual double VEa(int)=0;
  virtual double  rad(char)=0;
  virtual double  ncontact(char)=0;
  virtual double  Vrad(char)=0;
  virtual double rad1(char)=0;
  virtual double inc(char)=0;
  virtual double Vrad1(char)=0;
  virtual Visiblite is_visible()=0;
  virtual Diffuseur&  operator ++(int)=0;
  virtual Diffuseur&  operator --(int)=0;
  virtual double intersect(Param_Inter parag,Point *I)=0;
  virtual void show()=0;
  virtual ~Diffuseur(){}; 
};

//Class DiffReel: mere des diff normaux
class DiffReel:public Diffuseur{
protected:
  double Eabs[2],E1a[2],E2a[2]; //2 valeurs=ordre 0 et ordres>1 de rediffusion
  Primitive *prim;
  Actop *opti;
  long vis_max;
  long vis_cur;
  //Alea rambo;
public:
  DiffReel()
  { cout<<"DiffReel [constructeur] Warning : no parameters!\n";}
  DiffReel (Primitive *pprim,Actop *actop)
// : prim(p),opti(actop)      {vis_max=0;vis_cur=0; }
  { //cout<<"DiffReel[cstructor]debut\n"; pprim->show();
    opti=actop; vis_max=0; vis_cur=0;
    prim=pprim;
    //cout<<"DiffReel[cstructor]fin\n";
  }
  Primitive& primi()
  { return *prim; }
  // int nb_som(){ return prim->nb_sommet(); }
  bool isreal()
  {return true;}
  bool actif()
  {return true;}
  //  virtual  Param_Inter interact(Param_Inter&,bool,int )=0;
  virtual  double lux(Param_Inter&)=0;
  virtual bool isopaque()=0;
  virtual void maj_stat(int)=0;
  virtual void calc_stat(int,double Stoit)=0;
  virtual void norm_surf(int,double)=0;
  virtual double nrj(int)=0;
  virtual double VEa(int)=0;
  //  virtual double  rad(char)=0;
  // virtual double  Vrad(char)=0;
  Visiblite is_visible()
  { if(vis_cur==0) return pasvu;
  return ((vis_cur<vis_max)?ptet:vu);
  } 
  DiffReel& operator ++(int)
  {vis_max++; vis_cur++; return *this;}
  DiffReel&  operator --(int)
  {vis_cur--; return *this;}
  double intersect(Param_Inter parag,Point *I)
  {return (prim->intersect(parag,I));}
  void show() {prim->show(); 
  cout <<"curvis = "<< vis_cur<< " - maxvis = "<< vis_max<<endl;
  }
  //void fini() {opti->fini();}
};

class DiffO :public DiffReel{
private:
  double Einc,E1i,E2i;
  double  radio,radio2,radio1;
  double  B1st,B1st2,B1st1;
  double nk; //nbre de rayon reemis
public:
  ~DiffO() {delete prim;}// On ne delete pas *opti car il est partage avec d'autres diff
  DiffO(Primitive *p,Actop *actop): DiffReel(p,actop)
  { Eabs[0]=Eabs[1]=Einc=0.0;
  E1a[0]=E1i=E2a[0]=E2i=0.0;
  E1a[1]=E2a[1]=0.0;
  radio=radio1=radio2=0.0;
  B1st=B1st1=B1st2=0.0;
  nk=0.0;
  }
  Param_Inter interact(Param_Inter &parag,bool,int);
  double lux(Param_Inter& );
  bool isopaque() {return true;}
  void maj_stat(int);
  void calc_stat(int,double Stoit);
  void norm_surf(int,double);
  double nrj(int ordre) {return Eabs[ordre];}
  double VEa(int ordre) {return E2a[ordre];}
  double rad(char face) {face=0; return radio;}
  double  ncontact(char face){face=0;return nk;}
  double Vrad(char face) {face=0;return radio2;}
  double rad1(char face) {face=0;return B1st;}
  double inc(char face){face=0;return Einc;}
  double Vrad1(char face) {face=0;return B1st2;}
};

class DiffT :public DiffReel{
private:
  Actop * optback;
  double Einc[2],E1i[2],E2i[2]; // 0 : face sup - 1 : face inf 
  double   radio[2],radio1[2],radio2[2];      // 0 : face sup - 1 : face inf
  double   B1st[2],B1st1[2],B1st2[2];      // 0 : face sup - 1 : face inf
  double nk[2]; //nbre de rayon reemis
  Actop *popt(Face side){
    return (side==sup)? opti:optback;
  }
public:
 
  ~DiffT() {delete prim;}
  DiffT(Primitive *p,Actop *actop, Actop * backopt): DiffReel(p,actop){
    Eabs[0]=Eabs[1] =Einc[0]=Einc[1]=0.0;
    E1a[0]=E1a[1]=E2a[0]=E2a[1]=0.;
    E1i[0]=E1i[1]=E2i[0]=E2i[1]=0.0;
    radio[0]=radio1[0]=radio2[0]=0.0;
    radio[1]=radio1[1]=radio2[1]=0.0;
    B1st[0]=B1st1[0]=B1st2[0]=0.0;
    B1st[1]=B1st1[1]=B1st2[1]=0.0;
    nk[0]=nk[1]=0.0;
    // printf("DiffT::cstructor (%.0lf):: avant actop->koi()",p->name());
    //opti->koi();
    optback=backopt;
    //printf("DiffT::cstructor (%.0lf):: avant optback ->koi()",p->name());
    //optback->koi();
  }
  // virtual  Param_Inter interact(Param_Inter&,bool,int )=0;
  Param_Inter interact(Param_Inter&,bool,int);
  double lux(Param_Inter&);   
  bool isopaque() {return false;}
  void maj_stat(int nb_ray);
  void calc_stat(int nb_iter,double Stoit);
  void norm_surf(int nb_ray,double);
  double nrj(int ordre) {return Eabs[ordre];}
  double VEa(int ordre) {return E2a[ordre];}
  double rad(char face) {return radio[face];}
  double inc(char face){return Einc[face];};
  double  ncontact(char face =0){return nk[face];}
  double Vrad(char face) {return radio2[face];}
  double rad1(char face) {return B1st[face];}
  double Vrad1(char face) {return B1st2[face];}
};


/***********************************************************
 *   class Diff8 : pour traiter l'infini sur les bords!    *
 *                                                         *
 *   la prim n'intervient que pour calc le pt d'intersect  *
 *   pour interact. lux-mat, on n'utilise que la notion de *
 *   de direction, inchange par translation, donc on peut  *
 *   utiliser le repere du diffuseur originel...           *
 ***********************************************************/

class Diff8 :public Diffuseur{
protected:
  Primitive *prim;
  Diffuseur * diff;
  long vis_max;
  long vis_cur;
  //Alea rambo;
public:
  ~Diff8(){delete prim;}
  Diff8()
  { cerr<<"Diff8 [constructeur] Warning : no parameters!\n";}
  Diff8 (Diffuseur * pdif,Vecteur & delta);
  Primitive& primi(){
    return *prim; }
  // int nb_som(){ return prim->nb_sommet(); }
  bool isreal(){
    return false;}
  bool actif()
  {return true;}
  Param_Inter interact(Param_Inter &parag,bool inside,int order){
    return diff->interact(parag,inside,order);}
  double lux(Param_Inter& parag) {
    return diff->lux(parag);}
  bool isopaque() {
    return diff->isopaque();}
  void maj_stat(int nb_ray) {nb_ray=0;}
  void calc_stat(int nb_iter,double Stoit) {nb_iter=0; Stoit=0;}
  void norm_surf(int nb_ray,double Stoit) {nb_ray=0; Stoit=0;}
  
  double nrj(int ordre){
    return diff->nrj(ordre);}
  double VEa(int ordre) {
    return diff->VEa(ordre);}
  double  rad(char face){
    return diff->rad(face);}
  double  inc(char face){
    return diff->inc(face);}
 double  ncontact(char face =0){
   return diff->ncontact(face);}
  double  Vrad(char face){
    return diff->Vrad(face);}
  double  rad1(char face){
    return diff->rad1(face);}
  double  Vrad1(char face){
    return diff->Vrad1(face);}
  Visiblite is_visible(){
    if(vis_cur==0) return pasvu;
    return ((vis_cur<vis_max)?ptet:vu); } 
  Diff8&  operator ++(int){
    vis_max++; vis_cur++; return *this;}
   Diff8&  operator --(int){
    vis_cur--; return *this;}
  //geometrie ->new prim
  double intersect(Param_Inter parag,Point *I){
    return (prim->intersect(parag,I));}
  void show() {
    prim->show();
    cout <<"curvis = "<< vis_cur<< " - maxvis = "<< vis_max<<endl;
  }
  //void fini() {opti->fini();}
};

/***********************************************************
 *   class DiffP(Licor) ,DiffP2(Solem) : pour simuler un capteur passif          *
 *                                                         *
 *   la prim n'intervient que pour stocker l'energie       *
 *   qui la traverse ; ele ne modifie ni le poids,         *
 *   ni la direction du rayon. Elle stocke l'energie       *
 *   comme AnaFlux et l'ecrit en append dans des fichiers  *
 *   cumul : pvf.cum et profpen.cum                        *
 ***********************************************************/

class DiffP:public Diffuseur{
protected:
  Primitive *prim;
  float Z;
public:
  //variables publiques
  double phisun;//pour bien calculer le phi des flux dir.
  Tabdyn<double,4> fluxm;//[2],inf/sup, [6] ordres, [4] 1/2plans,[18] teta;
  Tabdyn<double,2> phim;//[2],inf/sup, [6] ordres : valeur du flux hemispherique
  double pene;//[nb] couches : penetration du direct
  //fonctions-membres
  DiffP()
  { cout<<"DiffP [constructeur] Warning : no parameters!\n";}
  DiffP (Primitive *pprim);
  ~DiffP();
  Primitive& primi()
  { return *prim; }
  // int nb_som(){ return prim->nb_sommet(); }
  bool isreal()
  {return true;}
   bool actif()
  {return false;}
  Param_Inter interact(Param_Inter&,bool,int );
  double lux(Param_Inter&)
  {return 0.;}
  bool isopaque()
  {return false;} 
  void maj_stat(int){}
  void calc_stat(int i,double Stoit) {i=0; Stoit=0;}
  virtual void norm_surf(int i ,double x){i=0;x=0;}
  double nrj(int ordre){ordre=0; return -1;}
  double VEa(int ordre){ordre=0; return -1;}
  double  rad(char face){face=0;return 0;}
  double  ncontact(char face){face=0;return 0;}
  double  Vrad(char face){face=0;return 0;}
  double inc(char face){face=0;return 0;}
  double rad1(char face){face=0;return 0;}
  double Vrad1(char face){face=0;return 0;} 
  Visiblite is_visible()
  { return pasvu;
  } 
  DiffP& operator ++(int)
  {return *this;}
  DiffP&  operator --(int)
  {return *this;}
  double intersect(Param_Inter parag,Point *I)
  {return (prim->intersect(parag,I));}
  void show() {prim->show(); 
  }
  //void fini() {opti->fini();}
};

class DiffP2:public Diffuseur{
protected:
  Primitive *prim;
public:
  //variables publiques
  double *E;//Direct, Diffus
  int No;
  //fonctions-membres
  DiffP2()
  { cout<<"DiffP2 [constructeur] Warning : no parameters!\n";}
  DiffP2 (Primitive *pprim);
  DiffP2 (Primitive *pprim,int nbo);
  ~DiffP2();
  Primitive& primi()
  { return *prim; }
  // int nb_som(){ return prim->nb_sommet(); }
  bool isreal()
  {return true;}
   bool actif()
  {return false;}
  Param_Inter interact(Param_Inter&,bool,int );
  double lux(Param_Inter&)
  {return 0.;}
  bool isopaque()
  {return false;} 
  void maj_stat(int){}
  void calc_stat(int i,double Stoit) {i=0;Stoit=0;}
  virtual void norm_surf(int,double){}
  double nrj(int ordre){ordre=0; return 0;}
  double VEa(int ordre){ordre=0; return 0;}
  double  rad(char face){face=0;return 0;}
  double  ncontact(char face){face=0;return 0;}
  double  Vrad(char face){face=0;return 0;}
  double inc(char face){face=0;return 0;}
  double rad1(char face){face=0;return 0;}
  double Vrad1(char face){face=0;return 0;} 
  Visiblite is_visible()
  { return pasvu;
  } 
  DiffP2&  operator ++(int)
  {return *this;}
   DiffP2&  operator --(int)
  {return *this;}
  double intersect(Param_Inter parag,Point *I)
  {return (prim->intersect(parag,I));}
  void show() {prim->show(); 
  }
  //void fini() {opti->fini();}
};


#endif





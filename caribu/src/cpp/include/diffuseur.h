#ifndef R_DIFFUSEUR
#define R_DIFFUSEUR

#include "primitive.h"
#include "verbose.h"
#include "actop.h"
//#include "outils.h"
// DIFFUSEUR

class Diffuseur {
protected:
  Primitive *prim;
  Actop *opti;
public:
  static unsigned int idx;
  signed char acv;// A Cheval
  Diffuseur(){ //cout<<"Diffuseur [constructeur] Warning : no parameters!\n";
  }
  Diffuseur (Primitive *pprim,Actop *actop) : prim(pprim),opti(actop) {}
  Primitive& primi() { return *prim; }
  virtual  Vecteur normal() {return prim->normal();}
  Point centre() {return prim->centre();}
  double surface() {return prim->surface();}
  double name() {return prim->name();}
  Vecteur azi0(){return prim->azi();}//azimuth zero
  virtual bool isopaque()=0;
  virtual bool isreal()=0;
  double intersect(Param_Inter parag,Point *I)
  {return (prim->intersect(parag,I));}
  void  show(char *texte="",ostream& out=cout) // montre!
  { prim->show(texte,out);  }
  virtual  unsigned int num()=0;
  virtual  void togle_face()=0;
  virtual void active(Vecteur&)=0;
  virtual void active(unsigned char)=0;
  virtual void activ_num(unsigned int)=0;
  virtual unsigned char face()=0;
  virtual double rho()=0;
  virtual double tau()=0;
  // amie
  friend int maxE(Diffuseur*,Diffuseur*); // utilise par QuickSort (TabDyn, ListeD)  
  // renvoie -1 si E1>E2, 0 si E1=E2, 1 si E1<E2 (Ei delta energie du difuseur i
  friend ostream& operator << (ostream& out,Diffuseur & diff);
  // cas du patch voir si on fait une hierachie double (class Patch)
  void select_patch(int p) {};
  int nb_patch() {return 1;}
};//Diffuseur


class DiffO :public Diffuseur{
protected:
  unsigned int no;
public:
  DiffO() {}
  DiffO(Primitive *p,Actop *actop): Diffuseur(p,actop) { no=idx++; }
  bool isopaque() {return true;}
  bool isreal()   {return true;}
  double rho() {return  opti->rho();}
  double tau() {return opti->tau();} 
  unsigned int num() {return no;}
  void togle_face() {}
  void active(Vecteur &dir) {}
  void active(unsigned char cefa) {}
  void activ_num(unsigned int nummer) {}
  unsigned char face() { return 0;}
};//DiffO

class DiffP :public DiffO{ //Capteur virtuel

public:
  DiffP(Primitive *p,Actop *actop): DiffO(p,actop){ }
  DiffP(Primitive *p) {Actop *actop; actop=new Lambert(); prim=p; opti=actop;no=idx++;} 
  bool isreal(){return false;}
};//DiffP


class DiffT :public Diffuseur{
private:
  Face actif; // 0 : face sup - 1 : face inf
  unsigned int no[2];
  Actop * optback;
  Actop *popt(Face side)
  {return (side==sup)? opti:optback;}
public:
  DiffT(Primitive *p,Actop *actop, Actop * backopt): Diffuseur(p,actop){
    optback=backopt;
    no[0]=idx++;
    no[1]=idx++;
    actif=sup;//default comme ca!
  }
  bool isopaque() {return false;}
  bool isreal()   {return true;}
  Vecteur normal(){
    if(actif==sup) return prim->normal();
    else           return prim->normal()*-1.0;
  }
 
  double rho() {return  popt(actif)->rho();}
  double tau() {return popt(actif)->tau();} 
 
  unsigned int num() {return no[actif];}
  void togle_face() {actif =1-actif;}
  void active(Vecteur &dir) {
    if(  dir.prod_scalaire(prim->normal())<0)
      actif=sup;
    else actif = inf;
  }
  void active(unsigned char cefa) {
    actif=(cefa==1)?inf : sup;
  }
  void activ_num(unsigned int nummer) {
    if(nummer!=num()){
      togle_face();
      if(nummer!=num()) {
	Ferr<<" Stronzo faux nummer dans DiffT::activ_num"<<'\n';
	exit(30);}
    }
  }//activ_num()
  unsigned char face() {return actif;}
};//DiffT
#endif



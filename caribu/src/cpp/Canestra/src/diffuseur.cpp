#include "diffuseur.h"

//init du compteur de faces visibles (static varaible)
unsigned int Diffuseur::idx=0;


//-*************** face() ***********************************
inline Face face(double scal) {
  return ((scal>0.0)?inf:sup);
}

//-*************** operator << ( Diffuseur) *****************
ostream& operator << (ostream& out,Diffuseur & diff){
  /*  if(&diff)
        diff.show(out); 
  else  out << 0;*/
   return out; 
}//operator << ( Diffuseur)

#ifdef _NRJ
//-*************** friend Emax() ****************************
int maxE(Diffuseur *g,Diffuseur *d){
  register double grad,drad;

  /*  grad=g->dRad();
  drad=d->dRad();
  if  (grad>drad)  return -1;
  else
    if(grad==drad) return  0;
    else         */  return  1;
}//maxE()
  
//-*************** DiffO::Radmax() *****************************
double DiffO::Radmax(){
  return rad;
}//DiffO:Radmax()

//-*************** DiffT::Radmax() *****************************
double DiffT::Radmax(){
  if(rad[inf]>rad[sup])
    actif = inf;
  else
    actif = sup;
  return rad[actif];
}//DiffT:Radmax() 

//-*************** DiffT::Ecum() ************************

double DiffT::Ecum(){
  return (E[0]+E[1])*surface();
}// DiffT::Ecum() 

//-*************** DiffO::Ecum() ************************

double DiffO::Ecum(){
  return E*surface();
}// DiffO::Ecum() 

//-************** DiffO::calc_stat() ************************

void DiffO::norm_surf()
 { 
 // calcul de la moyenne et de la variance
 // calcul de la moyenne et de la variance
    double surf;
    
    surf=prim->surface();
    Eabs/= surf;
    Einc/= surf;

 }//DiffO::norm_surf()

//-************** DiffT::norm_surf() ************************

void DiffT::norm_surf(){
  // normalise par la surface [W] -> [W/m*m] 
    double surf;
    
    surf=prim->surface();
    Eabs   /= surf;
    Einc[0]/= surf;
    Einc[1]/= surf;
 }//DiffT::norm_surf()


#endif

















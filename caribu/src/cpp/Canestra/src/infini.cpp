/*                  Infini.C - MC96
   Infinitise une projection parallele par pavage recursif

   MC09: Bug fixed: bias on infinitise() for high zenithal angles and particular azimuthal angles.
 */

#include <iostream> //
using namespace std;

#include <cmath>

#define _INFINI
#include "infini.h"
#include "image.h"
#include "chrono.h"

// global variable
static Tabdyn<void *,2> Zdat0;
static Tabdyn<REELLE,2> Zbuf0;
static void ***Zdat8;
static REELLE **Zbuf8;
static double cdist;
static int **T, Ti,Tj, tr[8];
static bool duplik;
static Chrono chrono;

//proto
void lateral(int x,int y, char idx) ;
void pave(int x,int  y, char idx) ;

//local function
bool out(int i,int j) {
  register int k,l;
  bool dehors[2]={true,true};


  //MC09: zarbi /k : A revoir !!!
  //for(k=0;k<4 && dehors;k++) {
  //Toit dehors
  //Bug au zenith eleve'  - MC09
  // dehors[0]&=( (T[k][0]+i)>=Ti || (T[k][0]+i)<0);
  //dehors[1]&=( (T[k][1]+j)>=Tj || (T[k][1]+j)<0);
  //}//MC09
    
  //toute l'img tranlatee dehors    
  // - MC09
  // dehors[0]&= fabs(i)>Ti;
  // dehors[1]&= fabs(j)>Tj ;
  dehors[0]&= fabs(i)>=Ti;
  dehors[1]&= fabs(j)>=Tj ;
  // BUG BUG BUG BUG  - MC09
  //MC96   return dehors[0]||dehors[1]; // BUG BUG theta fort > 70° !!! - MC09
  // MC09 return (dehors[0] && dehors[1]); // Bug qd phi=0°, boucle infinie...
  return (dehors[0] && dehors[1]) || (i==0 && dehors[1]) || (dehors[0] && j==0) ;

 
}//out()

void zbuf(int i,int j) {
  int k,l,max[2]={Ti,Tj},min[2]={0,0};
  double d;
  
  if(i>=0) min[0]=i; else max[0]+=i; 
  if(j>=0) min[1]=j; else max[1]+=j;
  for(l=min[1];l<max[1];l++)
    for(k=min[0];k<max[0];k++) {
      d=Zbuf0(k-i,l-j)+cdist*j;
      if( d < Zbuf8[k][l] ) {
	Zbuf8[k][l]=(REELLE) d;
	Zdat8[k][l]= duplik?Zdat0(k-i,l-j):NULL;
      }
    }      
}//zbuf()

void lateral(int x,int y, char idx) {
 if(out(x,y)) return;
 if(verbose>1) printf("lateral (%d,%d)-",x,y);
 zbuf(x,y);
 lateral(x+tr[idx],y+tr[idx+1],idx);
}//lateral()


void pave(int x,int  y, char idx ) {
  if(out(x,y)) return;
  if(verbose>1)printf("\n pave(%d,%d) :",x,y);
  if( !(x==0 && y == 0))
    zbuf(x,y);
  lateral(x+tr[2],y+tr[3],2);      // gauche
  lateral(x+tr[4],y+tr[5],4);      // droite
  pave(x+tr[idx],y+tr[idx+1],idx); // tout droit
}//pave()

//-***  Exported Functions : infinitise()
void infinitise(void ***Zprim, REELLE **Zbuf,
		double cste_dist,int** roof,int Tx, int Ty,bool dupli) {
  int i,j;
  if(verbose>1){
    chrono.Start();
    cout<<"* infinitise(): DEBUT";
  }
  //init
  Zdat0.alloue(Tx,Ty);
  Zbuf0.alloue(Tx,Ty);
  // parameters --> global variables
  Zdat8=Zprim;
  Zbuf8=Zbuf;
  T=roof;
  duplik=dupli;
  cdist=cste_dist;
  Ti=Tx;Tj=Ty;
  for(j=0;j<Ty;j++)
    for(i=0;i<Tx;i++) {
      Zdat0(i,j)=Zprim[i][j];
      Zbuf0(i,j)=Zbuf[i][j];

    }
  i=j=0;

  /* for(int ii=0; ii<4; ii++)
     for(int jj=0;jj<2;jj++)
     printf("---->  T(%d,%d) = %d\n",ii,jj,T[ii][jj]);
  */
  /* //MC 1996 : Ok mais pb d'arrondi possible donc...
     tr[0]=T[3][0]-T[0][0];
     tr[1]=T[3][1]-T[0][1];
     tr[2]=T[1][0]-T[0][0];
     tr[3]=T[1][1]-T[0][1];
  */
  //    printf("\nOLD tr(0-3)=(%d, %d, %d, %d)\n", tr[0], tr[1], tr[2], tr[3]);

  //MC09: debug / fort angles zenitaux => ne chga erien mais plus rigoureux
  tr[0]=round( ( (T[3][0]-T[0][0]) +(T[2][0]-T[1][0]) )/2.);
  tr[1]=round( ( (T[3][1]-T[0][1]) +(T[2][1]-T[1][1]) )/2.);
  tr[2]=round( ( (T[1][0]-T[0][0]) +(T[2][0]-T[3][0]) )/2.);
  tr[3]=round( ( (T[1][1]-T[0][1]) +(T[2][1]-T[3][1]) )/2.);
  //fin dbug MC09
 
  tr[4]=-tr[2];
  tr[5]=-tr[3];
  tr[6]=-tr[0];
  tr[7]=-tr[1];

  // printf(" NEW tr(0-3)=(%d, %d, %d, %d)\n", tr[0], tr[1], tr[2], tr[3]);

  // pave tq le toit et ses translates pavent tte l'image 
  pave(i,j,0); //up
  i+=tr[6];
  j+=tr[7];
  pave(i,j,6); //down
  //Mr Propre
  Zdat0.free();
  Zbuf0.free();
  if(verbose>1){
    chrono.Stop();
    cout<<"\n::::>  Infinitisation en "<<chrono<<endl;
    fflush(stdout);
  }
}//infinitise()

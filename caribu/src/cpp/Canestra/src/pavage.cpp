
#include <iostream>
using namespace std;

#include <cmath>

#include "image.h"
void lateral(Image *pict,double x,double y, double dx, double dy, int col) ;
void pave(Image *pict, double x,double y,  cord, int col) ;

bool out(int i,int j) {
  if(i>-1 && i<10 && j>-1 && j<10)
    return false;
  else return true;
}//in()

void lateral(Image &pict,double x,double y, double dx, double dy, int col) {
  int i,j;
  i=(int) x;
  j=(int) y;
 if(out(i,j)) return;
 printf("lateral%i (%d,%d)-",(int)(dir[0]*100),i,j);
 pict.maj(i,j,col);
 // calcul de la direction
 while (((int) x==i) &&((int) y==j)) {
  x+=dx[0];
  y+=dy;
 }
 lateral(pict,x,y,dx,dy,,col+10);
}//lateral()


void pave(Image &pict,double x,double y, double dx, double dy, int col) {
  int i,j;
  double dirold [2],ori[2];
  
  ori [0] = x;
  ori [1]=  y;
  dirold [0] = dx;
  dirold [1]=  dy;  
  i=(int) x;
  j=(int) y;
 if(out(i,j)) return;
 printf("\n pave(%d,%d) :",i,j);
 pict.maj(i,j,col);
 col+=100;
 // calcul de la direction lat 1
 dir[0]=-dirold[1];
 dir[1]=dirold[0];
 while (((int) cord[0]==i) &&((int) cord[1]==j)) {
  cord[0]+=dir[0];
  cord[1]+=dir[1];
 }
 lateral(pict,dir,cord,col);
 // calcul de la direction lat 1
 cord[0]=ori[0];
 cord[1]=ori[1];
 dir[0]=dirold[1];
 dir[1]=-dirold[0];
 while (((int) cord[0]==i) &&((int) cord[1]==j)) {
  cord[0]+=dir[0];
  cord[1]+=dir[1];
 }
 lateral(pict,dir,cord,col);
 col-=90;
 // calcul de la direction ttdroit
 cord[0]=ori[0];
 cord[1]=ori[1];
 dir[0]=dirold[0];
 dir[1]=dirold[1];
 while (((int) cord[0]==i) &&((int) cord[1]==j)) {
  cord[0]+=dir[0];
  cord[1]+=dir[1];
 }
 pave(pict,dir,cord,col);
}//pave()
main() {
  double dir[2], coord[2],dir0[2], ori[2],norm;
  Image pict(10,10,"pave.ppm");
  pict.raz(0);

  cout<<" coord du point de depart x,y :\n ";
  cin>>ori[0]>>ori[1];
  cout<<"\n coord de la direction  : \n";
  cin>>dir0[0]>>dir0[1];
  norm=sqrt(dir0[0]*dir0[0]+dir0[1]*dir0[1]);
  dir0[0]/=norm;
  dir0[1]/=norm;
  dir[0]=dir0[0];
  dir[1]=dir0[1];
  coord[0]=ori[0];
  coord[1]=ori[1];
  pave(pict,dir,coord,30);
  dir0[0]*=-1;
  dir0[1]*=-1;
   dir[0]=dir0[0];
  dir[1]=dir0[1];
  coord[0]=ori[0];
  coord[1]=ori[1]; 
  pave(pict,dir,coord,60);
  pict.sauve();
}//main()
/*
to compile
 g++ -o pave -I ../bibliotek -I /home/chelle/rezo/dev/mathlib/include/ pavage.C  $EXE/outils.o $EXE/T_geometrie.o -lm

 */

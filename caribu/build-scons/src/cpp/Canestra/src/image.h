//                      image.H

#ifndef _Image
#define _Image
 
#include <math.h>
#include <stdlib.h> // pour exit
#include <iostream> //
#include <fstream>  //
#include <iomanip>  //

#include "T_utilitaires.h"
#include "outils.h"

// gestion d'image au format PPM visible par XV

class Image {
 private:
   Tabdyn<double,2> btm;
   double bitmax;
   char filename[80];
 public:
    Image(){}
    Image(int x,int y,char * nom="out.ppm")
      {init(x, y,nom);}

    // typage AH 02 2001
    void init(int x,int y,char * nom="out.ppm")
      { btm.alloue(x,y);strcpy(filename,nom);bitmax=-9.9e10;}

    Image & operator =(Image &img)
      { cout<<"Image[=]  DEBUT \n";cout.flush();
        Image  *pim= new Image(img.btm.maxi()[0],img.btm.maxi()[1],img.filename);
/*        strcpy(filename,img.filename); cout<<"Image[=] :"<<filename<<endl;
        btm.alloue(img.btm.maxi()[0],img.btm.maxi()[1]);
*/
    cout<<"Image[=]  avant copie des variables de btm \n";cout.flush(); 
       for(register int i=0;i<pim->btm.maxi()[0];i++)
         for(register int j=0;j<pim->btm.maxi()[1];j++)
          pim-> btm(i,j)=   img.btm(i,j);
        return *pim;
      } 
   double maxval(){
     double val;
     register int i,j;
     
     bitmax=-9.9e20;
     for(j=0;j<btm.maxi()[1];j++)
       for(i=0;i<btm.maxi()[0];i++)
       { val= btm(i,j);
       bitmax=(bitmax<val)? val : bitmax;
       }
     return bitmax;  
   }      
  void fixmax(double maxxi) {bitmax=maxxi;}
    
   double &val(int i,int j)
    { return btm(i,j);
    }   

   // typage AH 02 2001 
   void maj(int i,int j,double val)
    { //raus((val<0&&val>255)," Image[maj] val  non valide\n");
        btm(i,j)=val; bitmax=(bitmax<val)? val : bitmax;
    }

   // typage AH 02 2001 
   void inc(int i,int j,double val)
    { raus((val<0.0)," Image[inc] val  non valide\n");
      btm(i,j)+=val;
//      if(val>0.0 ) cout <<"IMAGE val on nulle\n";
      if(btm(i,j)>bitmax) bitmax=btm(i,j);
      // cout<<"Image[inc] i= "<<i<<" - j = "<<j<<": val ="<<val<<endl;
      }

   // typage AH 02 2001 
    void raz(int val=0.0)
     { btm.maj(val);bitmax=val;}

   // typage AH 02 2001 
    void charge(char *name=(char*)" ")
      { if(!strcmp(name," ")) name=filename;
        cout<<"Image[charge] Pas encore implemente !\n";
      }

   // typage AH 02 2001 
    void     sauve(char *name=(char*)" ")
      {
	// unsigned char bit;
	unsigned int pix ;
	char *ccTmp = new char[1] ;
        if(!strcmp(name," ")) name=filename;
        ofstream fic(name,ios::out);
        raus(!fic,"Image[sauve] Ouverture fichier impossible!"); 
        //cout<<"Image[sauve] ficname = "<<name<<" - bitmax = "<<bitmax<<endl;  
        fic<<"P5\n";
        fic<<"# MONTE CARLO - MC - 08/1994 \n#      BitMax = "<<bitmax<<endl;
        fic<<btm.maxi()[0]<<" "<<btm.maxi()[1]<<endl;
        fic<<"255\n";
        register int i,j;
        for(j=0;j<btm.maxi()[1];j++)
          for(i=0;i<btm.maxi()[0];i++) {
	    pix= (unsigned char) (btm(i,j)/bitmax*255);
	    pix = (pix>=255)?255:pix;
	    sprintf(ccTmp,"%ud", pix) ; 
	    fic.write(ccTmp,1);
	  }      
        fic.close();
	delete [] ccTmp ;
	//cout<<"Image[sauve] FIN\n";
      }//sauve
  int taille(int i)
  { return btm.maxi()[i];}
  void free() {btm.free();}
};// Image PPM

class RGB {
 private:
   Tabdyn<double,3> btm;
   double bitmax;
   char filename[80];
 public:
    RGB(){}
    RGB(int x,int y,char * nom="out.ppm")
      {init(x, y,nom);}
   // typage AH 02 2001 
    void     init(int x,int y,char * nom="out.ppm")
      { btm.alloue(x,y,3);strcpy(filename,nom);bitmax=-9.9e10;}
    RGB & operator =(RGB &img)
      { cout<<"RGB[=]  DEBUT \n";cout.flush();
        RGB  *pim= new RGB(img.btm.maxi()[0],img.btm.maxi()[1],img.filename);
/*        strcpy(filename,img.filename); cout<<"RGB[=] :"<<filename<<endl;
        btm.alloue(img.btm.maxi()[0],img.btm.maxi()[1]);
*/
    cout<<"RGB[=]  avant copie des variables de btm \n";cout.flush(); 
       for(register int i=0;i<pim->btm.maxi()[0];i++)
         for(register int j=0;j<pim->btm.maxi()[1];j++)
	   for(register int k=0;k<3;k++)
          pim-> btm(i,j,k)=   img.btm(i,j,k);
        return *pim;
      } 
   double maxval(int k)
    { double val;
      register int i,j;
    
      bitmax=-9.9e20;
      for(j=0;j<btm.maxi()[1];j++)
        for(i=0;i<btm.maxi()[0];i++)
	  { val= btm(i,j,k);
           bitmax=(bitmax<val)? val : bitmax;
	  }
      return bitmax;  
    }      
  void fixmax(double maxxi) {bitmax=maxxi;}
   double &val(int i,int j,int k)
    { return btm(i,j,k);
    }    

   // typage AH 02 2001 
    void    maj(int i,int j,int k,double val)
    { //raus((val<0&&val>255)," RGB[maj] val  non valide\n");
        btm(i,j,k)=val; bitmax=(bitmax<val)? val : bitmax;
    }

   // typage AH 02 2001 
    void   maj(int i,int j,double* tval)
    { //raus((val<0&&val>255)," RGB[maj] val  non valide\n");
        btm(i,j,0)=tval[0]; bitmax=(bitmax<tval[0])? tval[0] : bitmax;
	btm(i,j,1)=tval[1]; bitmax=(bitmax<tval[1])? tval[1] : bitmax;
	btm(i,j,2)=tval[2]; bitmax=(bitmax<tval[2])? tval[2] : bitmax;
    }
  
   // typage AH 02 2001 
    void    inc(int i,int j,int k,double val)
    { raus((val<0.0)," RGB[inc] val  non valide\n");
      btm(i,j,k)+=val;
//      if(val>0.0 ) cout <<"RGB val on nulle\n";
      if(btm(i,j)>bitmax) bitmax=btm(i,j);
      // cout<<"RGB[inc] i= "<<i<<" - j = "<<j<<": val ="<<val<<endl;
      }

   // typage AH 02 2001 
    void     raz(double val=0.0)
     { btm.maj(val);bitmax=val;}

   // typage AH 02 2001 
    void     charge(char *name=" ")
      { if(!strcmp(name," ")) name=filename;
        cout<<"RGB[charge] Pas encore implemente !\n";
      }

   // typage AH 02 2001 
    void     sauve(char *name=" ")
      { 
	//unsigned char bit;
	unsigned int pix ;
	char *ccTmp = new char[1] ;

        if(!strcmp(name," ")) name=filename;
        ofstream fic(name,ios::out);
        raus(!fic,"RGB[sauve] Ouverture fichier impossible!"); 
        cout<<"RGB[sauve] ficname = "<<name<<" - bitmax = "<<bitmax<<endl;  
        fic<<"P6\n";
        fic<<"# class RGB - MC - 1996 \n#      BitMax = "<<bitmax<<endl;
        fic<<btm.maxi()[0]<<" "<<btm.maxi()[1]<<endl;
        fic<<"255\n";
        register int i,j,k;
        for(j=0;j<btm.maxi()[1];j++)
          for(i=0;i<btm.maxi()[0];i++)
	    for(k=0;k<3;k++) {
	      /*
		bit= (unsigned char) (btm(i,j,k)/bitmax*255);
		bit=(bit>=255)?255:bit;	       
		// printf("btm(%d,%d,%d)=%d - %f\n",i,j,k,bit,btm(i,j,k));
		// fic.write(&bit,1);  //sizeof(int));
		const char *ccTmp = &bit ; // HA 11 2003
		fic.write(ccTmp,1);  
		/*/
	      pix= (unsigned int) (btm(i,j,k)/bitmax*255);
	      pix=(pix>=255)?255:pix;
	      sprintf(ccTmp,"%ud", pix) ; 
	      fic.write(ccTmp,1);  
	      /**/
             }      
        fic.close();
	cout<<"RGB[sauve] FIN\n";
	delete[] ccTmp ;
      }//sauve

    int taille(int i)
      { return btm.maxi()[i];}

    void free() {btm.free();}
    
};// RGB PPM

#endif

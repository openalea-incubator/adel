
#include <fstream>
using namespace std ;

#include <outils.h>
#include <arbo.h>
//raus() : si cond vraie alors affiche msg et ciao
void raus(bool cond, const char *msg)
 { if(cond)
   { Ferr<<msg<< '\n' ;
      //Ferr.flush(); // '\n' --> endl l'a déjà fait
      exit(22);
    }
 }//raus()

//itoa primitif
// alloc de la chaine fait par la fonction
char * itoa(int x)
 { int lg;
   char * ch;
   int i;
   int puis;
  
   if(x!=0)
   { lg=((int) log10(double(x)));    // cast double(x) pour gcc-3.2.2 
     ch=new char[lg+2];
     if(lg==0) puis=1;
     else puis=(int)pow(10.0,(double)lg);
     for(i=0;i<=lg ;i++)
      { ch[i]=48+(x/puis);
        x-=(x/puis)*puis;
        puis/=10;
      }
   }
   else 
   { lg=-1;
     ch=new char[2];
     ch[0]=48;
   } 
  ch[lg+2]=0;
  return ch;     
 }//itoa()









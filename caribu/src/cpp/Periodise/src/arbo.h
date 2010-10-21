/* permet de referencer des chemins absolu a l'installation du code
   via l'option -DROOT=/racine du compilo
   */
#ifndef _ARBO
#define _ARBO
#define DQ(A) #A
#define WEG(A,B) DQ(A/B)

/* utilisation dans le code eg 
   fopen(chemin(ROOT,fic.m),"r");
   */
#endif

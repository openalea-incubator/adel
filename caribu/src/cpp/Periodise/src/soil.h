/* soil.h : header du calcul de reflectance de sol
           selon le modele d'Hapke, modifie Jaquemoud

*/
#ifdef SOIL_
#define EXTR
#else
#define EXTR extern
#endif

/* initialise les parametre du modele d'Hapke */
EXTR void param_sol(double w_, double h_, double b_, double c_, double bb_, double cc_);

/* calcule la reflectance bidirectionnelle */
EXTR double refbd(double ts, double to, double psi);

/* calcule la reflectance directionnelle-hemispherique */
EXTR double refdh(double ts);

/* calcule la reflectance bihemispherique par integration num (Gauss)*/
EXTR  double refbh();

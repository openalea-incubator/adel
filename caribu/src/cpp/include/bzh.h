/* definit le repertoire ou sont stockes diag.bzh et nzero.bzh */

/* Definitions exportees dans RayCanTools/RayTools/include/system.h
 */


extern void hdmat_init(char*,char*);
extern void hdmat_majname(char*,char*);
extern void hd_calc_Bfar(VEC *,char *,Diffuseur **,double);

/* exportation depuis bzh.cpp */
extern char pcNzName[];
extern char pcDgName[];
extern char pcBfName[];

/** efface les fichiers de donnees persistantes de Canestra
* definition dans bzh.cpp
* fonction appele depuis Canestra si les donnes ne doivent pas 
* etre remanentes */
void EffaceMatrices(void) ;

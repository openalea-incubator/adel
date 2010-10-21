/** \file utils.c
 * \brief utilitaires
 *
 * On encapsule en améliorant une primitive C
 */

#include <cmath>
#include <cstring>
using namespace std ;

#include <utils.h>


// //////// / / / / / / / / CreateOutName / / / / / / / / / ///////////
void CreateOutName (char *pcDest, char* pcSource, char *ext) {
  /* création du nom du fichier de sortie 
   * recherche de la dernière occurence of .png (pcInName could be 
   * ./files.png/file.png
   * so we must get first the generic filename */
   /* Images are created as png exclusively for now */

  char *pcDotPlace = NULL,
  *pcNewDotPlace = NULL;

  strcpy ( pcDest, pcSource ) ;
  pcDotPlace = pcDest ;

  pcNewDotPlace = (char*)strrchr( pcDotPlace,   PATH_SEP) ;
  if (  pcNewDotPlace != NULL) {
    pcDotPlace = pcNewDotPlace ;
  }
  /* pcDotPlace points the begining of the generic name */  

  pcNewDotPlace = (char*)strrchr( pcDotPlace, '.') ;
  if (pcNewDotPlace !=  NULL ) {
    pcDotPlace = pcNewDotPlace ;
  } else {
    pcDotPlace = pcDest + strlen(pcDest) ;
  }
  strcpy (pcDotPlace, ext )  ; 
}
/* CreateOutName */

#define VERY_SMALL 1E-6
/** By innorance, I rewrote atan2.
 * \param fOpp float cote oppose
 * \param fAdj float cote adjacent
 * \todo use atan2 instead
 */

double dArcTangente( float fOpp, float fAdj) {

/** retour resultat */
double dAngle ;

if ( fabs(fAdj) <= VERY_SMALL ) {
  if(fOpp > VERY_SMALL )
    return M_PI / 2.0 ;
  else if ( fOpp < -VERY_SMALL )
    return -M_PI / 2.0 ;
    else
      return 0.0 ;
}
 dAngle = atan( fOpp / fAdj) ;
   if (fAdj < 0) {
     if (fOpp < 0) {
       dAngle -= M_PI ;
     }else{
       dAngle += M_PI ;
     }
   }
   return dAngle ;
}

// ////////////////////////////////
/** la meme, avec plein de sorties d'erreurs 
 * \param fOpp float cote oppose
 * \param fAdj float cote adjacent
 */

double dArcTangente4Can( float fOpp, float fAdj) {
/** Valeur de retour */
double dAngle ;
 fprintf (stderr,"Warning! dArcTangente4Can working\n");
 dAngle = atan( fOpp / fAdj) ;
   if (fAdj < 0) {
     if (fOpp < 0) {
       dAngle -= M_PI ;
     }else{
       dAngle += M_PI ;
     }
   }
   return dAngle ;
}
/* dArcTangente */

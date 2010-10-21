/** \file utils.h
 * \brief Utilitaires 
 *
 * utilitaires encapsulant (en ameliorant) une primitive C
 */

#include <cantools.h>

#ifndef __UTILS_H__
#define __UTILS_H__

/** Nom du fichier de sortie
 *
 * création du nom du fichier de sortie 
 * recherche de la dernière occurence of .png (pcInName could be 
 * ./files.png/file.png
 * so we must get first the generic filename 
 * Images are created as png exclusively for now.
 * Renvoie dans le parametre 1 la concatenation de 2 et 3
 * \param pcDest char* resultat. Doit pointer une zone de memoire valide
 * \param pcSource char* nom de base
 * \param ext char* extension du nom
 */
void CreateOutName (char *pcDest, char* pcSource, char *ext) ;

/** \brief reecriture stupide de atan2 */
double dArcTangente( float fOpp, float fAdj) ;

/** debugging purpose */ 
double dArcTangente4Can( float fOpp, float fAdj) ;

#endif

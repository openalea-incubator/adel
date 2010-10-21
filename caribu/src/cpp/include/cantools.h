/** \file cantools.h 
 *
 * \brief We define path separators a.s.o
 *
 * We define system wide variables, pathes, etc.
 * \todo Use RayTools' system.h
 */

#ifndef _CANTOOLS_H
#define _CANTOOLS_H

// Format des fichiers de configuration
#define BVIS_CONF_VERSION "0.9"
#define CANVIEW_CONF_VERSION "0.9"

#include <stdio.h>
#include "system.h"

#define NJEMAX 100       /* nbre de couleurs de bvis */

#define D3D 3
#define D4D 4  /*dim tableaux */

/*#define CIEL_CLAIR   0.81f, .92f, .98f, 1.0f */
/*#define CIEL_CLAIR   0.64f, .78f, .85f, 1.0f*/
/*#define CIEL_CLAIR   0.47f, .68f, .80f, 1.0f*/
/*#define CIEL_CLAIR   0.43f, .70f, .82f, 1.0f*/
/*#define CIEL_CLAIR   0.72f, .83f, .90f, 1.0f*/
#define CIEL_CLAIR   0.64f, .80f, .94f, 1.0f

/*#define CIEL_COUVERT  0.15f, 0.15f, 0.35f, 1.0f*/
#define CIEL_COUVERT  0.54f, 0.62f, 0.64f, 1.0f
//#define LIMON           0.09f,  0.075f,  0.072f 
#define LIMON           0.24f,  0.25f,  0.27f  /* biblio: 0.24, 0.25, 0.27*/
//# define SOL		0.49f,	0.39f,	0.2f
# define SOL		0.4f,	0.26f,	0.13f

// passage des valeurs Normalisées vers un octet CHAR
#define N2CHAR 255
// degres * DEG2RAD -> radians
#define DEG2RAD 0.017453
#define RAD2DEG 57.295773

/** The distance beyond which things are considered te be far away */
#define FAR_AWAY 99.

/** A false/true type.
 *
 * There's no STL in C : we define the bool type for the project. */
#ifndef __cplusplus
typedef enum {false, true} bool ;
#endif

/** stores a can node 
 *
 * Contains information to identify the specie which the node belongs,
 * so that we can know wich material (ambiant color and so..) is to be 
 * used to display it. */
typedef struct _CanData {
  int
    iTranspar,
    iEspeceId,
    iPlanteId ,   /* plante Id */
    iCodeOrg ,   /* organe code */
    iOrganeId ,   /*  ID */
    iPolygnId ;   /* Polygone Id */
    
} CanData ;

 /* Image */
void Get_Image(char * ImageBuffer, 
		int iViewWidth, 
	       int iViewHeight) ;
/* Copie le buffer pcImageBuffer, aimenté par bvis, dans ImageBuffer,
 * de taille iViewWidth * iViewWidth * taile_de_pixel_en_byte
 * chaque programme doit fournir sa propre version */


int      Quitter          (int code) ; /* chaque programme sa version */

//#ifndef BCC32
int yylex (void) ;
//#endif

void Decode_Data (CanData * res, double donnee ) ; /*label2int.c */
/* renvoie une structure permettant l'identificaton complète
 * de la primitive */

void Write_Png_File (FILE *out, 
		     int iViewHeight,
		     int iViewWidth ) ;

int Label_2_Index (CanData data) ;
/* crée un index univoque avec le label de la primitive passe dans
 * une structure Prim *
 * appelle Decode_Data */

int Label_2_Indexf (double donnee ) ;
/* crée un index univoque avec le label de la primitive passe au format
 * double
 * appelle Label_2_Index */

#endif

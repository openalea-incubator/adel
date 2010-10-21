
#ifndef __TRANSF_H__
#define __TRANSF_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <sys/types.h>
#ifndef WIN32
#include <sys/ipc.h>
#include <sys/shm.h>
#endif
	  
#ifndef __GNUG__
#include "bool.h"
#endif

//#define CLEF (key_t)1000
#define OFFS 1900
#define CANOPY_SHMID 68 // must be below 100 (-1<id<100)

//#define SEGSIZE 60000 // vautre sur bcgn et ultra2
// #define SEGSIZE 1000 Augmentation le 25/11/98 car plantage Run4 dans $Calg/dev/Canestra/Caribu/Test2/ avec 1025 Leaf
#define SEGSIZE 200000 

typedef struct{
  signed char t; /*type :  0 sol, -i transparent de l'esp. i, i opque de l'esp. i*/
  float P[3][3];
} Patch;

/* These functions are to compute IDs for shared memory segments */
void DecodeClefIn (int *Nt,int *iClefShmOut,int iClefShm); /* canestra */
void DecodeClefOut(int *Nt, int *iClefShmIn,int iClefShm); /* celle de S2V */

#endif

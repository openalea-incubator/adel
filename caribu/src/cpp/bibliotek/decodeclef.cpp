#include <transf.h>

  /* These trancoding codes are the ones picked out of Radioxity.cpp and 
   * s2v.cpp (MC98) so, *every* program amongst the RayTools should be able 
   *to use the same code when turning MC's clef_shm into MC's shared seg IDs */

void DecodeClefIn(int *Nt,int *iClefShmIn,int iClefShm) {
  /* Back from the field */
  
  *Nt=iClefShm/100;
  *iClefShmIn= (iClefShm%100) + OFFS ;
}

void DecodeClefOut(int *Nt, int *iClefShmOut, int iClefShm) { /*  S2V's */
  /* going out to the field */
  *Nt=iClefShm/100;
  *iClefShmOut = iClefShm%100; 
}


// Defining TRUE and FALSE is usually a Bad Idea,
// because you will probably be inconsistent with anyone
// else who had the same clever idea.
// Therefore:  DON'T USE THIS FILE.

#ifndef _SGI

#ifndef _bool_h
#define _bool_h 1

//#include   <_G_config.h>

//enum bool { FALSE = 0, false = 0, TRUE = 1, true = 1 };
#define bool char
#define FALSE 0
#define false 0
#define TRUE 1
#define true 1
#endif
#endif



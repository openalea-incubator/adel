/*  Mersenne Twister random number generator */

/*
init_genrand(seed) initializes the state vector by using
one unsigned 32-bit integer "seed", which may be zero.

init_by_array(init_key, key_length) initializes the state vector 
by using an array init_key[] of unsigned 32-bit integers
of length key_kength. If key_length is smaller than 624,
then each array of 32-bit integers gives distinct initial
state vector. This is useful if you want a larger seed space
than 32-bit word. 
*/

void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);

/*
After initialization, the following type of pseudorandom numbers
are available. 

genrand_int32() generates unsigned 32-bit integers.
genrand_int31() generates unsigned 31-bit integers.
genrand_real1() generates uniform real in [0,1] (32-bit resolution). 
genrand_real2() generates uniform real in [0,1) (32-bit resolution). 
genrand_real3() generates uniform real in (0,1) (32-bit resolution).
genrand_res53() generates uniform real in [0,1) with 53-bit resolution. 
*/

unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

double Moro_NormSInv(double u);


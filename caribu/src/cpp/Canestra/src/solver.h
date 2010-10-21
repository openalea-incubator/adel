#ifdef _solver
#define EXTR
#else
#define EXTR extern 
#endif
EXTR void hd_mgcr(VEC *x,VEC *b, Diffuseur **TabDiff,double tol,int krylov,int limit, int *steps);
EXTR void print_hd_mat(Diffuseur **TabDiff);

#ifdef _NFF
#define EXTR
#else
#define EXTR extern 
#endif

//protos utilise dans Canopy

//partie commune
EXTR void proj_ortho(Diffuseur* E, reel *);
EXTR void init_proj(Diffuseur * diffR );
EXTR void stat_NFF();

//partie differente entre version RAM et HD
#ifdef _HD
EXTR void NFF(int i_sup,int i_inf,VEC **Cfar);
EXTR int ADS_val(int n);
EXTR void ADS_maj(int n,int val);
EXTR void init_NFF(char * envname,double *Esource,int nbpr,int nbf,double &Rsph,bool bias);
#else
EXTR void NFF(SPMAT*FF,int &i_sup,int &i_inf,VEC **Cfar);
EXTR void init_NFF(char * envname,double *Esource,double &Rsph,bool bias);
#endif

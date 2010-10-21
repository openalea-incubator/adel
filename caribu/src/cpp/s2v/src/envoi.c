#include "transf.H"

extern int errno;
shminfo seginf;

main(int argc, char **argv){
  Patch * T,*Ts;
  int n,t;
  double id;
  FILE * fout;
  register int i,j,k;
  char cmd[50];
  int shmid;
  bool prog=false;


  printf("<!>Taille du segment a reserve : %do\n\n",atoi(argv[1])*sizeof(Patch));
  //printf("<!>SHMMAX=%d, shmmax=%d\n\n",(int)SHMMAX, seginf.shmmax);
  T=new Patch[atoi(argv[1])];
  n=atoi(argv[1]);
  for(i=0;i<n;i++){
    T[i].t=(i%2)*2-1;
    printf("envoi: T(%d).t=%d\n",i,T[i].t);
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	T[i].P[j][k]=drand48();
  }//for n T

  //appel de s2v avec  segment partage
  shmid=shmget((key_t)68,SEGSIZE*sizeof(Patch),IPC_CREAT|0666);
  printf("=> shmid = %d\n",shmid);
  if(shmid==-1){
    fprintf(stderr,"Creation du segment partage impossible : erreur no. %d\n",errno);
    return (-1);
  } 
  Ts=(Patch *) shmat(shmid,0,NULL);;
  for(i=0;i<n;i++){
    if(T[i].t==0)
      Ts[i].t=T[i].t;
    else{
      if(T[i].t<0)
	Ts[i].t=T[i].t;
      else
	Ts[i].t=T[i].t;
    }
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	Ts[i].P[j][k] =T[i].P[j][k];
  }
    
  sprintf(cmd,"s2v %d 2 1.0 test.8 nir par ",68+n*100);
  printf("Appel de %s\n %c",cmd,7);
  system(cmd);
    
  shmctl(shmid,IPC_RMID,0);
  shmdt((void*)shmid);
  
  delete [] T;
}//main()


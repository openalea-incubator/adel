/**************************************************************************
**
** Copyright (C) 1993 David E. Stewart & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
====> Adaptation du solver MGCR de Leyk a la matrice specifique des FF
      dans le couvert, stoke sur disque hd_mgcr - MC 1996

***************************************************************************/

#include <iostream> //.h>
using namespace std;

#include <cmath>
#include <cstdio>

#include "canopy.h"
#include "outils.h"

#define _solver
#include "solver.h"
#include "bzh.h"


/*  PROTOTYPES_IN_STRUCT deja defini via canopy.H */
extern "C" {
#ifndef PROTOTYPES_IN_STRUCT
#define PROTOTYPES_IN_STRUCT
#endif
#include "sparse.h"
#include "iter.h"
#include "matrix2.h"
}

static double nff;

/* hd_mv_mlt -- hard disk sparse matrix/dense vector multiply
   -- result is in out, which is returned unless out==NULL on entry
   --  if out==NULL on entry then the result vector is created */
VEC *hd_mv_mlt(SPMAT* A,VEC *x,VEC* out) {
  int	tamp[2],i,is, j_idx,j, n,iff,nd,*diag;//ip : indice prim, is indice face
  Real	sum, *x_ve,sum2;
  double dff;//pour passer d'un FF en pixel a un FF reel
  double rho[2],tau[2],po;
  Diffuseur **TabDiff=(Diffuseur **)A;
  FILE *fic;
  char transp;
  
  fic=fopen(pcDgName,"rb");
  
  fread(&n,sizeof(int),1,fic);
  fread(&nd,sizeof(int),1,fic);
  fread(&dff,sizeof(double),1,fic);
  if ( ! A || ! x )
    error(E_NULL,(char*)"hd_mv_mlt");
  if ( x->dim != n ){
    Ferr <<"in "  << pcDgName<<" A=" << n<<"^2 but x.dim="<<x->dim<<" !\n" ;
    error(E_SIZES,(char*)"hd_mv_mlt");
  }
  if ( ! out || out->dim < n )
    out = v_resize(out,n);
  if ( out == x )
    error(E_INSITU,(char*)"hd_mv_mlt");
  x_ve = x->ve;
  //chargement des indices de la diago
  diag=new int[nd+1];//nd= nb prim + 1
  fread(diag,sizeof(int),nd+1,fic);
  fclose(fic);
  //Produit fait ligne a ligne
  fic=fopen(pcNzName,"rb");
  /*  FILE *fica,*f2;
      fica=fopen("nzero.txt","r");
      f2=fopen("nzero.dbx","w");
      */
  is=0;
  for ( i = 0; i < nd; i++ ){
     
     transp=(diag[i+1]>0)?0:1;
     //printf(".. i=%d, is=%d, diag=%d, n=%d, nd=%d ",i,is,diag[i+1],n,nd); fflush(stdout);
     sum=x_ve[is];//la diago vaut 1
     TabDiff[is]->activ_num(is);
    rho[0]= -TabDiff[is]->rho();
    if(transp) {
      sum2=x_ve[is+1];
      tau[1]=-TabDiff[is]->tau();
      TabDiff[is]->togle_face();
      rho[1]= TabDiff[is]->rho();
      tau[0]=TabDiff[is]->tau();
      TabDiff[is]->togle_face();
    }
    TabDiff[is]->activ_num(is);
    for (j_idx = (int) fabs(double(diag[i]))-1; j_idx<fabs(double(diag[i+1]))-1; j_idx++) {
      fread(tamp,sizeof(int),2,fic);
      j=tamp[0]; iff=tamp[1];
      /*if(j>=2000) 
	Ferr <<"Solver : i="  << i<<" : j="  << j<<" - ji-= "  
	<< diag[i]<<" -ji+="  << diag[i+1]<<"- iff="  << iff<<" \n" ;
	fprintf(f2,"%d %d  ",j,iff);
	fscanf(fica,"%d %d",&j,&iff);      
	if(j>=2000) 
	Ferr <<"Sol ASCII: i="  << i<<" : j="  << j<<" - ji-= "  
	<< diag[i]<<" -ji+="  << diag[i+1]<<"- iff="  << iff<<" \n" ;
       */
      if(iff>0)
	po=rho[0];
      else
	po=tau[0];
      sum += x_ve[j]*iff*dff*po;
      if(transp) {
	if(iff>0)
	  po=tau[1];
	else
	  po=rho[1];
	sum2 += x_ve[j]*iff*dff*po;
      }
    }
    out->ve[is] = sum;
     
    if(transp){
      out->ve[++is] = sum2;
      //Ferr <<"T\n" ;
    }
    //else
    //Ferr <<"O\n" ;

    //fflush(stdout);
      is++;
    // fprintf(f2,"\n");
  }//for i (ligne)
  fclose(fic);
  //fclose(fica); fclose(f2);
  return out;
}


void hd_mgcr(VEC *x,VEC *b, Diffuseur **TabDiff,double tol,int krylov,int limit, int *steps) {
  //mise en forme de la structure iterative
  ITER *ip;
  
  ip = iter_get(0,0);
  ip->Ax = (Fun_Ax) hd_mv_mlt;
  ip->A_par =TabDiff ;
  ip->Bx = (Fun_Ax) NULL;
  ip->B_par = NULL;
  ip->k = krylov;
  ip->limit = limit;
  ip->info = (Fun_info) NULL;
  ip->b = b;
  ip->eps = tol;
  ip->x = x;
   
  //resolution proprement dite
  static VEC *As, *beta, *alpha, *z;
  static MAT *N, *H;
  VEC *rr, v, s;  /* additional pointer and structures */
  Real nres;      /* norm of a residual */
  Real dd;        /* coefficient d_i */
  int i,j;
  int done;      /* if TRUE then stop the iterative process */
  int dim;       /* dimension of the problem */
   
  /* ip cannot be NULL */
  if (ip == INULL) error(E_NULL,(char*)"hd_mgcr");
  /* Ax, b and stopping criterion must be given */
  if (! ip->Ax || ! ip->b || ! ip->stop_crit) 
    error(E_NULL,(char*)"hd_mgcr");
  /* at least one direction vector must exist */
  if ( ip->k <= 0) error(E_BOUNDS,(char*)"mgcr");
  /* if the vector x is given then b and x must have the same dimension */
  if ( ip->x && ip->x->dim != ip->b->dim)
    error(E_SIZES,(char*)"hd_mgcr");
  if (ip->eps <= 0.0) ip->eps = MACHEPS;
   
  dim = ip->b->dim;
  As = v_resize(As,dim);
  alpha = v_resize(alpha,ip->k);
  beta = v_resize(beta,ip->k);
   
  MEM_STAT_REG(As,TYPE_VEC);
  MEM_STAT_REG(alpha,TYPE_VEC);
  MEM_STAT_REG(beta,TYPE_VEC);
   
  H = m_resize(H,ip->k,ip->k);
  N = m_resize(N,ip->k,dim);
   
  MEM_STAT_REG(H,TYPE_MAT);
  MEM_STAT_REG(N,TYPE_MAT);
   
  /* v and s are additional pointers to rows of N */
  /* they must have the same dimension as rows of N */
  v.dim = v.max_dim = s.dim = s.max_dim = dim;
   
   
  done = FALSE;
  for (ip->steps = 0; ip->steps < ip->limit; ) {
    (*ip->Ax)(ip->A_par,ip->x,As);         /* As = A*x */
    v_sub(ip->b,As,As);                    /* As = b - A*x */
    rr = As;                               /* rr is an additional pointer */
      
    /* norm of the residual */
    nres = v_norm2(rr);
    dd = nres;                            /* dd = ||r_i||  */
      
    /* check if the norm of the residual is zero */
    if (ip->steps == 0) {                
      /* information for a user */
      if (ip->info) (*ip->info)(ip,nres,As,rr); 
      ip->init_res = fabs(nres);
    }
    if (nres == 0.0) { 
      /* iterative process is finished */
      done = TRUE; 
      break;
    }
      
    /* save this residual in the first row of N */
    v.ve = N->me[0];
    v_copy(rr,&v);
      
    for (i = 0; i < ip->k && ip->steps < ip->limit; i++) {
      ip->steps++;
      v.ve = N->me[i];                /* pointer to a row of N (=s_i) */
      /* note that we must use here &v, not v */
      (*ip->Ax)(ip->A_par,&v,As); 
      rr = As;                        /* As = A*s_i */
      
      if (i < ip->k - 1) {
	s.ve = N->me[i+1];         /* pointer to a row of N (=s_{i+1}) */
	v_copy(rr,&s);                   /* s_{i+1} = B*A*s_i */
	for (j = 0; j <= i-1; j++) {
	  v.ve = N->me[j+1];      /* pointer to a row of N (=s_{j+1}) */
	  /* beta->ve[j] = in_prod(&v,rr); */      /* beta_{j,i} */
	  /* modified Gram-Schmidt algorithm */
	  beta->ve[j] = in_prod(&v,&s);  	         /* beta_{j,i} */
	  /* s_{i+1} -= beta_{j,i}*s_{j+1} */
	  v_mltadd(&s,&v,- beta->ve[j],&s);    
	}
	    
	/* beta_{i,i} = ||s_{i+1}||_2 */
	beta->ve[i] = nres = v_norm2(&s);     
	if ( nres <= MACHEPS*ip->init_res) { 
	  /* s_{i+1} == 0 */
	  i--;
	  done = TRUE;
	  break;
	}
	sv_mlt(1.0/nres,&s,&s);           /* normalize s_{i+1} */
	    
	v.ve = N->me[0];
	alpha->ve[i] = in_prod(&v,&s);     /* alpha_i = (s_0 , s_{i+1}) */
	    
      }
      else {
	for (j = 0; j <= i-1; j++) {
	  v.ve = N->me[j+1];      /* pointer to a row of N (=s_{j+1}) */
	  beta->ve[j] = in_prod(&v,rr);       /* beta_{j,i} */
	}
	    
	nres = in_prod(rr,rr);                 /* rr = B*A*s_{k-1} */
	for (j = 0; j <= i-1; j++)
	  nres -= beta->ve[j]*beta->ve[j];

	if (sqrt(fabs(nres)) <= MACHEPS*ip->init_res)  {
	  /* s_k is zero */
	  i--;
	  done = TRUE;
	  break;
	}
	if (nres < 0.0) { /* do restart */
	  i--; 
	  ip->steps--;
	  break; 
	}   
	beta->ve[i] = sqrt(nres);         /* beta_{k-1,k-1} */
	    
	v.ve = N->me[0];
	alpha->ve[i] = in_prod(&v,rr); 
	for (j = 0; j <= i-1; j++)
	  alpha->ve[i] -= beta->ve[j]*alpha->ve[j];
	alpha->ve[i] /= beta->ve[i];                /* alpha_{k-1} */
	    
      }
	 
      set_col(H,i,beta);

      /* other method of computing dd */
      /* if (fabs((double)alpha->ve[i]) > dd)  {     
	 nres = - dd*dd + alpha->ve[i]*alpha->ve[i];
	 nres = sqrt((double) nres); 
	 if (ip->info) (*ip->info)(ip,-nres,VNULL,VNULL);  	
	 break;     
	 }  */
      /* to avoid overflow/underflow in computing dd */
      /* dd *= cos(asin((double)(alpha->ve[i]/dd))); */
	 
      nres = alpha->ve[i]/dd;
      if (fabs(nres-1.0) <= MACHEPS*ip->init_res) 
	dd = 0.0;
      else {
	nres = 1.0 - nres*nres;
	if (nres < 0.0) {
	  nres = sqrt((double) -nres); 
	  if (ip->info) (*ip->info)(ip,-dd*nres,VNULL,VNULL);  	
	  break;
	}
	dd *= sqrt((double) nres);  
      }

      if (ip->info) (*ip->info)(ip,dd,VNULL,VNULL);     
      if ( ip->stop_crit(ip,dd,VNULL,VNULL) ) {
	/* stopping criterion is satisfied */
	done = TRUE;
	break;
      }
	 
    } /* end of for */
      
    if (i >= ip->k) i = ip->k - 1;
      
    /* use (i+1) by (i+1) submatrix of H */
    H = m_resize(H,i+1,i+1);
    alpha = v_resize(alpha,i+1);
    Usolve(H,alpha,alpha,0.0);       /* c_i is saved in alpha */
      
    for (j = 0; j <= i; j++) {
      v.ve = N->me[j];
      v_mltadd(ip->x,&v,alpha->ve[j],ip->x);
    }
      
      
    if (done) break;              /* stop the iterative process */
    alpha = v_resize(alpha,ip->k);
    H = m_resize(H,ip->k,ip->k);
      
  }  /* end of while */
   
  /* return ip->x;                return the solution */
  //Extraction des resultats dans ip
  x = ip->x;
  if (steps) *steps = ip->steps;
  ip->shared_x = ip->shared_b = TRUE;
  iter_free(ip);
}//hd_mgcr()

void print_hd_mat(Diffuseur **TabDiff) {
int	i,is, j_idx,j, n,iff,nd,*diag=NULL;//i : indice prim, is indice face
  double dff;//pour passer d'un FF en pixel a un FF reel
  double rho[2],tau[2],po;
  Tabdyn <double,2>ligne;
  FILE *fic;
  char transp;

  Ferr <<"*  print_hd_mat() : Debut\n";
  fic=fopen(pcDgName,"rb");
  Ferr <<"\t-> Lecture de LORZ/diag.bzh\n";
  fread(&n,sizeof(int),1,fic);
  Ferr <<"n = "<<n;  
  ligne.alloue(n,2);
  fread(&nd,sizeof(int),1,fic);
  Ferr <<"\nnd = "<<nd; 
  fread(&dff,sizeof(double),1,fic);
  Ferr <<"\ndff = "<<dff;
  //chargement des indices de la diago
  diag=new int[nd+1];//nd= nb prim + 1
  if(diag==NULL) {
    Ferr <<" Impossible d allouer diag["  << nd+1<<"]!\n" ;
    exit(17);
 }
  fread(diag,sizeof(int),nd+1,fic);
  fclose(fic);
  //Produit fait ligne a ligne
  fic=fopen(pcNzName,"rb");
  Ferr <<"\n\t-> Lecture de LORZ/nzero.bzh\n";
  is=0;
  for ( i = 0; i < nd; i++ ){
    //Ferr <<" ligne "<<i<<endl;
    ligne.maj(0);
    transp=(diag[i+1]>0)?0:1;
    //la diago vaut 1
    ligne(is,0)=1;
    TabDiff[is]->activ_num(is);
    rho[0]= -TabDiff[is]->rho();
    if(transp) {
      tau[1]=-TabDiff[is]->tau();
      TabDiff[is]->togle_face();
      rho[1]= TabDiff[is]->rho();
      tau[0]=TabDiff[is]->tau();
      TabDiff[is]->togle_face();
      ligne(is+1,1)=1;
    }
    TabDiff[is]->activ_num(is);
    for (j_idx = (int) fabs(double(diag[i])); j_idx<fabs(double(diag[i+1])); j_idx++) {
      fread(&j,sizeof(int),1,fic);
      //Ferr <<"j = "<<j<<endl;
      fread(&iff,sizeof(int),1,fic);
      if(iff>0)
	po=rho[0];
      else
	po=tau[0];
      ligne(j,0)=iff*dff*po;
      if(transp) {
	if(iff>0)
	  po=tau[1];
	else
	  po=rho[1];
	ligne(j,1)=iff*dff*po;
      }
    }
    //Ferr <<"n = "<<n<<endl;
    for(j=0;j<n;j++) 
      Ferr  << ligne(j,0)<<" " ;
    Ferr <<"\n" ;
    if(transp) {
      for(j=0;j<n;j++)
	Ferr  << ligne(j,1)<<" " ;
      Ferr <<"\n" ;
      is++;
    }
    is++;
  }//for i (ligne)
  fclose(fic);
  Ferr <<"*  print_hd_mat() : Fin\n";
}// print_hd_mat()



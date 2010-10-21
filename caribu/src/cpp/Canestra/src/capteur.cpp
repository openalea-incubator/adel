//           Capteur.C,  me, today

#include "capteur.h"

/* Le fichier .cam contient les parametres de vue ie
 - oeil : pt de vue
 - vers : vecteur de direction du regard
 - vert : vecteur donnant la verticale
 - angle: angle de vue
 - resol: resolution de l'ONU
*/

//*************** Camera::init() *************************************

void Camera::init(char * ficname, char * imgname)
 { char buffer[80];
   ifstream fic(ficname,ios::in);
   raus(!fic,"Camera[init] Ouverture fichier impossible!");    
   Vecteur vers,verti;
   Point P; 
   double angle_vue;
   int res[2];
   Liste<Point> Lpkt;

 fic>>buffer; // cout <<" 1er buffer "<<buffer<<endl;
   fic>> oeil[0]>>oeil[1]>>oeil[2];
 fic>>buffer;//cout <<" 2e  buffer "<<buffer<<endl;
   fic>>vers[0]>>vers[1]>>vers[2];
   fic>>buffer;  
   fic>>verti[0]>>verti[1]>>verti[2];
   fic>>buffer; 
   fic>>angle_vue;
   fic>>buffer; 
   fic>> res[0]>>res[1];
   fic>>buffer; 
   fic>>foc;
  cout<<"parametres de vue OK\n"; 
  cout<<"oeil "<<oeil[0]<<" "<<oeil[1]<<" "<<oeil[2]<<endl;
  cout<<"vers "<<vers[0]<<" "<<vers[1]<<" "<<vers[2]<<endl;
  cout<<"verti "<<verti[0]<<" "<<verti[1]<<" "<<verti[2]<<endl;
  cout<<" angle "<<angle_vue<<" - resolution "<<res[0]<<" "<<res[1]<<endl;
  cout <<"focal "<<foc<<endl;
   fic.close();
// fin de lecture du fichier .cam
   angle_vue*=M_PI/360.0;
   img=new Image(res[0],res[1],imgname);
   img->raz(0); 
   moy=new Image(res[0],res[1],"moy.ppm");
   moy->raz(0); 
   sigma=new Image(res[0],res[1],"sigma.ppm");
   sigma->raz(0); 

   vers.normalise();
   verti=vers.prod_vectoriel(verti);
   verti=verti.prod_vectoriel(vers);
   verti.normalise(); 
   cout<<" up ortho "<<verti[0]<<" "<<verti[1]<<" "<<verti[2]<<endl;

   E = P = oeil+ vers*foc; 
   foc=  foc*tan(angle_vue);
   P=P+verti*foc;
   P=P+vers.prod_vectoriel(verti)*foc;
   Lpkt.ajoute(P);
   taille=foc*=2.0;
   P=P-verti*foc;
   Lpkt.ajoute(P);
   P=P-vers.prod_vectoriel(verti)*foc;   
   Lpkt.ajoute(P);
   P=P+verti*foc;
   Lpkt.ajoute(P);  
   prim= new Polygone(Lpkt,2436);
   Primitive& primi=*prim;     
   primi.show("Camera[init] ");
   u.formation_vecteur(primi[1],primi[2]);
   u.normalise();u.show();
   v.formation_vecteur(primi[1],primi[0]);
   v.normalise(); v.show();
 }// init

//*************** Camera::capte() ***********************************

void Camera::capte(Rayon &strahl,Diffuseur *pinter,int nbray)
 {  Vecteur dir;
    Param_Inter oldpar;
    Diffuseur *pdif=NULL;
    //   cout<<"Camera[capte]";
    oldpar=strahl.par();
    // oldpar.direct().show();
    // cout<<"Camera[capte] depart rayon : ori  "; strahl.par().origin().show();
    // cout<<"Camera[capte] arrivee rayon: oeil "; oeil.show();    
    dir.formation_vecteur(strahl.par().origin(),oeil);
    dir.normalise(); 
    if (pinter->isopaque())
      if( dir.prod_scalaire(pinter->primi().normal())<0.0)
       { // cout <<"Camera[capte] Pas visble\n";
	 // pinter->primi().show();
	 // pinter->primi().normal().show();
         strahl.init(oldpar); 
         return; // pas visible 
       }	
      
//    cout<<"Camera[capte] direct rayon normee : ";dir.show() ;    
     strahl.par().change_direction(dir); // FAUX!... Ah le C++
     strahl.init_anamatides();;
    do	 
     { 
       if((pdif=strahl.intersect(pinter))!=NULL)
         { strahl.init(oldpar); 
           return;
         }        
        else
          strahl.promenade();
     }while(strahl.in());
    Point I; 
//    cout<<"Camera[capte]  rayon out : "; strahl.par().direct().show();
//    I=strahl.sortir();

   if( prim->intersect(strahl.par(),&I) < ER_MAEZ)
       { //cout<<"Camera[capte] dans la boite! taille "<<taille;I.show();
	   //oldpar.direct().normalise(); 
         Point P(oldpar.direct());
	 strahl.par().change_origine(P);
         dir.formation_vecteur( prim->operator[](0),I);
//    dir.show();cout <<" x= "<<dir.prod_scalaire(u)<<" - y = "<<-dir.prod_scalaire(v)<<endl;
         double weight=pinter->lux(strahl.par());
         weight*=taille*taille* fabs( strahl.par().direct().prod_scalaire(prim->normal()))/I.dist2(oldpar.origin());

         //maj de l'image
         register int i,j;
         i= img->taille(0)-1-(int)fabs(dir.prod_scalaire(u)/taille*(img->taille(0)-1));
         j=(int)fabs(dir.prod_scalaire(v)/taille*(img->taille(1)-1));
         img->inc(i,j, weight ); 
    } 
       
    strahl.init(oldpar); 
//cout <<"Camera[capte] FIN poids = "<<strahl.par().poid()<<endl;
 } // capte

//void Camera::developpe() {}// developpe

//*************** Camera::capte_visi() ******************************

void Camera::capte_visi(Rayon &strahl,Diffuseur *inter,int nbray)
 {  Vecteur dir;
    Param_Inter oldpar;
    Point I;
   
    oldpar=strahl.par();
    // cout<<"Camera[capte_visi] depart rayon : ori  "; strahl.par().origin().show();
    //    cout<<"Camera[capte_visi] arrivee rayon: oeil "; oeil.show();    
    dir.formation_vecteur(strahl.par().origin(),oeil); 
    dir.normalise();   
    if (inter->isopaque())
      if( dir.prod_scalaire(inter->primi().normal())<0.0)
       { //cout <<"Camera[capte_visi] Pas visble\n";
         return; // pas visible 
       }	
      
    //    cout<<"Camera[capte_visi] direct rayon normee : ";dir.show() ;    
    strahl.par().change_direction(dir); 
    strahl.init_anamatides(); 
    //  cout<<"Camera[capte_visi] rayon out : "; strahl.par().direct().show();
    //    I=strahl.sortir();
    if( prim->intersect(strahl.par(),&I) < ER_MAEZ)
       { //cout<<"Camera[capte_visi] dans la boite! taille "<<taille;I.show();
         Point P(oldpar.direct());
         strahl.par().change_origine(P);
         dir.formation_vecteur( prim->operator[](0),I);
         //    dir.show();cout <<" x= "<<dir.prod_scalaire(u)<<" - y = "<<-dir.prod_scalaire(v)<<endl;
        
         // calcul de la valeur du pixel  
         double weight=inter->lux(strahl.par()); 
         weight*=taille*taille* fabs(strahl.par().direct().prod_scalaire(prim->normal()))/I.dist2(oldpar.origin());   

         //maj de l'image
         register int i,j;
         i= img->taille(0)-1-(int)fabs(dir.prod_scalaire(u)/taille*(img->taille(0)-1));
         j= (int)fabs(dir.prod_scalaire(v)/taille*(img->taille(1)-1));
         img->inc(i,j, weight ); 
       } 
    strahl.init(oldpar); 
    //cout <<"Camera[capte] FIN poids = "<<strahl.par().poid()<<endl;
 } // capte_visi()

//*************** Camera::colorie_triangle() ************************
#define A Pp[i]
#define B Pp[j]
#define C Pp[k]
#define D Pp[3]

void Camera::colorie_triangle(Diffuseur * pdif,Point &a,Point &b,Point &c, Tabdyn<double, 2>& Zbuf, Tabdyn<Diffuseur *, 2>& Zprim,Image *pimg)
// Vorsicht : Zbuf, Zprim, pdiff passes en variable globale (-lourd)
 { double penteL,penteR,xL,xR,zL,yrel,vab,vac,vbc,z,dz; // L Left, R  Right
   register int i,j, iR,iL, deby,finy,sens;
   Diffuseur *exdif;
   double K[2];

//            a           b ______ c 
//           /\             \    /
//          /  \      ou     \  /
//        b/____\c            \/ a
/*                                    */
   K[0]=(img->taille(0)-1)/taille;
   K[1]=(img->taille(1)-1)/taille;
   sens   = (b[1] > a[1])? 1 : -1;         
   penteL = (b[0] - a[0]) / (b[1] - a[1]);   
   penteR = (c[0] - a[0]) / (c[1] - a[1]);  

                                         
   if(sens>0)
     { deby   = int(a[1]*K[1]+ 0.5);
       finy   = int(b[1]*K[1]- 0.5);
     }
    else
     { deby   = int(b[1]*K[1]+ 0.5);
       finy   = int(a[1]*K[1]- 0.5);
     }
 
   if( finy<0 || deby>=img->taille(1) || deby>finy )
      return;  
   
   deby=T_max(0,deby);
   finy=T_min(img->taille(1)-1,finy);

//	calcul de grandeurs fixes utilisees pour obtenir l'altitude des pixels
   vab = (b[2] - a[2])/(b[1] - a[1]);
   vbc = (c[2] - b[2])/(b[1] - a[1]);

//	calcul de grandeurs qui seront incrementes dans la boucle sur les lignes
   yrel=((((double)deby+0.5)/K[1])-a[1]);
   xL=penteL*yrel+a[0];
   xR=penteR*yrel+a[0];
   zL= a[2] + yrel*vab;

   penteL/=K[1];
   penteR/=K[1];
   vab/=K[1]; 
   K[1]=-sens/K[1];
   cout<<" Camera[colorie_triangle] : "<<pdif->primi().name()<<endl;
   cout<<"\t A : ";a.show();cout<<"\t B : "; b.show();cout<<"\t C : "; c.show();
   cout<<"\t deby = "<<deby  <<" - finy = "<<finy<<" - yrel = "<<yrel<<endl;
   cout<<"\t xL   = "<<xL    <<" - xR   = "<<xR<<" -zL="<<zL<<endl;
   cout<<"\tpenteL= "<<penteL<<" -penteR= "<<penteR<<" - sens = "<<sens<<" - vab= "<<vab<<endl;

//	boucle sur les lignes 
   for(j=deby;j<=finy;j++)
    { iL=T_min(img->taille(0)-1,T_max(0,int(xL*K[0]+0.5)));    // centre de gravite 
      iR=T_min(img->taille(0)-1,T_max(0,int(xR*K[0]-0.5)));    // dans le triangle
//         cout<<"LOOP\t xL= "<<xL    <<" - xR= "<<xR<< "(iL,iR)="<<iL<<", "<<iR<<" -j= "<<j<<" -zL="<<zL<< endl;
      if (iL<=iR) 
       { dz= yrel*vbc/(xR-xL);
	 z = zL;
	 for(i=iL;i<=iR;i++)
	  { //cout<<"Capteur[colorie-triangle] (i,j) = "<<i<<", "<<j<<endl;
  //cout <<"Camera [colorie_triangle](i,j) ("<<i<<","<<j<<") z = "<<z<<"- Z = "<<Zbuf(i,j)<<endl;
            if(z<Zbuf(i,j))
	     {
		 //    cout <<"Camera [colorie_triangle](i,j) ("<<i<<","<<j<<") z = "<<z<<"- Z = "<<Zbuf(i,j)<<endl;
               exdif=Zprim(i,j);
	       if(exdif!=NULL) {(*exdif)--; }
               (*pdif)++;
	       Zprim(i,j)=pdif; 
               Zbuf(i,j)=z;
               if(pimg!=NULL)pimg->maj(img->taille(0)-1-i,img->taille(1)-1-j,pdif->primi().name());
	     } 
	    z = z+dz;
	  }//for i
       }//if ligne pleine
      yrel+=K[1];
      xL+=penteL;  
      xR+=penteR;
      zL+=vab;
    }// for j	
//   cout<<"FIN\t i = "<<i  <<" - j = "<<j<<" - yrel = "<<yrel<<endl;
//   cout<<"\t xL   = "<<xL    <<" - xR   = "<<xR<<endl;
  }// Camera::colorie_triangle()

//*************** Camera::calc_visi() *******************************

void Camera::calc_visi(Liste <Diffuseur*>& Ldiff)
 { Tabdyn<double, 2> Zbuf(img->taille(0),img->taille(1));
   Tabdyn<Diffuseur *, 2> Zprim(img->taille(0),img->taille(1));
   Point P[3],Pp[4];
   Liste<Point> lpt; Point Pt1,Pt2; ; Primitive *pprim;
   Param_Inter parag;
   Diffuseur *pdiff;
   register int i,j,k,l;
   double distZ,pente;
   Vecteur &w=prim->normal();cout<<" w=vers? ";w.show();
   Vecteur SvE=E,EvI; // SvE : Scene vers Ecran, EvI : Ecran vers Image
   bool up,down,pastoutvu;

   Image imgvisi(img->taille(0),img->taille(1), "visi.ppm");
   imgvisi.raz();
/* Produit vecteur -matrice :pas implemente cf IF 94
   Matrice<double> SvE;

   for(i=0;<3;i++)
      { SvE(i,0)=u[i];
        SvE(i,1)=v[i];
        SvE(i,2)=w[i];
      }
*/
   Zbuf.maj(9.9E20);
   Zprim.maj(NULL);
   EvI.formation_vecteur(E,(*prim)[1]);
   EvI=EvI.chgt_base(u,v,w);
   //cout<<"Camera[calc_visi]EvI doit z=0";EvI.show();
   for(Ldiff.debut();! Ldiff.finito();Ldiff.suivant())
     { pdiff=Ldiff.contenu();
/*     if(pdiff->primi().name()==0) //maj visi sol
         { (*pdiff)++; (*pdiff)++; (*pdiff)--;
            break;
         }      
*/   if(pdiff->primi().name()==0) pdiff->primi().show("Camera[calc_visi] ");
       pastoutvu=false;
       //cout <<"Camera[calc_visi] primitive = "<<pdiff->primi().name()<<endl;
       for(i=0;i<3;i++) // Cas des triangles
	 { P[i]=pdiff->primi()[i];
           P[i]-=SvE;
           P[i]=P[i].chgt_base(u,v,w);
	   cout <<"Camera[calc_visi] P{Re} = ";P[i].show();
           if( (P[i][2]<0.0) || pastoutvu)  // Prim  PARTIELLEMENT pas vue
	       if( P[i][2]<0.0) //
		   {  pastoutvu=true; 
                      //cout <<"Camera[calc_visi] Z neg\n";
		   }
               else
		   { (*pdiff)++; (*pdiff)++; (*pdiff)--;
                     //cout <<"Camera[calc_visi] P vu mais pastoutvu actif!\n";
		   }
             else 
	      { Pp[i][0]=P[i][0]*foc/( P[i][2]+foc);
                Pp[i][1]=P[i][1]*foc/( P[i][2]+foc);
	        Pp[i][2]=0.0;//P[i][2];
	        //cout <<"Camera[calc_visi] Pp{Re}["<<i<<"]  = ";Pp[i].show();
                // Chgt d'origine E->O (coinBG)
               // distZ=Pp[i].dist2(P[i]);
                Pp[i]-=EvI;
	        Pp[i][2]=P[i][2];//distZ;
	        cout <<"Camera[calc_visi] Pp{Rimage}["<<i<<"]  = ";Pp[i].show();
              }//else P[i][2]<0.0) || pastoutvu)
        }// for points triangle
       if(!pastoutvu) 
        {// Tri sommets tq Pp[i][1]<<Pp[j][1]<<Pp[k][1] ie A[1] < B[1] < C[1]
         j = (Pp[1][1]>Pp[2][1])? 1: 2; // calc intermed
         k = (Pp[0][1]>Pp[j][1])? 0: j; // indice max pour coord y
         i = (k+1)%3; j= (i+1)%3;
         i = (Pp[i][1]<Pp[j][1])? i : j; // indice min pour coord y
         j = 3- i-k;
//       cout<<" (i,j,k) = "<<i<<j<<k<<endl;
         if ((i!=k)&& !((A[0]==B[0])&&(B[0]==C[0]))&& !((A[1]==B[1])&&(B[1]==C[1])))
         { // Pts A,B,Cpas  alignes selon les axes Xou Y
           // tri Ok
           up=down=false;
            if (A[1]==B[1]) // up 
	       { i = (A[0] <B[0])?i:j;
	         j = 3-i-k; 
                 up=(A[0]==B[0])?false: true;
	       }//if up
	      else
               { if(B[1] == C[1]) // down 
	          { k = (B[0] <C[0])?k:j; 
	            j = 3-i-k;
	            down =(B[0]==C[0])?false: true;
                  }//if down
                 else 
	           {
                     D[1] = B[1];
                     pente=(D[1]-A[1])/(C[1]-A[1]);
	             D[0] = pente*(C[0]-A[0])+A[0];
                   // D[2] = (A[2]*(C[1]-D[1]) + C[2]*(D[1]-A[1]))/(C[1]-A[1]);
	           // D[2] = pente*(C[2]-A[2])+A[2];
                   D[2] = 0.0;
                 // calcul du vrai Z de D
                   lpt.ajoute(P[i]);
                   lpt.ajoute(P[j]);
                   lpt.ajoute(P[k]);
                   pprim= new Triangle(lpt,-1);
                   lpt.free_liste();
                   Pt1=oeil;
                   Pt1-=SvE;
                   Pt1=Pt1.chgt_base(u,v,w);
                   Pt2=D+EvI;
                   Vecteur dir(Pt1,Pt2);
                   dir.normalise();
                   parag.change_origine(Pt1);
                   parag.change_direction(dir);
                  cout <<"param k de l'intersection = "<< pprim->intersect(parag,&Pt1) <<endl;
                   delete pprim;
                   cout<<" D = ";Pt1.show();
                   D[2]=Pt1[2];
                    up=down=(B[0]==D[0])?false: true;
                    if(D[0]>B[0]) { l=i; i=j; j=3;}
	              else        { l=i; i=3;     }
                   }//else cas quelconque, ni up , ni down 
	      }
     //     cout<<"2 (i,j,k) = "<<i<<j<<k<<endl;
            if(up)
	       colorie_triangle(pdiff,C,A,B, Zbuf,Zprim,&imgvisi);
            if(down)
	     { if (up){ k=j; j=i; i=l; }
	       colorie_triangle(pdiff,A,B,C, Zbuf,Zprim,&imgvisi);
             }// if down
        }//if pas un triangle plat
       }//if !pastoutvu 
     }// for liste diffuseurs
    imgvisi.sauve();
 }//calc_visi() 


//*************** Camera::calc_stat() *******************************

void Camera::maj_stat()
 { double x;

   for(register int i=0;i<img->taille(0);i++)
     for(register int j=0;j<img->taille(1);j++)
         { x=img->val(i,j);
           moy->inc(i,j,x);
           x*=x;
           sigma->inc(i,j,x);
           img->maj(i,j,0.0);
         }           
    
 }//Camera::maj_stat()


//*************** Camera::calc_stat() *******************************

void Camera::calc_stat(unsigned int niter=0)
 { raus(niter==0,"Capteur[calc_stat] Doit avoir le nb d'iteration en parametre!");
 // calcul de la moyenne et de la variance
   double x,test1=0, test2=0;
   for(register int i=0;i<img->taille(0);i++)
     for(register int j=0;j<img->taille(1);j++)
         { x=moy->val(i,j)/(double)niter;
           moy->maj(i,j,x);
           img->maj(i,j,x);
           sigma->maj(i,j,(sigma->val(i,j)/(double)niter-x*x)/((double)niter-1.0));
         }           
   img->maxval();
   sigma->maxval();   
 }//Camera::calc_stat()
 











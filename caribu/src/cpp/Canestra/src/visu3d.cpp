/************************************************************
*                 visu3d.C  - mike 95                       *
*       main() du plaquage de donnees sur une structure 3d  *
*************************************************************/

#include <unistd.h>

#include "canopy.h"
//#include "rayon.h"
//#include "lumiere.h"
#include "outils.h"
#include "image.h"
 
#include <GetOpt.h>

#include "chrono.h"

inline void beep(const char *msg="M'enfin ...",int nbeep=1)
{ cout<<(char) 7 <<msg<<endl;
  for(register int i=1;i<1;i++) cout<<(char) 7<<endl;
}//beep()

#define PAUSE(msg)  printf(msg);printf("- Taper la touche Any");getchar();

char verbose=0;

void FileOk(char * nfic){
  if(access(nfic, R_OK)!=0){
    fprintf(stderr,"<!> Arret: fichier %s pas accessible !\n",nfic);
    exit(-1);
  }
}//FileOk()

int main(int argc,char **argv){
  register unsigned int i,j;
  // Gestion des options
  bool sol=false,infty=false, abs=false;
  double max=-1;
  char *imgname=NULL,*maqname=NULL,*envname=NULL,*lightname="std.light",*name8=NULL,dataname[80],*optname=NULL;
  GetOpt option(argc,argv,"hsa8:M:o:l:m:"); 
  int c;
  while((c=option())!=EOF)
    switch(c) {
    case 's' : sol=true;                        break;// ajoute un sol
    case '8' : infty=true; name8=option.optarg;   FileOk(name8); break;// infini  
    case 'M' : max=atof(option.optarg);         break;//  
    case 'o' : imgname=option.optarg;           break;
    case 'm' : maqname=option.optarg;           break;//maquette .can
    case 'l' : lightname=option.optarg;         break;   
    case 'a' : abs=true;                        break;//cas de flux par feuilles (2 faces)  
    case 'h' : Ferr <<"usage : "<<argv[0]<<" [hsa8<name>o<name>m<name>l<name>]\n" ; return;
    default  : Ferr <<"usage : "<<argv[0]<<" [hsa8<name>o<name>m<name>c<name>l<name>]\n"; return ;
    }// switch
  if(imgname==NULL) 
    { imgname= new char[10]; strcpy(imgname,"out.ppm");}      
  if(maqname==NULL) 
    { maqname= new char[10]; strcpy(maqname,"test.can");}      
   if(optname==NULL) 
   { optname= new char[10]; strcpy(optname,"fir.opt");}
   cout <<"\n Fichier maquette  :: "<<maqname;
   cout <<"\n Fichier image ppm :: "<<imgname;
   cout <<"\n Fichier sources   :: "<<lightname; 
   if(sol)
     cout<<   " Avec Sol          :: oui\n"; 

   //Test d'acces au fichier
   FileOk(maqname);
   FileOk(optname);

     


 //Fin Gestion des Options

  //PAUSE("fin des options");
  char * nsolem=NULL;
  reel bornemin[3]={99999999.0,99999999.0,99999999.0};
  reel bornemax[3]={-99999999.0,-99999999.0,-99999999.0};
  Chrono chrono;	
  Canopy scene;
  Diffuseur **TabDiff;
  //PAUSE("fin des declarations");
  
  //Chargement de la scene
  cout<<"Avant Chrono.start\n";
  chrono.Start();
   cout<<"Apres  Chrono.start et avant parse_can() \n";
  scene.parse_can(maqname,optname,name8,(reel *)bornemin, (reel *)bornemax,sol,nsolem,TabDiff);
  cout <<"-> scene.radim = "<<scene.radim<<endl;
  chrono.Stop();
  //PAUSE("charge scene");
  for(j=0; j<3; j++)
    cout<<" Bornemin["<<j<<"] = "<<bornemin[j]<<"  Bornemax["<<j<<"] = "<<bornemax[j]<<endl;

  //    calcul de visibilite (purely geometric)
  Vecteur visee;
  int col,tx,ty,nbp=0;
  double **data,fond[3];
  char nomR[200], nomG[200], nomB[200];
  FILE *fR,*fB,*fG,*fic;
  if(!abs){
    printf(" Valeur par feuille, et non par face (1/0) :");
    scanf("%d",&col);
    if(col!=0) abs=true;
  }
  printf(" mode couleur  (1/0) :");
  scanf("%d",&col);
  printf("\nDirection de visee");
  scanf("%f %f %f",&(visee[0]),&(visee[1]),&(visee[2]));
  printf("\nTaille image en x = ");
  scanf("%d",&tx);
  printf("\nTaille image en y = ");
  scanf("%d",&ty);
  printf("\n couleur du fond & fichiers = ");
  if(col) {
    scanf("%lf %lf %lf",&(fond[0]),&(fond[1]),&(fond[2]));
    scanf("%s %s %s",nomR,nomG,nomB);
    printf("fics = %s-%s-%s\n",nomR,nomG,nomB); 
    FileOk(nomR); FileOk(nomG); FileOk(nomB); 
  }
  else {
    scanf("%lf",&(fond[0]));
    scanf("%s",nomR);
    printf("fic =%s-\n",nomR);
    FileOk(nomR); 
  } 
 
  printf("\n");
  long int **Zno;
  scene.data3d(tx,ty,visee,infty,Zno);
      
  data=new (double*)[scene.radim];
  if (abs) fic=fopen("out.dat","w");  
  fR=fopen(nomR,"r");
  if(col) {
     fG=fopen(nomG,"r");
     fB=fopen(nomB,"r");
  }
  register int ii,im;
  for(i=0;i<scene.radim;){
    im=(TabDiff[nbp]->isopaque())? 1: 2;
    for(ii=0;ii<im;ii++,i++){
      //printf("i=%d, im=%d, ii=%d\n",i,im,ii);
      if(col) {
	data[i]=new double[3];
	if(abs && ii==1){
	  data[i][0]=data[i-1][0];
	  data[i][1]=data[i-1][1];
	  data[i][2]=data[i-1][2]; 
	}
	else{
	  fscanf(fR,"%lf",&(data[i][0]));
	  fscanf(fG,"%lf",&(data[i][1]));
	  fscanf(fB,"%lf",&(data[i][2]));
	}
      }//if col
      else {
	data[i]=new double[1];
	if(abs && ii==1){
	  //printf("Cas abs||ii==1 VRAI\n");
	  data[i][0]=data[i-1][0];
	}
	else{
	  //printf("Cas abs||ii==1 FAUX\n");
	  fscanf(fR,"%lf",&(data[i][0]));
	}
	//printf("data(%d) = %.8lf\n",i,data[i][0]);
	if(abs) fprintf(fic,"%.8lf\n",data[i][0]);
      }
    }//for ii
  }//for scene.radim
  if(abs) fclose(fic);
  if(col) {
    RGB imgD (tx,ty,imgname);
    for(i=0;i<tx;i++)
      for(j=0;j<ty;j++) {
	if(Zno[i][j]==-1)
	  imgD.maj(tx-1-i,ty-1-j,fond);
	else
	  imgD.maj(tx-1-i,ty-1-j,data[(Zno[i][j])]);
      }
    if(max!=-1)
      imgD.fixmax(max);
    imgD.sauve();
    imgD.free();
    fclose(fR); fclose(fG); fclose(fB);
  }
  else {
   Image imgD (tx,ty,imgname);
    for(i=0;i<tx;i++)
      for(j=0;j<ty;j++) {
	if(Zno[i][j]==-1)
	  imgD.maj(tx-1-i,ty-1-j,fond[0]);
	else
	  imgD.maj(tx-1-i,ty-1-j,data[(Zno[i][j])][0]);
      }
    if(max!=-1)
      imgD.fixmax(max);
    imgD.sauve();
    imgD.free();
    fclose(fR);
  }  

  return 0;
}












/* Routines for reading a background scene from an input text file and then
   drawing it after each cpfg redraw (called from iGlFinishUp() in irisGL.c). 
   Handles properly translucent objects.
   Cannot instantiate.

   Author: Radomir Mech    August 1995
           Michael Chelle  Nov 1997
   */

#include <iostream>
using namespace std ;
#include <ferrlog.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#ifndef WIN32
#include <unistd.h>
#endif
#include <malloc.h>
#ifdef _SGI
#include <ctype.h>
#endif
#include <ctype.h>
#include "parse_geometry.h"

struct syntax_item {
  int flag;
  char *keyword;
};


#define POLYGON     1
#define MESH        2
#define PRISM       3
#define CONE        4
#define SPHERE      5
#define MATERIAL    6
#define PUSHMATRIX  7
#define POPMATRIX   8
#define TRANSLATE   9
#define ROTATE     10
#define SCALE      11
#define MULTMATRIX 12
#define RECTANGLE  13
#define POLYGONUV  14
#define CYLINDER   15

#define NUM_FLOATS  50
#define PITEM       6    /*x,y,z,Nx,Ny,Nz*/ /* [11] size of one vertex defini dans cpfg:control.h*/
#define NORMAL          1
#define NORMAL_X        3
#define NORMAL_Y        4
#define NORMAL_Z        5

struct syntax_item syn[]= {
  { POLYGON,    "polygon"},
  { POLYGONUV,  "polygonuv"},
  { RECTANGLE,  "rectangle"},
  { MESH,       "mesh"},
  { PRISM,      "prism"},
  { PRISM,      "box"},         /* ala prism */
  { CONE,       "cone"},
  { CYLINDER,   "cylinder"},
  { SPHERE,     "sphere"},
  { MATERIAL,   "material"},
  { PUSHMATRIX, "pushmatrix"},
  { POPMATRIX,  "popmatrix"},
  { TRANSLATE,  "translate"},
  { ROTATE,     "rotate"},
  { SCALE,      "scale"},
  { MULTMATRIX, "multmatrix"},
  { -1, NULL}                 /* must be the last one */
};


/*************************************************************************/

#ifdef OBSO
void SetTurtle(TURTLE *tu){
  float pt[3], ptc[3];
  float rotmat[16];
  int i;
  
  GetMatrix(rotmat);
  
  /* position */
  pt[0]=pt[1]=pt[2] = 0;
  
  Transform3Point(pt, rotmat, ptc);
  for(i=0;i<3;i++)
    tu->position[i] = ptc[i];
  
  /* heading */
  pt[0]=0; pt[1]=1; pt[2]=0;
  
  Transform3Point(pt, rotmat, ptc);
  for(i=0;i<3;i++)
    tu->heading[i] = ptc[i] - tu->position[i];
  
  /* left */
  pt[0]=-1; pt[1]=0; pt[2]=0;
  
  Transform3Point(pt, rotmat, ptc);
  for(i=0;i<3;i++)
    tu->left[i] = ptc[i] - tu->position[i];
  
  /* up */
  pt[0]=0; pt[1]=0; pt[2]=-1;
  
  Transform3Point(pt, rotmat, ptc);
  for(i=0;i<3;i++)
    tu->up[i] = ptc[i] - tu->position[i];
}
#endif


/*************************************************************************/
/* draws triangles */
void OGLPrim::DrawTriangle(float *p1, float *p2, float *p3){
  register short i;
  float u[3],v[3],w[3],Gt[3],surf;
  
  for(i=0;i<3;i++){
    T[nbT][0][i]=p1[i];
    T[nbT][1][i]=p2[i];
    T[nbT][2][i]=p3[i];
    u[i]=p2[i]-p1[i];
    v[i]=p3[i]-p1[i];
  }
  /* compute the isobarycentre of the triangle */
  for(i=0;i<3;i++){
    Gt[i]=(p1[i]+p2[i]+p3[i])/3.;
  }
  //if(1)fprintf(stderr,"DrawTriangle: Gt=%.4f,  G=%f : ", Gt[1], G[1]);
  
  /* compute the area of the triangle through "produit vectoriel" */    
  w[0] = v[2]*u[1] - v[1]*u[2];
  w[1] = v[0]*u[2] - v[2]*u[0];
  w[2] = v[1]*u[0] - v[0]*u[1];
  
  surf=sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2])/2.0;
  if(surf<1e-8){
    //MC05: elimination des triangles de surface nuls
    Ferr<<">>> parse_geometry.cpp:OGLPrim::DrawTriangle() : elimination d'un triangle aberrant..."<<'\n';

  }
  else{
    //MC05: Triangle bon pour le service
    area+=surf;
    for(i=0;i<3;i++) G[i]+=Gt[i]*surf;
    nbT++;
  //if(1)fprintf(stderr," surf=%f,G(%.4f,%.4f,%.4f) => area[%d]=%f\n",
  //				 surf,G[0],G[1],G[2],nbT,area);
  }
}

/*************************************************************************/
/* draws a rectangle with 2 triangles */
void OGLPrim::DrawRectangle(float *p1, float *p2, float *p3, float *p4){
  float vec1[3], vec2[3], norm[3], normc[3];
  float zero[3]={0}, zeroc[3];
  int j;
  float pt[3][PITEM];
  float rotmat[16];
  
  for(j=0; j<3 ; j++) {
    vec1[j] = p3[j] - p1[j];
    vec2[j] = p2[j] - p1[j];
  }
  /* get the normal */
  norm[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  norm[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  norm[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
  
  /* if(dr->output_type == TYPE_RAYSHADE)
    for(j=0;j<3;j++)
      norm[j] = -norm[j];
  */
  GetMatrix(rotmat);
  Transform3Point(norm, rotmat, normc);
  Transform3Point(zero, rotmat, zeroc);
  
  for(j=0;j<3;j++) {
    pt[j][NORMAL_X] = normc[0]-zeroc[0];
    pt[j][NORMAL_Y] = normc[1]-zeroc[1];
    pt[j][NORMAL_Z] = normc[2]-zeroc[2];
  }
  
  Transform3Point(p1, rotmat, pt[0]);
  Transform3Point(p2, rotmat, pt[1]);
  Transform3Point(p3, rotmat, pt[2]);


  DrawTriangle(pt[0], pt[1], pt[2]);

  Transform3Point(p4, rotmat, pt[1]);

  DrawTriangle(pt[0], pt[2], pt[1]);
}

/*************************************************************************/

short OGLPrim::prim2triangle(){
  int i,j,x,y,z,c;
  short answer=0;
  float vec1[3], vec2[3], vec[2], norm[3], normc[3];
  float zero[3]={0}, zeroc[3];
  float rotmat[16];
  float P[2][2][2][3];
  char draw = 1, translucent;
  float rotate90[4]= {90,1,0,0};
  float translate[3] ={0};
  float pt[4][PITEM];        /* for drawing triangles */

#ifdef OBSO
  TURTLE dummy_turtle;       /* for positioning spheres */

  /* set the dummy turtle */ 
  InitializeTurtle(dr,&dummy_turtle);
  dummy_turtle.color_index = clp.colormap;
#endif

  switch(primitive.flag) 
  {
  case MATERIAL:       break;
  case MESH:
    if(primitive.num_items >= 4*3){
      for(j=0; j<= primitive.num_items-4*3; j+=2*3) {
	DrawRectangle(&primitive.data[j] , 
		      &primitive.data[j+6],
		      &primitive.data[j+9], 
		      &primitive.data[j+3] );
      }
      answer=1;
    }
    else
      Ferr << "Caribu<prim2triangle> At least four points must be specified"
	" for mesh." << "\n";
    break;
    
  case POLYGON:
    if(primitive.num_items >= 9) {
      for(j=0; j<3 ; j++) {
	vec1[j] = primitive.data[j] - primitive.data[3+j];
	vec2[j] = primitive.data[6+j] -primitive.data[3+j];
      }
      /* get the normal */
      norm[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
      norm[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
      norm[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
      
      /* if(dr->output_type == TYPE_RAYSHADE)
	for(j=0;j<3;j++)
	  norm[j] = -norm[j];
      */
      GetMatrix(rotmat);
      Transform3Point(norm, rotmat, normc);
      Transform3Point(zero, rotmat, zeroc);
      
      for(j=0;j<3;j++) {
	pt[j][NORMAL_X] = normc[0] - zeroc[0];
	pt[j][NORMAL_Y] = normc[1] - zeroc[1];
	pt[j][NORMAL_Z] = normc[2] - zeroc[2];
      }
      
      Transform3Point(&primitive.data[0], rotmat, pt[0]);
      Transform3Point(&primitive.data[3], rotmat, pt[1]);
      Transform3Point(&primitive.data[6], rotmat, pt[2]);
      
      DrawTriangle(pt[0], pt[1], pt[2]);
      
      for(j=9; j<= primitive.num_items-3; j+=3) { 
	for(c=0;c<3;c++) {
	  pt[1][c] = pt[2][c];
	}
	Transform3Point(&primitive.data[j], rotmat, pt[2]);
	
	DrawTriangle(pt[0], pt[1], pt[2]);
      }
      answer=1;
      /* DEBUG: test de l'effet des matrices stackees */
      vec2[0] = pt[1][0]-pt[0][0];
      //fprintf(stderr,"POLYGON P2x-P1x: anant=%.3f, apres=%.3f\n",vec1[0],vec2[0]);
      
    }
    else
      Ferr <<"POLYGON P2x-P1x: anant="<<vec1[0]<<", apres="<<vec2[0]<<"\n" 
	   <<"Caribu<prim2triangle>Polygon At least three points"
	   <<" must be specified for polygon.\n";
    break;
  case POLYGONUV:
    if(primitive.num_items >= 9*2) {	  
      GetMatrix(rotmat);
      Transform3Point(zero, rotmat, zeroc);
      
      Transform3Point(&primitive.data[0], rotmat, pt[0]);
      Transform3Vector(&primitive.data[3], rotmat, pt[0]+3);
      Transform3Point(&primitive.data[6], rotmat, pt[1]);
      Transform3Vector(&primitive.data[9], rotmat, pt[1]+3);
      Transform3Point(&primitive.data[12], rotmat, pt[2]);
      Transform3Vector(&primitive.data[15], rotmat, pt[2]+3);
      
      /* if(dr->output_type == TYPE_RAYSHADE)
	for(j=0;j<3;j++) {
	  pt[0][3+j] = -pt[0][3+j];
	  pt[1][3+j] = -pt[1][3+j];
	  pt[2][3+j] = -pt[2][3+j];
	}
     */ 
      DrawTriangle(pt[0], pt[1], pt[2]);
      
      for(j=18; j<= primitive.num_items-6; j+=6) { 
	for(c=0;c<6;c++) {
	  pt[1][c] = pt[2][c];
	}
	Transform3Point(&primitive.data[j], rotmat, pt[2]);
	Transform3Point(&primitive.data[j+3], rotmat, pt[2]+3);
	
	/* if(dr->output_type == TYPE_RAYSHADE)
	  for(c=0;c<3;c++)
	    pt[2][3+c] = -pt[2][3+c];
        */
	DrawTriangle(pt[0], pt[1], pt[2]);
      }
    }
    else
      Ferr <<"Caribu<prim2triangle> At least three points with"
	   <<" normals must be specified for polygonuv\n";
  break;
 case RECTANGLE:
   if(primitive.num_items >= 2) {
     vec[0]= vec[1] = 0;
     
     if(primitive.num_items >= 4) {
       vec[0] = primitive.data[2];
       vec[1] = primitive.data[3];
       
       for(j=0;j<4;j++) {
	 pt[j][2] = 0;
       }
       
       for(j=0;j<2;j++) {
	 pt[0][j] = vec[j];
       }
       pt[1][0] = vec[0] + primitive.data[0];
       pt[1][1] = vec[1];
       pt[2][0] = vec[0] + primitive.data[0];	   
       pt[2][1] = vec[1] + primitive.data[1];	   
       pt[3][0] = vec[0];
       pt[3][1] = vec[1] + primitive.data[1];
       
       DrawRectangle(pt[0],pt[1],pt[2],pt[3]);
       
     }
     answer=1;
   }
   else
     Ferr <<"At least two values must be specified for rectangle.\n";
   break;
   
 case PUSHMATRIX:
   PushMatrix();
   if(primitive.num_items >= 1)
     Ferr <<"Caribu<prim2triangle> Warning: pushmatrix doesn't need a parameter.\n" ;
   break;
 case POPMATRIX:
   PopMatrix();
   if(primitive.num_items >= 1)
     Ferr <<"Caribu<prim2triangle> Warning: pushmatrix doesn't need a parameter.\n" ;
   break;
 case TRANSLATE:
   if(primitive.num_items >= 3)
     Translate(primitive.data);
   else
     Ferr <<"Caribu<prim2triangle> Warning: not enough parameters for translate!\n" ;
   break;
 case ROTATE:
   if(primitive.num_items >= 4)
     Rotate(primitive.data);
   else
     Ferr <<"Caribu<prim2triangle> Warning: not enough parameters for rotate!\n" ;
   break;
 case SCALE:
   if(primitive.num_items >= 3)
     Scale(primitive.data);
   else
     Ferr <<"Caribu<prim2triangle> Warning: not enough parameters for scale!\n" ;
   break;
 case MULTMATRIX:
   if(primitive.num_items >= 16)
     MultMatrix(primitive.data);
   else
     Ferr <<"Caribu<prim2triangle> Warning: not enough parameters for scale!\n" ;
   break;
 case SPHERE:
   Ferr <<"Caribu<prim2triangle> sphere not yet implemeted\n" ;
   /* if(primitive.num_items >=1) {
      SetTurtle(&dummy_turtle);
      (*dr->tdd->Sphere)(&dummy_turtle, dr, &viewparam, 
      primitive.data[0]*2);
      }
      else
      Ferr <<"Caribu<prim2triangle> Warning: not enough parameters for sphere!\n" ;
   */
   break;
 case CYLINDER:
   Ferr <<"Caribu<prim2triangle> cylinder not yet implemeted\n" ;
   /*if(primitive.num_items >=2) {
     SetTurtle(&dummy_turtle);
     dummy_turtle.line_width = 2*primitive.data[0];
     (*dr->tdd->StartNode)(&dummy_turtle, dr, &viewparam,
     primitive.data[1],'\0');
     PushMatrix();
     Rotate(rotate90);
     SetTurtle(&dummy_turtle);
     (*dr->tdd->Circle3D)(&dummy_turtle, dr, &viewparam,
     primitive.data[0]*2);
     PopMatrix();
     PushMatrix();
     translate[1] = primitive.data[1];
     Translate(translate);
     SetTurtle(&dummy_turtle);
     (*dr->tdd->EndNode)(&dummy_turtle, dr, &viewparam,'\0');
     Rotate(rotate90);
     SetTurtle(&dummy_turtle);
     (*dr->tdd->Circle3D)(&dummy_turtle, dr, &viewparam,
     primitive.data[0]*2);
     PopMatrix();
     dummy_turtle.position[1] = 0;  // back default 0 
     }
     else
     Ferr <<"Caribu<prim2triangle> Warning: not enough parameters for cylinder!\n" ;
   */
   break;
   
 case CONE:
   Ferr <<"Caribu<prim2triangle> cone not yet implemeted\n" ;
   /* if(primitive.num_items >=3) {
      SetTurtle(&dummy_turtle);
      
      dummy_turtle.line_width = 2*primitive.data[0];
      
      (*dr->tdd->StartNode)(&dummy_turtle, dr, &viewparam,
      primitive.data[2],'\0');
      
      PushMatrix();
      Rotate(rotate90);
      SetTurtle(&dummy_turtle);
      (*dr->tdd->Circle3D)(&dummy_turtle, dr, &viewparam,
      primitive.data[0]*2);
      PopMatrix();
      
      PushMatrix();	  
      translate[1] = primitive.data[2];
      Translate(translate);
      SetTurtle(&dummy_turtle);
      dummy_turtle.line_width = 2*primitive.data[1];
      
      (*dr->tdd->EndNode)(&dummy_turtle, dr, &viewparam,'\0');
      
      Rotate(rotate90);
      SetTurtle(&dummy_turtle);
      
      (*dr->tdd->Circle3D)(&dummy_turtle, dr, &viewparam,
      primitive.data[1]*2);
      PopMatrix();
      
      dummy_turtle.position[1] = 0;  // back default 0 
      }
      else
      Ferr <<"Caribu<prim2triangle> Warning: not enough parameters for cylinder!\n" ;
   */
   break;
 case PRISM:
   Ferr <<"Caribu<prim2triangle> prism not yet implemeted\n" ;
   /* if(primitive.num_items >= 3) {
      for(x=0;x<2;x++) 
      for(y=0;y<2;y++) 
      for(z=0;z<2;z++) {
      P[x][y][z][0] = x * primitive.data[0];
      P[x][y][z][1] = y * primitive.data[1];
      P[x][y][z][2] = z * primitive.data[2];
      }
      
      DrawRectangle(P[0][0][0], P[0][0][1],  P[0][1][1], P[0][1][0]);
      DrawRectangle(P[1][0][0], P[1][1][0],  P[1][1][1], P[1][0][1]);
      
      DrawRectangle(P[0][0][0], P[0][1][0],  P[1][1][0], P[1][0][0]);
      DrawRectangle(P[0][0][1], P[1][0][1],  P[1][1][1], P[0][1][1]);
      
      DrawRectangle(P[0][0][0], P[0][0][1],  P[1][0][1], P[1][0][0]);
      DrawRectangle(P[0][1][0], P[1][1][0],  P[1][1][1], P[0][1][1]);
      }
      else
      Ferr <<"Caribu<prim2triangle> Warning: not enough parameterss for prism!\n" ;
   */
   break;  
 default:
   Ferr <<"Caribu<prim2triangle> Warning: unknown primitive.\n" ;
}
/* Normaliastion by area of the isobarycentre of the prim */
//if(1)fprintf(stderr,"prim2triangle: Gi=%f =>",  G[1]);
if(area>0 ) {
  for(i=0;i<3;i++) G[i]/=area;   
}
//if(1)fprintf(stderr,": area=%f,  Gz= %f,  nbT=%d\n",
//					 area,  G[1], nbT);
return answer;
}


/*************************************************************************/
int OGLPrim::parse_line(char *line)
{
//line <= CSGetString(), T[nbT] return the triangles describing the analysed primitive
struct syntax_item *ptr;
register int i=0,ir=0;
char real[500],c,*str=NULL;
short record=0,lstr;
int answer=0;


//fprintf(stderr,"parse_geom_line: emty=%d, Gz=%f,  str=%s...\n",empty, G[1], line);
if(empty) 
  raz();

if(!isalpha(line[0]))
{
  str=line;
}//if alpha(line[0])
 else{
   if(analysing){//all the parameters are stored => dealing the corresponding primitive
     answer=prim2triangle();  
     empty=1;
     //fprintf(stderr,"parse_geom_line: apre prim2T emty=%d", empty);
   }
   //Start of the analysis of the following primitive
   ptr = syn;
   while(ptr->flag != -1) {
     lstr=strlen(ptr->keyword);
     if(strncmp(line, ptr->keyword,lstr) == 0) {
       primitive.flag = ptr->flag;
       primitive.num_items = 0;
       str=line+lstr;
       //fprintf(stderr,"OGLPrim.parse_line: keyword=%s\n",ptr->keyword);
       analysing=1;
       break; /* out of the loop going through all keywords */
     }
     ptr++;
   }    
   if(ptr->flag == -1){
     //  fprintf(stderr, "Unknown keyword: %s.\n", line);
     analysing=0;
   }
 }
 if(analysing){
   // Read str of numbers and store their in primitive.data[]
   for(i=0;;i++) {
     c=str[i];
     if(isspace(c) || c==0){
       if(record){
	 real[ir]=0;
	 //fprintf(stderr,"parse_geom_str: data[num_item=%d]=",primitive.num_items); 
	 primitive.data[primitive.num_items]=atof(real);
	 //fprintf(stderr,"%f\n",primitive.data[primitive.num_items]);
	 primitive.num_items++;
	 record=ir=0;
       }//if recording
     }// if space
     else{
       record=1;
       real[ir++]=c;
     }//else space
     if(c==0) break;
   }// for str[i]    
 }// if analysing
 //fprintf(stderr,"parse_geom_line: END  empty=%d\n",empty);
 return answer;
}

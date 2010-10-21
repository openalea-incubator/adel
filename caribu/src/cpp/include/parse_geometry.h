/* Class OGLPrim dealing the Open-GL primitives flow, coming from Radek's CommLib

 MC97*/

extern "C" {
#include "matrix.h"
}

#define NB_TRIANGLE_PRIM 100
#define NUM_FLOATS 50

struct primitive_type {
  int  flag;
  int  num_items;
  float data[NUM_FLOATS]; 
};
typedef struct primitive_type primitive_type;



class OGLPrim{
public:
  void raz(){area=G[0]=G[1]=G[2]=nbT=empty=0;}
  OGLPrim(){analysing=0; raz();  InitializeMatrixStack();}
  //contains the triangles resulting from the triangulation of the primitive
  float T[NB_TRIANGLE_PRIM][3][3];
  int nbT; // Number of triangles describing the primitive
  float area; //approximate surface computed by the sum of triangles surface
  float G[3]; //Gravity center of the primitive coming from the G of each triangles
  signed char analysing;

  //Parse a line of an OpenGL primitive flow
  int parse_line(char *line);

private:
  primitive_type primitive;
  signed char empty;
 
  short prim2triangle();
  void DrawRectangle(float *p1, float *p2, float *p3, float *p4);
  void DrawTriangle(float *p1, float *p2, float *p3);
};//OGLPrim

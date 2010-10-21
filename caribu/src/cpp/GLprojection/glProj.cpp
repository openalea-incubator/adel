/**********************************************************************/
/*         glProj.cpp                                    
	   projette une liste de triangles avec le Zbuffer d'OpenGL;
	   la couleur des faces des triangles etant définie par un indice entierunique 
	   pour chaque face de  triangle
*/
/**********************************************************************/
//         MC09
/**********************************************************************/
//compile: 
//MacOSX
// g++  -framework GLUT -framework OpenGL -Wall -o glproj glProj.cpp  -L"/System/Library/Frameworks/OpenGL.framework/Libraries" -lGL -lGLU -lm -lobjc -lstdc++ 

// Run:  echo "1 600 600 10" | glproj

// TODO list 9/1/09 MC
//todo voir gl liste pour ne pas recharger la scene: DONE
// todo: tester chaargement de la scene avnt le gluLookat et de chger le gluLookAt ss recharger la scene : DONE grace au Display list
// to do encapsuler le taf OpenGL dans une classe: DONE
// Todo: mettre les couleur de la lut et compter les occurences des idx de la LUT: done
//Todo: sortir proprement sans exit (exception?): A chiader...
//Todo: version sans gluMainLoop: DONE
// Todo: recuperer le depth buffer pour infinitise: 
// tester sur d'autres macheine (cartes graphiques ) et OS

using namespace std;
// #include <GL/glut.h>
//Sous MacOSX
//#include <GLUT/glut.h>

#if defined(__GNUC__) && (defined(__APPLE_CPP__) || defined(__APPLE_CC__))
 #include <GLUT/glut.h>
#else
 #include <GL/glut.h>
#endif

#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <exception>

#include <vector>

using namespace std;


// typedef void (*PFV)(void);

//      PFV set_terminate(PFV); //conflit compilo avec 

void t() {
  cerr << "terminate() called" << endl;
  exit(1);
}

class myexception: public exception{
  virtual const char* what() const throw()
  {
    return "My exception happened";
  }
} myex;


class Pixel{
public:
  GLfloat Tpix[4]; //R G B A
  //definir surcharge operator [] pour faire comme si un tableau
  Pixel(){}
  Pixel(float r, float g, float b){
    init(r,g,b,1.);
  }
  Pixel(float r, float g, float b, float a){
    init(r,g,b,a);
  }
  Pixel(float *P){
    init(P[0], P[1], P[2]);
  }
  void init(float r, float g, float b, float a=1){
    Tpix[0] = r;
    Tpix[0] = g;
    Tpix[2] = b;
    Tpix[3] = a;
  }

  GLfloat& operator [] (int i){
    return Tpix[i];
  }

  GLfloat * ptr(){return Tpix;}
  /*  Pixel & operator = (Pixel& p){
      for (int i=0; i<4; i++)
      Tpix[i]=p.Tpix[i];
      return *this;
      }
  */

};//class Pixel

  //class LUT
class LUT{
  //class enabling to convert an integer index in a unique RGB color pixel, 
  // and the inverse, maximizing the color contrast
  //M. Chelle - INRA -  january 2009
  //chelle@grignon.inra.fr
private:
  //GLfloat **  Tlut;
  //  GLfloat Tlut[200001][4];
    vector <Pixel> Tlut;
 
 int dpix,Npix;
  float my_epsi;
public:
  int N; 
  void init(int N);
  void idx2col(int i, Pixel & col);
  int col2idx(Pixel &col);
  void status();
  ~LUT(){Tlut.clear();}
};

int test_LUT(int Nmax){
  //test que la classe LUT fonctionne bien ie bijection entre indice entier et palette de couleur
  //MC09

  //init de la LUT
  LUT lut;
  int N, no, bug=0;
  Pixel pix;  
  //  for (N=4097; N<4098; N++){
  //for (N=9; N<10; N++){
  for (N=2; N<Nmax; N++){
    bug=0;
    if( N%1000 == 0) printf(">>>>   N=%d\n", N);
    lut.init(N);
    for (int i=0; i<N;i++){
      lut.idx2col(i,pix);
      no=lut.col2idx(pix);
      // printf("i=%d, no=%d, pix=(%.2f, %.2f, %.2f )\n",i,no,pix[0], pix[1], pix[2]);

      if(no!=i) {bug++; 
	printf("i=%d, no=%d, pix=(%.2f, %.2f, %.2f )\n",i,no,pix[0], pix[1], pix[2]);
	lut.status(); 
	return -2;}
    }
    //    lut.status();
    if(bug>0) {
      printf("<!> N=%d bug=%d !!: ", N, bug);
      lut.status();
    }
    // printf("--------\n\n");
  }
  return 0;
}


class Scene{
  // classe beta a remplacer par la classe Canopy dans la radiosité mixte...
  // construit juste une liste de triangles 
public:
  float Ttri[10][3][3];
  int N;
  Scene(){
    float  triangle[3][3];
    double  dx, dy, dz;
    double a=.1;
    N=10;
    triangle[0][0]=-a;  triangle[0][1]=-a;  triangle[0][2]=-a; 
    triangle[1][0]=a;  triangle[1][1]=-a;  triangle[1][2]=a;
    triangle[2][0]=-a;  triangle[2][1]=a;  triangle[2][2]=0;
  
    for(int i=0; i<N; i++){
      dx=1-drand48()*2;
      dx=1-drand48()*2;
      dx=1-drand48()*2;
    
      Ttri[i][0][0]=triangle[0][0]+dx*2; 
      Ttri[i][0][1]=triangle[0][1]+dy*2; 
      Ttri[i][0][2]=triangle[0][2]+dz*2;
 
      Ttri[i][1][0]=triangle[1][0]+dx*2; 
      Ttri[i][1][1]=triangle[1][1]+dy*2; 
      Ttri[i][1][2]=triangle[1][2]+dz*2;
 
      Ttri[i][2][0]=triangle[2][0]+dx*2; 
      Ttri[i][2][1]=triangle[2][1]+dy*2; 
      Ttri[i][2][2]=triangle[2][2]+dz*2; 
    }

    Ttri[0][0][0]=-1;
    Ttri[0][0][1]=-1;
    Ttri[0][0][2]=0;
 
    Ttri[0][1][0]=1; 
    Ttri[0][1][1]=-1; 
    Ttri[0][1][2]=0;
 
    Ttri[0][2][0]=0;
    Ttri[0][2][1]=1;
    Ttri[0][2][2]=0;

    Ttri[1][0][0]=0;
    Ttri[1][0][1]=-1;
    Ttri[1][0][2]=-1;
 
    Ttri[1][1][0]=0; 
    Ttri[1][1][1]=1; 
    Ttri[1][1][2]=-1;
 
    Ttri[1][2][0]=0;
    Ttri[1][2][1]=0;
    Ttri[1][2][2]=1;

    Ttri[2][0][0]=-1;
    Ttri[2][0][1]=0;
    Ttri[2][0][2]=-1;
 
    Ttri[2][1][0]=1; 
    Ttri[2][1][1]=0; 
    Ttri[2][1][2]=-1;
 
    Ttri[2][2][0]=0;
    Ttri[2][2][1]=0;
    Ttri[2][2][2]=1;
    
  }//init()

}; //class Scene

  float BigTab[5000][5000][3];

class GLscene{
  // Class permettant le calcul du Zbuffer par OpenGL
  // à partir d'une liste de couleur et utilisant la classe LUT pour realiser la bijection etntre indice du triangle et couleur du Zbuffer
  //MC09
private: 
  Pixel Black;
  GLuint id; //display list
  LUT lut;
  float ***pixels;
  
  void load(Scene &scene);
  //Analyse le contenu du COLOR bugger d'OpenGL
  void ScreenGrabLUT(); 

public:
  bool debug;
  int N, width, height;
  void init(Scene  &scene, int resX, int resY);
  void project(double t, double p);
  ~GLscene(){ glDeleteLists(id, 1) ;}
};//GLscene


// fonctions utiles a glutMainLoop
 //Truc pour passer une fction-membre à un gestionnaire de signal non C++ eg GLUT
// La fonction enveloppe utilise une variable comme objet:
double Theta=0, Phi=0;
int nrot=0, nrot_max;
GLscene* object_which_will_handle_signal;
void GLscene_project_wrapper(){
  object_which_will_handle_signal->project(Theta,Phi);
}
 
void idle(){
  if(nrot<nrot_max){
    Theta += M_PI/180.*5 ;
    Phi += M_PI/180.*5 ;
    printf("> idle: nrot=%d theta=%.1f\n",nrot++, Theta*180/M_PI);
    glutPostRedisplay() ;
  }
}

void clavier( unsigned char touche, int x, int y ){
  switch ( touche )  {
  case 'q' :
  case 27  :
    exit( 0 );
  }
  glutPostRedisplay();
}
// fin des fonction utile a glutMainLoop


//Fonctions-membres


  void GLscene::init(Scene &scene, int resX, int resY){
    //MC09
    Black.init(0,0,0,1.);
    id=1; //display list
    N=scene.N;

    //init de la LUT
    // test_LUT(200000); sucessfull !!! : MCjanv 2009
    lut.init(2*N);
    printf("GLscene.init(): LUT.N=%d\n",lut.N);
    int argc=0; 
    char **argv; 
    glutInit( &argc, argv);
    //glutInit( &argc, argv );
    /* The X Window System specific options parsed by glutInit are as follows:
       -display DISPLAY
       Specify the X server to connect to. If not specified, the value of the DISPLAY environment variable is used.
       -geometry W x H + X + Y
       Determines where window's should be created on the screen. The parameter following -geometry should be formatted as a standard X geometry specification. The effect of using this option is to change the GLUT initial size and initial position the same as if glutInitWindowSize or glutInitWindowPosition were called directly.
       -iconic
       Requests all top-level windows be created in an iconic state.
       -indirect
       Force the use of indirect OpenGL rendering contexts.
       -direct
       Force the use of direct OpenGL rendering contexts (not all GLX implementations support direct rendering contexts). A fatal error is generated if direct rendering is not supported by the OpenGL implementation.
       If neither -indirect or -direct are used to force a particular behavior, GLUT will attempt to use direct rendering if possible and otherwise fallback to indirect rendering.
    */
    glutInitDisplayMode( GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE );

    int screenX    =  glutGet(GLUT_SCREEN_WIDTH)*.9;
    int screenY    =  glutGet(GLUT_SCREEN_HEIGHT)*.9;
    printf(">> resolution ecran = X=%d, Y=%d\n", screenX, screenY);
    //  cout <<"resX [1-"<<screenX<<"] "; cin>>resX; 
    //cout<<"resY= [1"<<screenY<<"] "; cin>>resY;
    if(resX>screenX){
      resY=resY/(double) resX*screenX;
      resY=(resY>screenY)?screenY:resY;
      resX=screenX;
      printf("<!> taille fenetre > taille ecran => resolution recalc0ulée... \n"); 
    }else  
      if(resY>screenY){
	resX=resX/(double) resY*screenY;
	resY=screenY;
	printf("<!> taille fenetre > taille ecran => resolution recalc0ulée... \n"); 
      }
    printf(">> resolution recalculee /ecran = resX=%d, resY=%d\n", resX, resY);
   
    //maj variable classe GLscene
    width=resX;
    height=resY;
    pixels=(float ***)BigTab;
    glutInitWindowPosition(0, 0 );
    glutInitWindowSize( resX, resY);
    glutCreateWindow( "GLscene" );

    // Initialisation d'OpenGL
    glClearColor( 0, 0,.0, .0 );

    // Version avec display List : declarer une fois la scene avant les calculs de projection, +optimisé
    id = glGenLists(1);
    load(scene);
    int error = glGetError();
    if (error != GL_NO_ERROR) {
      std::cout << "An OpenGL error has occured: " << gluErrorString(error) << std::endl;
    }
    //  Version GLUT pour debug
    if(debug){
      //Mise en place des fonctions evenements 
      object_which_will_handle_signal=this;

      glutDisplayFunc( GLscene_project_wrapper );
      glutKeyboardFunc( clavier   );
      //glutReshapeFunc(  reshape   );
      glutIdleFunc(idle) ;
      
      //Entree dans la boucle principale
      glutMainLoop();
    }
 

  }



  void GLscene::project(double theta, double phi )
  {
    //Efface le buffer des couleurs et le Z buffer
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable( GL_LIGHTING );
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);

    glDisable(GL_CULL_FACE);
    glEnable( GL_DEPTH_TEST );
    //Type de rendu
    glShadeModel( GL_FLAT );   //Ombrage facettes planes
    // Force le mode de detremination des normales aux polygones
    // glFrontFace(GL_CW);
    glFrontFace(GL_CCW);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,1) ; // 0 pour front-side only
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT,  Black.ptr());

    // pas necessaire : glDisable(GL_BEND);
  
    //Type de projection et dimension clipping 
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho(-2,2,-2,2,-100,100);
   
    // Positionnt de la camera
    glMatrixMode( GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt( 2*sin(theta)*cos(phi), 2*sin(theta)*sin(phi),2*cos(theta), // oeil
	       0.0, 0.0, 0.0,  // pt vise'
	       0.0, 0, 1.0 );  // verticale de la camera 
    glCallList(id);
    glutSwapBuffers();
    int error = glGetError();
    if (error != GL_NO_ERROR) {
      std::cout << "An OpenGL error has occured: " << gluErrorString(error) << std::endl;
    }
    // Analyze des couleurs de l'images
    ScreenGrabLUT();
  }

void GLscene::ScreenGrabLUT() {

  // get some info about the screen size
  int sw      =  glutGet(GLUT_WINDOW_WIDTH);
  int sh      =  glutGet(GLUT_WINDOW_HEIGHT); // ne renvoit pas toujours la bonne taille !!
  // unsigned bpp     =  glutGet(GLUT_WINDOW_RGBA) ? 4 : 3;
  //	GLenum   format  =  (bpp==4) ? GL_RGBA : GL_RGB;

  printf("ScreenGrabLUT : sw=%d, sh=%d \n", sw, sh);
  if(sw != width || sh != height){
    printf(" <!> taille fenetre reelle != taille fenetre calcul=> arret immediat\n");
    exit(1);
  }

  // read from front buffer
  glReadBuffer(GL_FRONT);
   
  //  float pixels[sw][sh][3]; // bug !!! 
    
  //    glReadPixels(0,0,sw,sh,GL_RGB,GL_FLOAT,pixels);
  /* glReadPixels returns values from each	pixel with lower left
     corner at (x + i, y +	j) for 0<i<width  and 0<j<height.
     This pixel is	said to	be the ith pixel in the	jth row.
     Pixels are returned in row order from	the lowest to the
     highest row, left to right in	each row.
  */
    
  GLfloat *Tpix, *ppix;
  std::vector<GLfloat> Vpix;
  bool all_pix=true, avec_stl=true;
  if(avec_stl){
    Vpix.resize(sw* sh*3);
    glReadPixels(0,0,sw,sh,GL_RGB,GL_FLOAT,&Vpix[0]); 
  }else
    if(all_pix){
      ppix=Tpix=new GLfloat[sw*sh*3];
      glReadPixels(0,0,sw,sh,GL_RGB,GL_FLOAT,Tpix); 
    }
    
  int idx, fond=0, badpix=0;

  int *Tcol;
  Pixel pixel;
  Tcol=new int[lut.N];
  for (int ii=0; ii<lut.N; ii++)
    Tcol[ii]=0;
  printf("ScreenGrabLUT(): LUT.N=%d\n", lut.N);
  int id=0;
  for ( int x=0; x<sw; x++)
    for ( int y=0; y<sh;y++ ){
      //idx=lut.col2idx(pixels[x][y]);
      pixel[0]=Vpix[id++];
      pixel[1]=Vpix[id++];
      pixel[2]=Vpix[id++];
      idx=lut.col2idx(pixel);

      /*
      //avec stl: a tester
      std::vector<BYTE> m_pix;
      m_pix.resize(sw* sh);
      //glReadPixels( 0, 0, sw, sh, GL_RGB, GL_UNSIGNED_BYTE, reinterpret_cast<LPVOID>( &m_pix[0] ) );
      glReadPixels( 0, 0, m_nWidth, m_nHeight, GL_RGB, GL_UNSIGNED_BYTE, reinterpret_cast<GLvoid *>( &m_pix[0] ) );
      ou encore 
      vector<GLubyte> store;
      store.resize(Wide * Height * 4);
      glReadPixels(0, 0, Wide, Height, GL_RGBA, GL_UNSIGNED_BYTE, &store[0])
      */
	
	
      if(idx>=0)
	Tcol[idx]++;//Data
      else if(idx==-1)
	fond++;
      else
	//Erreur pixel pas correct
	badpix++;
    }
  printf("=> nb pix=%d : fond =%d, bad pixels=%d \n",
	 sw*sh, fond, badpix);
  for( int k=0; k<lut.N;k++)
    printf("    col(%d) =%d\n",k, Tcol[k]);

  if(avec_stl)
    Vpix.clear();
  else if(all_pix)
    delete [] Tpix;

  int midx=(int) (sw/2.), midy=(int) (sh/2.);
	
  //    printf(">>> Pixmap : width=%d, height=%d => pixel(%d, %d)=%.4f, %.4f,%.4f\n",  sw, sh, midx, midy, pixels[midx][midy][0], pixels[midx][midy][1], pixels[midx][midy][2] );

}




  void GLscene::load(Scene &scene){
    Pixel col;
    
    glNewList(id, GL_COMPILE);
    //glPushMatrix();
    
    for(int i=0; i<GLscene::N; i++){
      lut.idx2col(2*i,col);
      glMaterialfv( GL_BACK, GL_EMISSION, col.ptr() );
      lut.idx2col(2*i+1,col);
      glMaterialfv( GL_FRONT, GL_EMISSION, col.ptr() );
      glBegin (GL_POLYGON);

      glVertex3f (scene.Ttri[i][0][0],scene.Ttri[i][0][1],scene.Ttri[i][0][2]);
      glVertex3f (scene.Ttri[i][1][0],scene.Ttri[i][1][1],scene.Ttri[i][1][2]);
      glVertex3f (scene.Ttri[i][2][0],scene.Ttri[i][2][1],scene.Ttri[i][2][2]);
      glEnd();
    }
    //glPopMatrix();
    glEndList(); /* stop */
    return; 
  }



  void LUT::init(int Nval ){
    // Construit la LUT entre un indice compris entre  0 et N et un triplet de GLfloat RBV  
    // LUT=new GLfloat*[N];
    //  for (int i=0; i<N; i++) LUT[i]=new GLfloat[4];
    N=Nval;
    Tlut.resize(N);
    my_epsi=1/300.;
    Npix=ceil(pow(N,1/3.)); // Nb de pts au sein d'une composante couleur;
    dpix=255./ (Npix-1) ; //ecart de valeur entre 2 pts au sein d'une composante couleur

  
    if(pow(Npix,3.) < N) 
      printf("<!> pb....\n");

    if(1)  printf(">> LUT::init(): N=%d, Npix=%d, dpix=%d \n", N, Npix, dpix);
    int iR=255, iG=255, iB=255;
    int nR=0, nG=0, nB=0;
    for(int i=0;i<N;i++){
      if(nB>=Npix){
	iG-=dpix; iB=255;
	nG++; nB=0;
      } 
      if(nG>=Npix){
	iR-=dpix; iG=255; iB=255;
	nR++, nG=nB=0;
      }
      if(nR>=Npix)
	printf("<!> N=%d, i=%d: Npix=%d : Error depassement de la LUT\n", N, i, Npix);
      
      if(0) printf(">>>  LUT(%d)=(%d, %d, %d)\n", i,iR, iG, iB);
      
      Tlut[i][0] = (1 + iR/255. *254) /255.;
      Tlut[i][1] = (1 + iG/255. *254) /255.;
      Tlut[i][2] = (1 + iB/255. *254) /255.;
      
      iB-=dpix;
      nB++;
    }
  }
  
  void LUT::idx2col( int i, Pixel &col){
    // Convertir un indice compris entre  0 et N en un triplet de GLfloat RBV 
    col[0]=Tlut[i][0];
    col[1]=Tlut[i][1];
    col[2]=Tlut[i][2];
    col[3]=1;
    if(0) printf(">> idx2col()i=%d,  col=(%.2f, %.2f, %.2f )\n",i,col[0], col[1], col[2]);

  }

  int LUT::col2idx( Pixel& col){
    int idx, icol[3];
    double ipix[3];

    if( col[0] < my_epsi && col[1] < my_epsi && col[02] < my_epsi )
      return -1; //pixel du fond
  
    //cas general
    for(int i=0; i<3; i++){
      ipix[i]=(255*col[i]-1)*255./254;
      // le round est importabt a cause des erreurs num...
      icol[i]=  round((255-ipix[i])/(double) dpix); 
      if(0) printf(">>> i=%d: ipix=%.f, icol=%.2f, (int)=%d\n",i,ipix[i], (255-ipix[i])/(double) dpix, icol[i]);
    }
    idx = icol[0] * Npix*Npix   
      + icol[1] * Npix 
      + icol[2] ;
    if(0)  printf("> col2idx(): iR=%d, iG=%d, iB=%d => idx=%d => ", 
		  icol[0], icol[1], icol[2], idx);
    
    //Verif / a un bug de couleur / OpenGL
    //ie une couleur pas dans la LUT 
    Pixel pix;  
    idx2col(idx,pix);
    if( (fabs(pix[0]-col[0])<my_epsi) && (fabs(pix[1]-col[1])<my_epsi) && (fabs(pix[2]-col[2])<my_epsi) ) //Ne marche pas ??? => Ok pb d'epsi var globale quelque part :-( => chger en my_epsix
      //  if( pix[0]==col[0] && pix[1]==col[1] && pix[2]==col[2] )//temporaire ,mais marche !
      return idx;
    else{
      printf(">>>error:  LUT(%d)=(%.2f, %.2f, %.2f) - my_epsi=%.5f\n",idx, pix[0],pix[1],pix[2],my_epsi);
      printf("fabs(%.6f, %.6f, %.6f)\n", fabs(pix[0]-col[0]), fabs(pix[1]-col[1]), fabs(pix[2]-col[2]));
      return -2; //Error
      //a gerer plus tard avec les exceptions - MC09
    }
  }

  void LUT::status(){
    printf("> class LUT: N=%d, N**1/3=%.2f,  dpix=%d, Npix=%d\n", N, pow(N,1/3.), dpix, Npix);
  }


  int main(){
    printf("-> step 1\n");
    Scene soupe;
  
    int isglut, resX, resY; 
    printf("-> step 2\n");

    //version debug avec glutMainLoop 0: non !=0: oui
    cout<<"is glut? (0/1) :"; cin>>isglut; 
    //resX et resY : dimension de l'image pour le Zbuffer
    cout <<"resX= "; cin>>resX; 
    cout<<"resY= "; cin>>resY;
    // nb de chgt de point de vue de la camera eg 10
    cout<<"nrot_max= "; cin>>nrot_max;
    cout<<"-------\n"; 
  
    try{
      printf("-> step 3\n");

      GLscene glsoup; 
      double theta=0, phi=0;
 
      glsoup.debug=(isglut==0)?false:true;
      printf("-> step 4\n");

      glsoup.init(soupe, resX, resY);
  
      for (int ii=0;ii<nrot_max && !isglut;ii++){
	theta += M_PI/180.*5 ;
	phi += M_PI/180.*5 ;
	printf("> loop %d:  theta=%.1f\n",ii, theta*180/M_PI);
	glsoup.project(theta, phi);
      }
    }//fin try
    
    /* catch (exception& e)
       {
       cout << e.what() << endl;
       printf("on fait koi maintenant ...\n");
       }
    */
    catch (...)
      {
	printf("on fait koi maintenant ...\n");
      }
  
    return 0;
  }

  /* Notes - MC09

     1) Memes resultats (en terme de nb de pixels colore's)en onscreen (avec glutMainLoop) et en offscreen (isglut=0) 
     2) resolution ne peut pas exceder 670 x 670 sur le Mac !!!! => grosse limitation
 
  */

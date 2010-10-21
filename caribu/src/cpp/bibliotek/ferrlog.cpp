#include "ferrlog.h"
#include <cstdlib>
#include <iostream>

using namespace std ;

void ferrlog::open(char *filename)
{
#ifdef DEBUG_OBJ
  clog << "Objet: reçu demande d'ouverture de "<<filename<<endl ;
#endif

  // La destruction du fichier précédent n'est possible que si
  // aucun autre process ne l'utilise, l'ouverture est soumise
  // aux memes conditions et "resette" l'ancien ==> on le laisse.

  char *pcTmpName=NULL;
  pcTmpName=GetAllFileName(filename);
  if(pcTmpName) { 
    out = new ofstream(pcTmpName, ios::out) ;
    if (!out->good())
      {
	clog << "Pas pu ouvrir " << filename<< endl ;
	out = (ofstream*) NULL ;
      } 
    free (pcTmpName);
  }
    
#ifdef DEBUG_OBJ
  clog << "Le flux " <<  filename <<" est ouvert" << endl ;
#endif
  
} ;
//constructeur ferrlog

ferrlog::ferrlog(char *filename){
	open(filename);
      } 
// //  //   //    //   //  //  // //
ferrlog::~ferrlog()
{
#ifdef DEBUG_OBJ
  *out << "Destruction du flux Ferr" << endl ;
  clog << "Destruction du flux Ferr" << endl ;
#endif

    if (out->good()) {
      *out << "ferrlog stream close by ~ferrlog()" << endl ;
      *out << "\t(may be abnormal)" << endl ;
      out->flush() ;
      out->close() ;
    }
} ;
// destructeur ferrlog

// //  //   //    //   //  //  // //
ferrlog & ferrlog::flush (void) 
{
  clog.flush() ;
  if (out != NULL)
    out->flush() ;
  return *this ;
}
// flush

// //  //   //    //   //  //  // //
ferrlog & ferrlog::operator << (char *msg)
{
  clog << msg ;
if (out != NULL)
  *out << msg ;
  return *this ;
} ;
// operateur ferrlog << char *

// //  //   //    //   //  //  // //
ferrlog & ferrlog::operator << (const char *msg)
{
  if (msg != (char*) NULL) {
    if ( *msg == '\n' ) { // palliatif pour endl 
      clog << endl ;
      if (out != NULL)
	*out << endl ;
    } else {
      clog << msg ;
      if (out != NULL)
	*out << msg ;
    }
  }
  return *this ;
} ;
// operateur ferrlog << char *

// //  //   //    //   //  //  // //
ferrlog & ferrlog::operator << (char msg)
{
  if ( msg == '\n' ) { // palliatif pour endl 
    clog << endl ;
    if (out != NULL)
      *out << endl ;
  } else {
    clog << msg ;
    if (out != NULL)
      *out << msg ;
  }
  return *this ;
} ;
// operateur ferrlog << const char 

// //  //   //    //   //  //  // //
ferrlog & ferrlog::operator << (int msg)
{
  clog << msg ;
if (out != NULL)
  *out << msg ;
  return *this ;
} ;
// operateur ferrlog << int

// //  //   //    //   //  //  // //
ferrlog & ferrlog::operator << ( long int msg) 
{
  clog << msg ;
if (out != NULL)
  *out << msg ;
  return *this ;
} ;
// operateur ferrlog << long int

// //  //   //    //   //  //  // //
ferrlog & ferrlog::operator << ( unsigned int msg) 
{
  clog << msg ;
if (out != NULL)
  *out << msg ;
  return *this ;
} ;
// operateur ferrlog << unsigned int

// //  //   //    //   //  //  // //
ferrlog & ferrlog::operator << (float msg)
{
  clog << msg ;
if (out != NULL)
  *out << msg ;
  return *this ;
} ;
// operateur ferrlog << float

// //  //   //    //   //  //  // //
ferrlog & ferrlog::operator << (double msg)
{
  clog << msg ;
if (out != NULL)
  *out << msg ;
  return *this ;
} ;
// operateur ferrlog << double

// //  //   //    //   //  //  // //
ferrlog & ferrlog::operator << (void * msg)
{
  clog << msg ;
if (out != NULL)
  *out << msg ;
  return *this ;
} ;
// operateur ferrlog << void*

ferrlog & ferrlog::operator << ( ostream & other) 
{
  clog << other ;
  if (out != NULL) {
    *out << other ;
  }
  return *this ;
  
} ;

ferrlog & ferrlog::operator << ( string msg) 
{
  clog << msg ;
  if (out != NULL) {
    *out << msg ;
  }
  return *this ;
} ;

void ferrlog::close(void) {
    if (out->good()) {
      *out << "ferrlog stream close() called." << endl ;
      out->flush() ;
      out->close() ;
    }
}

// operateur ferrlog << string

// Historique !!
// operator ferrlog << ostream &
//void prog_terminate (int code) 
//{
//  Ferr << "Exit appelé avec valeur " << code << "\n" ;
//  delete (ferr) ;
//  exit (code) ;
//}

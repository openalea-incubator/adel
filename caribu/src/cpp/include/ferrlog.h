/*
 * Flux C++ de logs; les messages gérés par ce flux sont recopiés sur 
 * cerr ET dans le fichier dont le nom est fourni au constructeur
 * de l'objet ferrlog.
 * Bug ? 
 * aucune gestion des accès au fichier n'est réalisée
 *
 */
#ifndef __FERRLOG_H
#define __FERRLOG_H

#include <fstream>
#include <string>
using namespace std;
#include <system.h>

class ferrlog
{
 private:
  ofstream *out ;
 public:
  ferrlog( char *filename) ;
  void open(char *filename);
  ~ferrlog() ;
  ferrlog &operator << ( char * msg) ;
  ferrlog &operator << ( const char *msg ) ;
  ferrlog &operator << ( char msg) ;
  ferrlog &operator << ( int msg) ;
  ferrlog &operator << ( long int msg) ;
  ferrlog &operator << ( unsigned int msg) ;
  ferrlog &operator << ( float msg) ;
  ferrlog &operator << ( double msg) ;
  ferrlog &operator << ( void *msg) ;
  ferrlog &operator << ( ostream & other) ;
  ferrlog &operator << ( string msg) ;
  ferrlog &flush (void) ;
  void close(void);

} ;

//#define Ferr *ferr

extern ferrlog Ferr ;               //ferrlog.cpp

//void prog_terminate (int code) ;
/* assure la destruction du ferrlog avant exit;
 * puis apelle exit(code)
 */

#endif

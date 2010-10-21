#include <iostream>
using namespace std ;
#include <boost/python.hpp>
using namespace boost::python ;
#include <cstring>


int s2v(int argc, char **argv);
int s2v_wrap(boost::python::list args ) {
  int iSize = extract<int> (args.attr("__len__")()),
    iIb ;
  char ** argv = new char*[iSize+1] ;
  
#define NOM "s2v"

  argv[0] = new char[strlen(NOM)+1] ;
  strcpy (argv[0], NOM);

  for (iIb=0; iIb < iSize ; iIb ++) {
    const char * pcArg ;
    extract <const char*> earg(args[iIb]);
    if (earg.check()) {
      pcArg= earg() ;
      argv[iIb+1] = new char[strlen(pcArg)+1] ;
      strcpy (argv[iIb+1], pcArg);
    }
  } 

  int result = s2v (iSize+1, argv ) ;

  for (iIb = 0 ; iIb < iSize ;iIb ++) {
    delete [] argv[iIb] ;
  }
  delete [](argv) ;
  return result ;
}


int deux() {return 2;}

BOOST_PYTHON_MODULE(s2py) {
  def ("s2v", s2v_wrap ) ;
  def ("deux", deux ) ;
}

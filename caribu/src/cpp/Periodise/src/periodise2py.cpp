#include <iostream>
using namespace std ;
#include <boost/python.hpp>
using namespace boost::python ;
#include <cstring>


int periodise(int argc, char **argv);
int periodise_wrap(boost::python::list args ) {
  int iSize = extract<int> (args.attr("__len__")()),
    iIb ;
  char ** argv = new char*[iSize+1] ;
  
#define NOM "periodise"

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

  int result = periodise (iSize+1, argv ) ;

  for (iIb = 0 ; iIb < iSize ;iIb ++) {
    delete [] argv[iIb] ;
  }
  delete [](argv) ;
  return result ;
}


int deux() {return 2;}

BOOST_PYTHON_MODULE(per2py) {
  def ("periodise", periodise_wrap ) ;
  def ("deux", deux ) ;
}

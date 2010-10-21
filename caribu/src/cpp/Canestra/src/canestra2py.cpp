#include <iostream>
using namespace std ;
#include <boost/python.hpp>
using namespace boost::python ;
#include <cstring>

int canestra(int argc, char **argv);
int canestra_wrap(boost::python::list args ) {
  int iSize = extract<int> (args.attr("__len__")()),
    iIb ;
  char ** argv = new char*[iSize] ;
  
  for (iIb=0; iIb < iSize ; iIb ++) {
    const char * pcArg ;
    extract <const char*> earg(args[iIb]);
    if (earg.check()) {
      pcArg= earg() ;
      argv[iIb] = new char[strlen(pcArg)+1] ;
      strcpy (argv[iIb], pcArg);
      // Ferr << "Argument "<<iIb<< " : "<< argv[iIb] << endl ;
    }
  } 

  int result = canestra (iSize, argv ) ;

  for (iIb = 0 ; iIb < iSize ;iIb ++) {
    delete [] argv[iIb] ;
  }
  delete [](argv) ;
  return result ;
}

BOOST_PYTHON_MODULE(canestra2py) {
  def ("canestra", canestra_wrap ) ; 
} 

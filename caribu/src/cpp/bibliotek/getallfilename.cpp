#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>

#include <cstdlib>

using namespace std ;

#include <system.h>

char *GetAllFileName(char *nom) {
  char *pcTmpName=NULL;
  pcTmpName = (char*) malloc(MAX_PATH_LEN) ;
  
  if(true){
    if (pcTmpName == NULL){
      clog<<__FILE__<<" : "<<__LINE__<<" : Plus de mémoire"<<endl;
      exit (1);
    }
    // fin de chaîne
    pcTmpName[0]=0;
    
    // this #def is to be set for every Win32 compiler to run
#ifndef WIN32
    // TEMP est une macro dans stytem.h
    string sTmpName =TEMP ;
#else

    // On Win32, the %TEMP% variable _is_ variable.
    // I don't trust Windows functions when compiled by Bcc32, so 
    //I find out %TEMP% by myself  Ferr << __FILE__ " : "<< __LINE__ << '\n' ;
    // , assuming C:\ is allways read/write

    string  sTmpName;
    // on cherche TEMp dans l'environnement
    char ccCommande[MAX_PATH_LEN]= "echo %TEMP% > " ;
    
    strcpy(pcTmpName,"C:\\");
    strcat (pcTmpName,nom);
    strcat (ccCommande,pcTmpName);
    
    int iTest =0;
    iTest = system (ccCommande);
    if (iTest != 0) {
      cerr <<__FILE__<< ": Pas pu executer "<<ccCommande<<endl ;
      exit (20);
    };
    ifstream fin (pcTmpName, ios::in);
    
    // est-ce que FLUX >> string reserve la memoire ?
    fin >> sTmpName ;
    //clog <<"stmpname = "<< sTmpName ;
    // Attention : XP ne renvoie rien si une variable n'est pas définie
    // (style posix) HA 07 2004
    if (( sTmpName[0] == '%' ) || ( sTmpName.size() == 0 )) {
      istringstream isIn (".");
      isIn >>sTmpName ;
      clog << "No %TEMP% environmental found.\n";
      clog << "From now on the temp dir is going to be .\\."<<endl ;
    }
    fin.close();
    sprintf(ccCommande,"%s %s", RM, pcTmpName) ;
    system(ccCommande); 
    
    sTmpName.append("\\");
#endif

    sTmpName.append(nom);
    
    istringstream ssIn(sTmpName);
    ssIn >> pcTmpName ;
  }
  // sprintf(pcTmpName,"C:\\Windows\\Temp\\");
  return pcTmpName ;
}

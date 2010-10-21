//#include "utils.h"
#include "ferrlog.h"
#include "system.h"

extern void Stop(void) ; // dans utils.h pr chque version de crb

/****************************************************************************/
void RunCommand(char *cmd){
  int retour;
  
  //debug
  FILE *pfTmp ;
  pfTmp = fopen ("RunCmd", "at") ;
  fprintf (pfTmp,"%s\n",cmd);
  fclose (pfTmp) ;
  // debug
  
  Ferr<<"Running ==> "<<cmd<<'\n';
  retour=system(cmd);
  if(retour!=0){
    Ferr <<"Aborted: return code= "<<retour<<'\n';
    Crash(__FILE__, __LINE__ ) ;
  }/*if plante*/
  else
    Ferr<<" Done successfully"<<'\n';
}
// RunCommand


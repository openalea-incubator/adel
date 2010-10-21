#include <math.h>
#include <cantools.h>

#define RANG_E 1E11
#define RANG_P 1E6
#define RANG_O 1E5
#define RANG_F 1E3


////////// / / / / / / / / DECODE_DATA / / / / / / / / / ///////////
void Decode_Data (CanData * res, double donnee ) 
{
/* Si donnee <= 1001 (1000 par defaut dans les fichiers caribu&Co), 
  * on force les valeurs de especeId, planteId, iCodeOrg, OrganeId 
  * et iPolygnId comme si 1 100011205001 avait été lu :
  * espece 1 , plante 11, organe de type 2 et d'ID 5
  * HA 10/2004 */
  double tmp ;

  if(donnee > 1001.) { // 1001 pour pbms de test egalite en flottants
    res->iTranspar = donnee >= 0 ? 1 : -1 ;     /* transparence */
    donnee = fabs(donnee ) ;
    res->iEspeceId = (int) (donnee / RANG_E) ; /* espece */
    tmp = donnee - res->iEspeceId * RANG_E ; 
    res-> iPlanteId = (int) (tmp / RANG_P) ;   /* id plante */
    tmp -= res -> iPlanteId * RANG_P ;
    res -> iCodeOrg = (int)(tmp / RANG_O) ;    /* code organe */
    tmp -= res -> iCodeOrg *  RANG_O ;
    res -> iOrganeId = (int)(tmp / RANG_F) ;   /* id organe */
    tmp -= res -> iOrganeId *  RANG_F ;
    res -> iPolygnId = (int)tmp ;              /* id polygone */
  } else {
    res->iEspeceId=1;
    res->iPlanteId=11;
    res->iCodeOrg=2;
    res->iOrganeId=5;
    res->iPolygnId=0;
  }
#ifdef DEBUG_ON
  fprintf (stderr, "%012.0lf -> %d %d %d %d %d\n", 
	   donnee, res->iEspeceId,  res-> iPlanteId, 
	   res -> iCodeOrg,res -> iOrganeId, res -> iPolygnId );
  fflush (stderr) ;
  fflush (stdout) ;
#endif

}

int Label_2_Index (CanData data) {
  return 10* data.iEspeceId + data.iCodeOrg ;
}

int Label_2_Indexf (double donnee ) {
  CanData data ;
  Decode_Data( &data, donnee ) ;
  return 10* data.iEspeceId + data.iCodeOrg ;
}

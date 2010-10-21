/*------------------------------------------------------------------------------
 *                                                                              
 *        Alinea.Nema : Nema Model
 *                                                                              
 *        Copyright 2010 INRA - CIRAD
 *                                                                              
 *        File author(s): Elmer Ccopa Rivera <elmer.ccopa-rivera@avignon.inra.fr>
 *                        Christophe Pradal <christophe.pradal@cirad.fr>    
 *                        Christian Fournier <Christian.Fournier@supagro.inra.fr>
 *                                                                              
 *        Distributed under the CeCILL-C License.                               
 *        See accompanying file LICENSE.txt 
 *                                                                              
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr                    
 *       
 *        $Id: $
 *                                                                       
 *-----------------------------------------------------------------------------*/

#include "export_plant.h"

#include "nema/Plant.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>
#include <iostream>
#include <assert.h>

#include "nema_util.h"

using namespace boost::python;
using namespace boost;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Wrapping methods
////////////////////////////////////////////////////////////////////////////////

// Get/Set Params (Grain,Root,Entity)

dict GetGrainParams( Grain& grain )
{
    dict d;
    vector<string> names(grain.GetParamNames());
    vector<double> values( grain.GetParams() );

    assert( names.size() == values.size() );

    int i = 0;
    for( i= 0; i < names.size(); i++ )
        d[names[i]] = values[i];
    return d;
}

void SetGrainParams( Grain& grain, dict kwds )
{
    double p = extract<double>(kwds["p"]);
    double tt = extract<double>(kwds["tgrain"]);
    double g = extract<double>(kwds["ggrain"]);
    double sink = extract<double>(kwds["sink"]);
    double alpha = extract<double>(kwds["alpha"]);
    double beta = extract<double>(kwds["beta"]);
    double TTexp = extract<double>(kwds["TTexp"]);
    double TTinit = extract<double>(kwds["TTinit"]);
    grain.SetParamGrain(p,tt, g, sink, alpha, beta, TTexp,TTinit);
}

dict GetRootParams( Root& root )
{
    dict d;
    vector<string> names(root.GetParamNames());
    vector<double> values( root.GetParams() );

    assert( names.size() == values.size() );

    int i = 0;
    for( i= 0; i < names.size(); i++ )
        d[names[i]] = values[i];
    return d;
}

void SetRootParams( Root& root, dict kwds )
{
  double Umax;
  if (kwds.has_key("Umax"))
    Umax = extract<double>(kwds["Umax"]);
  else {
    cout << "Umax missing !! initilisation with 0 !!!" <<endl;
    Umax = 0;
  }
    double kroot1 = extract<double>(kwds["kroot1"]);
    double kroot2 = extract<double>(kwds["kroot2"]);
    double QonDmin = extract<double>(kwds["QonDmin"]);
    double beta1 = extract<double>(kwds["beta1"]);
    double beta2 = extract<double>(kwds["beta2"]);
    double NsoilMin = extract<double>(kwds["NsoilMin"]);
    double p = extract<double>(kwds["p"]);
    double sink = extract<double>(kwds["sink"]);
    double alpha = extract<double>(kwds["alpha"]);
    double beta = extract<double>(kwds["beta"]);
    double TTexp = extract<double>(kwds["TTexp"]);
    double degDMcoeff = extract<double>(kwds["degDMcoef"]);
    double degcoeff = extract<double>(kwds["degcoef"]);
    double TTinit = extract<double>(kwds["TTinit"]);
    root.SetParamRoot(Umax,kroot1,kroot2,QonDmin,beta1,beta2,NsoilMin,p,sink,alpha,beta,TTexp, degDMcoeff, degcoeff,TTinit);
}

dict GetEntityParams( Entity& entity )
{
    dict d;
    vector<string> names(entity.GetParamNames());
    vector<double> values( entity.GetParams() );

    assert( names.size() == values.size() );

    int i = 0;
    for( i= 0; i < names.size(); i++ )
        d[names[i]] = values[i];
    return d;
}

void SetEntityParams( Entity& entity, dict kwds )
{

    double p = extract<double>(kwds["p"]);
    double SynthRate = extract<double>(kwds["SynthRate"]);
    double DegRate = extract<double>(kwds["DegRate"]);
	double k1 = extract<double>(kwds["k1"]);
	double k2 = extract<double>(kwds["k2"]);
	double peff = extract<double>(kwds["peff"]);
	double pm1 = extract<double>(kwds["pm1"]);
	double pm2 = extract<double>(kwds["pm2"]);
	double pdeath = extract<double>(kwds["pdeath"]);
    double sink = extract<double>(kwds["sink"]);
    double alpha = extract<double>(kwds["alpha"]);
    double beta = extract<double>(kwds["beta"]);
    double TTexp = extract<double>(kwds["TTexp"]);
    double InsertionAngle = extract<double>(kwds["InsertionAngle"]);
	double DegDMRate = extract<double>(kwds["DegDMRate"]);
    entity.SetParamEntity(p,SynthRate,DegRate,k1,k2,peff,pm1,pm2,pdeath, sink, alpha, beta, TTexp,InsertionAngle,DegDMRate);
}

// Entity
DEFINE_GETSET(Entity,Nstructs0)
DEFINE_GETSET(Entity,Nphs0)
DEFINE_GETSET(Entity,Lengths0)
DEFINE_GETSET(Entity,Areas0)
DEFINE_GETSET(Entity,TTinits0)
DEFINE_GETSET(Entity,DMrems0)
DEFINE_GETSET(Entity,DMstructs0)

DEFINE_GETVECTOROFVECTOR(Entity,Nstructss)
DEFINE_GETVECTOROFVECTOR(Entity,Nphss)
DEFINE_GETVECTOROFVECTOR(Entity,TTinitss)
DEFINE_GETVECTOROFVECTOR(Entity,DMremss)
DEFINE_GETVECTOROFVECTOR(Entity,DMstructss)
DEFINE_GETVECTOROFVECTOR(Entity,AreaGreenss)
DEFINE_GETVECTOROFVECTOR(Entity,NphM2Maxss)
DEFINE_GETVECTOROFVECTOR(Entity,NphSynthRatess)
DEFINE_GETVECTOROFVECTOR(Entity,NphloemSynthRatess)
DEFINE_GETVECTOROFVECTOR(Entity,DMremSynthRatess)
DEFINE_GETVECTOROFVECTOR(Entity,Photoss)

//Plant

DEFINE_GETSET_PUBLIC(Plant,TTimes)
DEFINE_GETSET_PUBLIC(Plant,dTTs)
DEFINE_GETSET_PUBLIC(Plant,Nmobs)
DEFINE_GETSET_PUBLIC(Plant,DMGreenTots)
DEFINE_GETSET_PUBLIC(Plant,Imports)
DEFINE_GETSET_PUBLIC(Plant,ImportPots)
DEFINE_GETSET_PUBLIC(Plant,AreaGreenTots)
DEFINE_GETSET_PUBLIC(Plant,Productions)
DEFINE_GETSET_PUBLIC(Plant,remobilizedDMs)
DEFINE_GETSET_PUBLIC(Plant,Demands)
DEFINE_GETSET_PUBLIC(Plant,DMgrains)
DEFINE_GETSET_PUBLIC(Plant,Ngrains)
  



void InitialPlantExt(Plant& p, int nbsubsteps, double _Nmob, list _PARss,double _Nsoil)
{
  vector< vector <double> > par = listlist_vectorvector<double>( _PARss);
  p.InitialPlantExt(nbsubsteps,  _Nmob, par, _Nsoil);
}

void ActuPlantExt(Plant& p, int step, int nbsubsteps,double dTTstep,list _PARss)
{
  p.ActuPlantExt(step,nbsubsteps,dTTstep, listlist_vectorvector<double>( _PARss));
}

////////////////////////////////////////////////////////////////////////////////
// Wrapping methods
////////////////////////////////////////////////////////////////////////////////

void class_grain()
{
    class_< Grain > ("Grain", "Documentation of Grain class\n"
                              "                            ")
      //methods
        .def( init< string, string, string >("Module definition of Grain class") )
        .def("setParams", &Grain::SetParamGrain, "SetParams(p, tt, g, sink, alpha, beta, TTexp,TTinit)")
      //parameters
        .add_property( "p",&Grain::GetP, &Grain::SetP, 
                       "N fraction of the dry mass influx, p (dimensionless)" )
        .add_property( "tgrain",&Grain::GetTT, &Grain::SetTT,
                       "Duration of the endosperm cell division phase, tau (oCd)" )
        .add_property( "ggrain",&Grain::GetG, &Grain::SetG, "Relative rate of potential grain N filling during cell division, gamma (oCd-1)" )
        .add_property( "sink",&Grain::GetSink, &Grain::SetSink, "Relative sink strength of grains, delta_M_g, (oCd-1)" )
        .add_property( "alpha",&Grain::GetAlpha, &Grain::SetAlpha, "Parameter determining the shape of the Beta function for grains, alpha_g (dimensionless)" )
        .add_property( "beta",&Grain::GetBeta, &Grain::SetBeta, "Parameter determining the shape of the Beta function for grains, beta (dimensionless)" )
        .add_property( "TTexp",&Grain::GetTTExp, &Grain::SetTTExp, "Duration during which grains can accumulate dry mass, tt_g (oCd)" )
        .add_property( "TTinit",&Grain::GetTTinit, &Grain::SetTTinit, "Date of organ birth for sink function (oCd)" )
        .add_property( "params", GetGrainParams, SetGrainParams, 
                       "Get/Set the Parameters (p, tgrain, ggrain, sink, alpha, beta, TTexp)")
      //initial state
        .add_property( "N0",&Grain::GetNGrainInit, &Grain::SetNGrainInit, "Initial total N mass for grains, N_tot_g (g)" )
        .add_property( "DM0",&Grain::GetDMGrainInit, &Grain::SetDMGrainInit, "Initial total dry mass for grains, M_tot_g (g)" )

        ;
}
void class_root()
{
    class_< Root > ("Root", "Documentation of Root class\n")
      //methods
        .def( init< string, string >("Module definition of Root class") )
        .def("setParams", &Root::SetParamRoot)
      //parameters
        .add_property("Umax", &Root::GetUMAX, &Root::SetUMAX, "Theoretical maximum root N uptake at saturating soil N concentration, U_r,max, (g.m-3.oCd-1)")
        .add_property("kroot1", &Root::GetKroot1, &Root::SetKroot1, "Constant of the michaelis function reflecting HATS activity, beta_C, (g.m-3)")
        .add_property("kroot2", &Root::GetKroot2, &Root::SetKroot2, "Constant of the michaelis function reflecting HATS activity, beta_C, (g.m-3)")
        .add_property("QonDmin", &Root::GetQONDMIN, &Root::SetQONDMIN, "Rate constant of the linear function reflecting LATS activity, k_r,1, (g.m-3.oCd-1)")
        .add_property("beta1", &Root::GetBETA1, &Root::SetBETA1, "Coefficient for N availability effect on root N uptake, beta_N, (dimensionless)")
        .add_property("beta2", &Root::GetBETA2, &Root::SetBETA2, "Coefficient for C availability effect on root C uptake, beta_C (dimensionless)")
        .add_property("Nsoilmin", &Root::GetNsoilMin, &Root::SetNsoilMin, "Minimum nitrogen in the soil, Nsoilmin (g.m^3)")
        .add_property( "degDMcoef", &Root::GetDEGDMCOEF , &Root::SetDEGDMCOEF, "Relative degradation rates of remobilizable dry mass for roots, delta_M_r (oCd-1)" )
        .add_property( "degcoef", &Root::GetDEGCOEF , &Root::SetDEGCOEF, "Relative degradation rates of remobilizable N for roots, delta_N_r")
        .add_property( "p",&Root::GetP, &Root::SetP, "Proportion coefficient for N influx following dry mass influx into roots, p_r (dimensionless)" )
        .add_property( "sink",&Root::GetSink, &Root::SetSink, "Relative sink strength of roots, sigma_r (dimensionless)")
        .add_property( "alpha",&Root::GetAlpha, &Root::SetAlpha, "Parameter determining the shape of the Beta function for roots, alpha_r (dimensionless)")
        .add_property( "beta",&Root::GetBeta, &Root::SetBeta, "Parameter determining the shape of the Beta function for roots, beta_r (dimensionless)")
        .add_property( "TTexp",&Root::GetTTExp, &Root::SetTTExp, "Duration during which roots can accumulate dry mass, tt_r (oCd)" )
        .add_property( "TTinit",&Root::GetTTinit, &Root::SetTTinit, "Date of organ birth for sink function (oCd)" )
        .add_property( "params", GetRootParams, SetRootParams, 
                       "Get/Set the Parameters (Umax,kroot1,kroot2,QonDmin,beta1,beta2,Nsoilmin,p,sink,alpha,beta,TTexp, degDMcoef, degcoef,TTinit)")
       
      // initialisation
      .add_property( "DMrem0",&Root::GetDMremInit, &Root::SetDMremInit, "Initial remobilizable dry mass for roots, M_rem_r (g)" )
      .add_property( "DMstruct0",&Root::GetDMstructInit, &Root::SetDMstructInit, "Initial structural dry mass for roots, M_struct_r (g) " )
      .add_property( "Nrem0",&Root::GetNremInit, &Root::SetNremInit, "Initial root remobilizable N mass, N_rem_r (g)" )
      .add_property( "Nstruct0",&Root::GetNstructInit, &Root::SetNstructInit, "Initial structural N mass for roots, N_struct_r (g)" )
        ;
}
void class_environment()
{
    class_< Environment > ("Environment", "Documentation of Environment class\n")
        
        .def( init< string, string,string, string >("Module definition of Environment class") )
        .def("setParams", &Environment::SetParamEnv, "")
        .add_property("transmission", &Environment::GetTransmission, &Environment::SetTransmission, "Transmission, (dimensionless)")
        .add_property("coef_ext", &Environment::GetCoef_Ext, &Environment::SetCoef_Ext, "PAR extinction coefficient for vertical entities, k_vertical (m2.m-2)")
        ;
}
void class_entity()
{
    class_< Entity > ("Entity", "Documentation of Entity class\n")
        
        .def( init< string, string,int, string >("Module definition of Entity class") )
        .def("addEntity", &Entity::AddEntity, 
                "Add a new elemet to the entity object.\n"
                "addEntity(Nstruct,Nph,Length,Area,TTinit,DMrem,DMstruct)")
        .def("setParams", &Entity::SetParamEntity)
        .add_property( "p", &Entity::GetP, &Entity::SetP, "Proportion coefficient for N influx following dry mass influx into entities tp, delta_tp,i (g.g-1.oC.d-1)" )

        .add_property( "sink",&Entity::GetSink, &Entity::SetSink, "Relative sink strength of entities tp, delta_tp,i (dimensionless)" )
        .add_property( "alpha",&Entity::GetAlpha, &Entity::SetAlpha, "Parameter determining the sahpe of the Beta function for entities tp, alpha_tp (dimensionless)" )
        .add_property( "beta",&Entity::GetBeta, &Entity::SetBeta, "Parameter determining the sahpe of the Beta function for entities tp, beta_tp (dimensionless)" )
        .add_property( "TTexp",&Entity::GetTTExp, &Entity::SetTTExp, "Duration during which entities tp can accumulate dry mass, tt_tp (oCd)" )


        .add_property( "SYNTCOEF", &Entity::GetSYNTCOEF, &Entity::SetSYNTCOEF, "Relative rate of photosynthetic N synthesis associated to xylem influx for entities tp, delta_tp,i (g.g-1.oCd-1)" )
        .add_property( "DEGCOEF", &Entity::GetDEGCOEF, &Entity::SetDEGCOEF, "Relative degradation rates of remobilizable dry mass for entities tp, delta_tp (oCd-1)")
        .add_property( "k1", &Entity::GetK1, &Entity::SetK1, "Michaelis-Menten constants defining photosynthetic N synthesis associated to xilem influx for entities tp, k_tp,1 (g.g-1)" )
        .add_property( "k2", &Entity::GetK2, &Entity::SetK2, "Michaelis-Menten constants defining photosynthetic N synthesis associated to xilem influx for entities tp, k_tp,2 (J.m-2.d-1)" )
        .add_property( "peff", &Entity::GetPeff, &Entity::SetPeff, "Photosynthetic efficiency, epsilon_tp (g.J-1)")
        .add_property( "pm1", &Entity::GetPm1, &Entity::SetPm1, "Intercept of the linear relationship linking photosynthetic area, omega_tp,1 (g.m-2.d-1)")
        .add_property( "pm2", &Entity::GetPm2, &Entity::SetPm2, "Slope of the linear relationship linking photosynthetic area, omega_tp,2 (d-1)" )
        .add_property( "pdeath", &Entity::GetPdeath, &Entity::SetPdeath, "Proportion of maximum specific N mass at which tissues die for entities tp, d_tp (dimensionless)")


        .add_property( "params", GetEntityParams, SetEntityParams, 
                       "Get/Set the Parameters with dictionary")
        .def( "__len__", &Entity::GetNbRow )
        ADD_VECTOR_PROPERTY(Entity,Nstructs0,"Initial structural N mass for entity tp,i; N_struct_tp,i (g)" )
        ADD_VECTOR_PROPERTY(Entity,Nphs0,"Initial photosynthetic N mass for entity tp,i; N_ph_tp,i (g)" )
        ADD_VECTOR_PROPERTY(Entity,Lengths0,"Initial total length for entity tp,i; L_tot_tp,i (m) " )
        ADD_VECTOR_PROPERTY(Entity,Areas0,"Initial total area for entity tp,i; A_tot_tp,i (m^2) " )
        ADD_VECTOR_PROPERTY(Entity,TTinits0,"Initial thermal time, t (oCd)" )
        ADD_VECTOR_PROPERTY(Entity,DMrems0,"Initial remobilizable dry mass for entity tp,i; M_rem_tp,i (g)" )
        ADD_VECTOR_PROPERTY(Entity,DMstructs0,"Initial structuralremobilizable dry mass for entity tp,i; M_struct_tp,i (g)" )
     
      ADD_VECTOR_OUTPUT(Entity,Nstructss,"Structural N mass for entity tp,i; N_struct_tp,i (g)")
      ADD_VECTOR_OUTPUT(Entity,Nphss,"Photosynthetic N mass for entity tp,i; N_ph_tp,i (g)")
      ADD_VECTOR_OUTPUT(Entity,TTinitss,"Thermal time, t (oCd)")
      ADD_VECTOR_OUTPUT(Entity,DMremss,"Remobilizable dry mass for entity tp,i; M_rem_tp,i (g)")
      ADD_VECTOR_OUTPUT(Entity,DMstructss,"Structuralremobilizable dry mass for entity tp,i; M_struct_tp,i (g)")
      ADD_VECTOR_OUTPUT(Entity,AreaGreenss,"Exposed area for entity tp,i; A_exp_tp,i (m^2) ")
      ADD_VECTOR_OUTPUT(Entity,NphM2Maxss,"Photosynthetic N mass per unit area, P_max_tp (g.d-1) ")
      ADD_VECTOR_OUTPUT(Entity,NphSynthRatess,"Photosynthetic N synthesis rate from xylem N for entity tp,i; S_Nph,xylem_tp,i (g.d-1) ")
      ADD_VECTOR_OUTPUT(Entity,NphloemSynthRatess,"Photosynthetic N synthesis rate from phloem N for entity tp,i; S_Nph,phloem_tp,i (g.d-1)")
      ADD_VECTOR_OUTPUT(Entity,DMremSynthRatess,"Remobilizable N synthesis rate for roots, S_Nrem_r (g.d-1)")
      ADD_VECTOR_OUTPUT(Entity,Photoss,"Dry mass production by entity tp,i; P_tp,i(g.d-1)")

        ;


}

void class_plant()
{
    class_< Plant > ("Plant")
        .def( init< const Grain& , 
                    const Root& , 
                    const Environment& , 
                    const Entity& , 
                    const Entity& , 
                    const Entity& , 
                    const Entity& , 
                    const Entity&  >() )
        .def("InitialPlant", &Plant::InitialPlant )
        .def("InitialPlantExt", InitialPlantExt )
        .def("ResetPlant", &Plant::ResetPlant )
        .def("ActuPlant", &Plant::ActuPlant)
        .def("ActuPlantExt", ActuPlantExt)
        .def("Display", &Plant::Affich)
        .def( "GetGrain", &Plant::GetGrain, return_internal_reference<1>() )
        .def( "GetRoot", &Plant::GetRoot, return_internal_reference<1>() )
        .def( "GetEnvironment", &Plant::GetEnvironment, return_internal_reference<1>() )
        .def( "GetLamina", &Plant::GetLamina, return_internal_reference<1>() )
        .def( "GetSheath", &Plant::GetSheath, return_internal_reference<1>() )
        .def( "GetInternode", &Plant::GetInternode, return_internal_reference<1>() )
        .def( "GetPeduncle", &Plant::GetPeduncle, return_internal_reference<1>() )
        .def( "GetChaff", &Plant::GetChaff, return_internal_reference<1>() )
        .def_readwrite( "Nsoil", &Plant::Nsoil )
    ADD_VECTOR_PROPERTY(Plant,TTimes,"Thermal time, t (oCd)")
ADD_VECTOR_PROPERTY(Plant,dTTs,"Thermal time within one day, t (oCd)")
ADD_VECTOR_PROPERTY(Plant,Nmobs,"Mobile N mass concentration, N_mob_c (g.g-1)")
ADD_VECTOR_PROPERTY(Plant,DMGreenTots,"")
ADD_VECTOR_PROPERTY(Plant,Imports,"Doc")
ADD_VECTOR_PROPERTY(Plant,ImportPots,"Potential root N uptake rate per unit root mass (g.g.-1.d-1)")
ADD_VECTOR_PROPERTY(Plant,AreaGreenTots,"Total area for entity tp,i; A_tot_tp,i (m^2)")
ADD_VECTOR_PROPERTY(Plant,Productions,"Doc")
ADD_VECTOR_PROPERTY(Plant,remobilizedDMs,"Remobilizable N degradation rate for roots, D_Nrem_r ()")
ADD_VECTOR_PROPERTY(Plant,Demands,"")
ADD_VECTOR_PROPERTY(Plant,DMgrains,"Total dry mass for grains, M_tot_g (g)")
ADD_VECTOR_PROPERTY(Plant,Ngrains,"Total N mass for grains, N_tot_g (g)")
     

      ;
};

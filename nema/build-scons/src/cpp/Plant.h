/**
 * @file Plant.h
 * @brief Declaration of class Plant
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifndef PLANT_H_
#define PLANT_H_

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamPlant.h"
#include "Grain.h"
#include "Root.h"
#include "Entity.h"
#include "Environment.h"

// System includes ////////////////////////////////////////// System include //
#include <vector>
#include<string>

class Plant:public ParamPlant
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	Plant();
	Plant(std::string fileNamePlant,std::string fileNameParamEnv,std::string fileNameNEnv, std::string fileNameTEnv, std::string fileNameMeteoEnv,
		std::string fileNameParamGrain, std::string fileNameVar0Grain,
		std::string fileNameParamLamina, std::string fileNameVar0Lamina, int nrowLamina,
		std::string fileNameParamSheath, std::string fileNameVar0Sheath, int nrowSheath,
		std::string fileNameParamInternode, std::string fileNameVar0Internode, int nrowInternode,
		std::string fileNameParamPeduncle, std::string fileNameVar0Peduncle, int nrowPeduncle,
		std::string fileNameParamChaff, std::string fileNameVar0Chaff, int nrowChaff,
		
		
		std::string fileNameRoot,
		//std::string fileNameRoot, std::string fileNameVar0Root, //ELMER
		
		std::string outputdir);

    Plant( const Grain& grain, 
           const Root& root, 
           const Environment& environment, 
           const Entity& lamina, 
           const Entity& sheath, 
           const Entity& internode, 
           const Entity& peduncle, 
           const Entity& chaff);

    virtual ~Plant();

// Public Methods /////////////////////////////////////////// Public Methods //
public:

	void InitialPlant(int nbsubsteps);
	void InitialPlantExt(int nbsubsteps,double _Nmob,std::vector<std::vector<double> > _PARss,double _Nsoil);
	void ResetPlant();
	//void ComputeFluxes(int step);
	void ActuPlant(int step, int nbsubsteps);
	void ActuPlantExt(int step, int nbsubsteps,double dTTstep,std::vector<std::vector<double> > _PARss);
	void Affich(std::string outputdir,int nbsubsteps);
//WQL-SA
	void Random(int randUG,ParamPlant* orig_Plant,ParamPlant* Sigmas_Plant,ParamEnv* orig_Env,ParamEnv* Sigmas_Env,ParamGrain* orig_Grain,ParamGrain* Sigmas_Grain,
		ParamRoot* orig_Root,ParamRoot* Sigmas_Root,ParamEntity* orig_lamina,ParamEntity* Sigmas_lamina,ParamEntity* orig_sheath,ParamEntity* Sigmas_sheath,
		ParamEntity* orig_internode,ParamEntity* Sigmas_internode,ParamEntity* orig_peduncle,ParamEntity* Sigmas_peduncle,
		ParamEntity* orig_chaff,ParamEntity* Sigmas_chaff);
	//Pick out all the parameters with None zero sigmas from each type of Param class
	//set vector NNP_params_size to record how many NNP params for each
	int SetNNP(ParamPlant* Sigmas_Plant,ParamEnv* Sigmas_Env,ParamGrain* Sigmas_Grain,ParamRoot* Sigmas_Root,
		ParamEntity* Sigmas_lamina,ParamEntity* Sigmas_sheath,ParamEntity* Sigmas_internode,ParamEntity* Sigmas_peduncle,ParamEntity* Sigmas_chaff);
	std::vector<int> NNP_params_size;
	std::vector<int> size_ladder;
	int NbNNP;
	//Always after use SetNNP function!//put all original NNP params value and Names together
	std::vector<double> SetAllParamsVector();	
	std::vector<double> AllParamVectorNNP;
	std::vector<std::string> AllParamNameVectorNNP;
	//Change the required parameter's value according to the index of NNP
	void Modify_Params(Plant* object,int j);
	//copy
	void Copy(Plant* source);
//WQL-SA

	// ADD CPL
	Grain* GetGrain( ) {return grain;};
	Root* GetRoot( ) {return root;};
	Environment* GetEnvironment( ) {return environment;};
	Entity* GetLamina( ) {return lamina;};
	Entity* GetSheath( ) {return sheath;};
	Entity* GetInternode( ) {return internode;};
	Entity* GetPeduncle( ) {return peduncle;};
	Entity* GetChaff( ) {return chaff;};

	
// Private Methods ///////////////////////////////////////// Private Methods //
private:
	int nrow(int type);
	//int nrow_exposed(int type);
	Entity* TypeToClass(int type);
	void ComputeDMGreenTot(std::vector<int> typs);
	//void ComputeDMGreenTot(int ntype);
	void ComputeNmob(std::vector<int> typs,double ConcNmob,std::vector<std::vector<double> > PARss, double Absorb, double NInfluxGrain, double NInfluxRoot);
	double ConcNmob(double Nmob,double DMgreen);
	void ComputeProduction(std::vector<int> typs,std::vector<std::vector<double> > PARss, int time, int nbsubsteps, double ttime);
	void ComputeRemobilizedDM(std::vector<int> typs, double dTT, double TT);
	void ComputeDemand(std::vector<int> typs, double ttime);
	void ComputeOrganHeightMaxFromSoil(std::vector<int> typs);
	double HeightPreviousOrgan(int type, int row);
	std::vector<double> HeightViss(int type);
	void ComputeVisibleAreas(std::vector<int> typs);
	void ActuRowExposed(std::vector<int> typs);
	double MeanHeightExposedPart(int type, int row);
	std::vector<double> GAI_cums(int type, int type2);
	void ComputeAreaGreenTot(std::vector<int> typs);
	std::vector<double> CumulGAI(int typ, std::vector<std::vector<double> > gaiss);

// Attributes /////////////////////////////////////////////////// Attributes //
private:
	/*ParamPlant* parameters;
	std::vector<double> paramgrains;
	std::vector<double> paramroots;
	std::vector<double> paramlaminas;
	std::vector<double> paramsheaths;
	std::vector<double> paraminternodes;
	std::vector<double> paramchaffs;
	std::vector<double> paramenvs;
	std::vector<std::string> paramNamegrains;
	std::vector<std::string> paramNameroots;
	std::vector<std::string> paramNamelaminas;
	std::vector<std::string> paramNamesheaths;
	std::vector<std::string> paramNameinternodes;
	std::vector<std::string> paramNamechaffs;
	std::vector<std::string> paramNameenvs;*/
//public:
	Grain* grain;
	Root* root;
	Environment* environment;
	Entity* lamina;
	Entity* sheath;
	Entity* internode;
	Entity* peduncle;
	Entity* chaff;
	//int ntype;
public:
	int Time;
	double TTime;
	std::vector<double> TTimes;
	std::vector<double> dTTs;
	double dTT;
	double plt_density;

	std::vector<int> types;
	double Nsoil;
	double Nmob;
	std::vector<double> Nmobs;
	double DMGreenTot;
	std::vector<double> DMGreenTots;
	std::vector<double> Imports;
	std::vector<double> ImportPots;
	std::vector<std::vector<std::vector<double> > > gai_cumsss;
	std::vector<std::vector<double> > PARss;
	std::vector<std::vector<std::vector<double> > > PARsss;
	std::vector<double> AreaGreenTots;
	//std::vector<double> Forcs;
	std::vector<std::vector<double> > HeightMaxss;
	std::vector<std::vector<double> > visibleAreass;
	std::vector<std::vector<int> > rowss_exposed;

	double Production; 
	std::vector<double> Productions;
	double remobilizedDM;
	std::vector<double> remobilizedDMs;
	double Demand;     // en fait il vaudrait mieux les avoir en vecteur, de façon ?pouvoir les gader
	std::vector<double> Demands;

	std::vector<double> DMgrains;
	std::vector<double> Ngrains;

};

#endif /* PLANT_H_ */

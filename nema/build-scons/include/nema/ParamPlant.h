/**
 * @file ParamPlant.h
 * @brief Declaration of class ParamPlant
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifndef PARAMPLANT_H_
#define PARAMPLANT_H_

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamGrain.h"
#include "ParamRoot.h"
#include "ParamEntity.h"
#include "ParamEnv.h"

// System includes ////////////////////////////////////////// System include //
#include <vector>
#include <string>

class ParamPlant:public ParamEnv,public ParamGrain,public ParamRoot
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	ParamPlant();
    ParamPlant(std::string fileNamePlant,std::string fileNameParamEnv,std::string fileNameNEnv, std::string fileNameTEnv, std::string fileNameMeteoEnv,
		std::string fileNameParamGrain, std::string fileNameVar0Grain,
		std::string fileNameParamLamina, std::string fileNameVar0Lamina, int nrowLamina,
		std::string fileNameParamSheath, std::string fileNameVar0Sheath, int nrowSheath,
		std::string fileNameParamInternode, std::string fileNameVar0Internode, int nrowInternode,
		std::string fileNameParamPeduncle, std::string fileNameVar0Peduncle, int nrowPeduncle,
		std::string fileNameParamChaff, std::string fileNameVar0Chaff, int nrowChaff,
		
		std::string fileNameRoot
		//std::string fileNameRoot, std::string fileNameVar0Root //ELMER
		
		);
    virtual ~ParamPlant();

// Public Methods /////////////////////////////////////////// Public Methods //
public:
	//Individual parameter values
	double GetConcNmobMin();
	/*ParamGrain* GetParamGrain();
	ParamRoot* GetParamRoot();
	ParamEntity* GetParamLamina();
	ParamEntity* GetParamSheath();
	ParamEntity* GetParamInternode();
	ParamEntity* GetParamPeduncle();
	ParamEntity* GetParamChaff();
	ParamEnv* GetParamEnv();*/

	//Set individual parameter values//WQL
	/*double SetConcNmobMin(double value);
	ParamGrain* SetParamGrain();
	ParamRoot* SetParamRoot();
	ParamEntity* SetParamLamina();
	ParamEntity* SetParamSheath();
	ParamEntity* SetParamInternode();
	ParamEntity* SetParamPeduncle();
	ParamEntity* SetParamChaff();
	ParamEnv* SetParamEnv();*/

	//Vector of parameter values and names
	/*std::vector<std::string> GetParamNames();
	std::vector<double> GetParamgrains();
	std::vector<std::string> GetParamNamegrains();
	std::vector<double> GetParamroots();
	std::vector<std::string> GetParamNameroots();
	std::vector<double> GetParamlaminas();
	std::vector<std::string> GetParamNamelaminas();
	std::vector<double> GetParamsheaths();
	std::vector<std::string> GetParamNamesheaths();
	std::vector<double> GetParaminternodes();
	std::vector<std::string> GetParamNameinternodes();
	std::vector<double> GetParamchaffs();
	std::vector<std::string> GetParamNamechaffs();
	std::vector<double> GetParamenvs();
	std::vector<std::string> GetParamNameenvs();*/
	
	//double GetInfluxStem(int tTime);
	//double GetAbsorb(int tTime);
//WQL-SA
	void SetConcNmobMin(double value);
	void SetParamPlant(double concNmobMin);
	//NNP
	void SetNNP(ParamPlant* Sigmas);
	//random
	void Random(int randUG,ParamPlant* orig,ParamPlant* Sigmas);
	//copy
	void Copy(ParamPlant* source);
//WQL-SA

// Private Methods ///////////////////////////////////////// Private Methods //
private:
	double GetFirstDouble(std::string line);
	//double InterpolLin(int t, std::vector<double> fileTps, std::vector<double> fileVal);


// Attributes /////////////////////////////////////////////////// Attributes //
private:
	//Individual parameter values
	double concNmobMin;
	//ParamGrain* paramGrain;
	//ParamRoot* paramRoot;
public:
	ParamEntity* paramLamina;
	ParamEntity* paramSheath;
	ParamEntity* paramInternode;
	ParamEntity* paramPeduncle;
	ParamEntity* paramChaff;
	
	//ParamEnv* paramEnv;
	//std::vector<double> TpsForcs;
	//std::vector<double> InfluxStems;
	//std::vector<double> Absorbs;

public:
	//Vector of parameter values
	std::vector<double> params;

	//Vector of parameter names
	std::vector<std::string> paramNames;

	//WQL-SA index vetor of None_zero_sigmas parameters
    std::vector<int> params_NNP;
	//std::vector<int> paramNames_NNP;
	//WQL-SA

};

#endif /* PARAMPLANT_H_ */

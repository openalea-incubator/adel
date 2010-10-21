/**
 * @file ParamPlant.cpp
 * @brief Implementation of class ParamPlant
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamPlant.h"
#include "mt_generator.h"

// System includes ////////////////////////////////////////// System include //
#include <fstream>
#include <sstream>
#include <iostream>

// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
ParamPlant::ParamPlant()
{
	paramLamina=new ParamEntity();
	paramSheath=new ParamEntity();
	paramInternode=new ParamEntity();
	paramPeduncle=new ParamEntity();
	paramChaff=new ParamEntity();
	//params.push_back(0);
	//paramNames.push_back("ConcNmobMin");
}
ParamPlant::ParamPlant(std::string fileNamePlant,std::string fileNameParamEnv,std::string fileNameNEnv, std::string fileNameTEnv, std::string fileNameMeteoEnv,
		std::string fileNameParamGrain, std::string fileNameVar0Grain,
		std::string fileNameParamLamina, std::string fileNameVar0Lamina, int nrowLamina,
		std::string fileNameParamSheath, std::string fileNameVar0Sheath, int nrowSheath,
		std::string fileNameParamInternode, std::string fileNameVar0Internode, int nrowInternode,
		std::string fileNameParamPeduncle, std::string fileNameVar0Peduncle, int nrowPeduncle,
		std::string fileNameParamChaff, std::string fileNameVar0Chaff, int nrowChaff,
		
		std::string fileNameRoot
		//std::string fileNameRoot, std::string fileNameVar0Root //ELMER
		
		):ParamEnv(fileNameParamEnv,fileNameNEnv, fileNameTEnv, fileNameMeteoEnv)
		,ParamGrain(fileNameParamGrain, fileNameVar0Grain)
		
		,ParamRoot(fileNameRoot)
		//,ParamRoot(fileNameRoot, fileNameVar0Root) //ELMER
		
{
	ifstream file;
	file.open(fileNamePlant.c_str());
	string line;
	getline(file, line);
	concNmobMin = GetFirstDouble(line);
	file.close();

	//vector of parameter values
	params.push_back(concNmobMin);
	//Vector of parameter names
	paramNames.push_back("ConcNmobMin");

	paramLamina = new ParamEntity(fileNameParamLamina, fileNameVar0Lamina, nrowLamina);
	paramSheath = new ParamEntity(fileNameParamSheath, fileNameVar0Sheath, nrowSheath);
	paramInternode = new ParamEntity(fileNameParamInternode, fileNameVar0Internode, nrowInternode);
	paramPeduncle = new ParamEntity(fileNameParamPeduncle, fileNameVar0Peduncle, nrowPeduncle);
	paramChaff = new ParamEntity(fileNameParamChaff, fileNameVar0Chaff, nrowChaff);

}

ParamPlant::~ParamPlant()
{
	if (paramLamina!=NULL) delete paramLamina;
	if (paramSheath!=NULL) delete paramSheath;
	if (paramInternode!=NULL) delete paramInternode;
	if (paramPeduncle!=NULL) delete paramPeduncle;
	if (paramChaff!=NULL) delete paramChaff;
}



////////////////////////////// Private Methods ////////////////////////////////

double ParamPlant::GetFirstDouble(string line)
{
	stringstream oss;
	oss << line;
	double tempValue;
	oss >> tempValue;
	return tempValue;
}

/////////////////////////////// Public Methods ////////////////////////////////

//Individual parameter values
double ParamPlant::GetConcNmobMin()
{
	return concNmobMin;
}
//WQL-SA
void ParamPlant::SetConcNmobMin(double value)
{
	concNmobMin=value;
}
void ParamPlant::SetParamPlant(double concNmobMin)
{
	SetConcNmobMin(concNmobMin);
}
//NNP
void ParamPlant::SetNNP(ParamPlant* Sigmas)
{
	int unsigned i;
	for(i=0;i<this->params.size();i++)
	{
		if (Sigmas->params[i]!=0)
		{
			this->params_NNP.push_back(i);
		}
	}
	//return params_NNP;
}
//random
void ParamPlant::Random(int randUG,ParamPlant* orig,ParamPlant* Sigmas)
{
	int unsigned i;
	switch(randUG)
  {
     case 1://uniform distribution
		 for(i=0;i<this->params.size();i++)
		 {
			 this->params[i]=orig->params[i]+((Sigmas->params[i]==0)?(0):(Sigmas->params[i]*2*(genrand_real3()-0.5)));
		 }
		 this->SetParamPlant(params[0]);
		 break;
	 case 0://Guauss distribution
		 for(i=0;i<this->params.size();i++)
		 {
			 this->params[i]=orig->params[i]+((Sigmas->params[i]==0)?(0):(Sigmas->params[i]*Moro_NormSInv(genrand_real3())));
		 }
		 this->SetParamPlant(params[0]);
		 break;
	 
  }
}

void ParamPlant::Copy(ParamPlant* source)
{

	//for father class copy
	this->ParamEnv::Copy((ParamEnv*)source);
	this->ParamGrain::Copy((ParamGrain*)source);
	this->ParamRoot::Copy((ParamRoot*)source);

	//Individual parameter values
	concNmobMin=source->concNmobMin;

	paramLamina->ParamEntity::Copy(source->paramLamina);
	paramSheath->ParamEntity::Copy(source->paramSheath);
	paramInternode->ParamEntity::Copy(source->paramInternode);
	paramPeduncle->ParamEntity::Copy(source->paramPeduncle);
	paramChaff->ParamEntity::Copy(source->paramChaff);


	//Vector of parameter values
	params=source->params;

	//Vector of parameter names
	paramNames=source->paramNames;

	//WQL-SA index vetor of None_zero_sigmas parameters
	params_NNP=source->params_NNP;
	//std::vector<int> paramNames_NNP;
	
}
//WQL-SA

// END OF File : ParamPlant.cpp

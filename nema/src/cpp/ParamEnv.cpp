/**
 * @file ParamEnv.cpp
 * @brief Implementation of class ParamEnv
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamEnv.h"
#include "mt_generator.h"

// System includes ////////////////////////////////////////// System include //
#include <fstream>
#include <sstream>
#include <iostream>

// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
ParamEnv::ParamEnv()
{
}
ParamEnv::ParamEnv(std::string fileNameParam,std::string fileNameN, std::string fileNameT, std::string fileNameMeteo)
{
	//Individual parameter values
	ifstream fileParam;
	fileParam.open(fileNameParam.c_str());
	string lineParam;
	getline(fileParam, lineParam);
	transmission = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	coef_ext = GetFirstDouble(lineParam); // il faut le virer par la suite car je ne l'utilise pas !!!
	fileParam.close();

	//Vector of parameter values
    params.push_back(transmission);
	params.push_back(coef_ext); // le mettre ailleurs par la suite car ce n'est plus un paramètre

	//Vector of parameter names
	paramNames.push_back("transmission");
	paramNames.push_back("coef_ext");

	//Environmental variables values through time
	ifstream fileN;
	fileN.open(fileNameN.c_str());
	string lineN;
	while(getline(fileN, lineN))
	{
		stringstream oss;
		oss << lineN;
		double tempValue;
		oss >> tempValue;
		Nsoils.push_back(tempValue);
	}
	fileN.close();
	ifstream fileT;
	fileT.open(fileNameT.c_str());
	string lineT;
	while(getline(fileT, lineT))
	{
		stringstream oss;
		oss << lineT;
		double tempValue;
		oss >> tempValue;
		TimeSoils.push_back(tempValue);
	}
	fileT.close();
	ifstream fileMeteo;
	fileMeteo.open(fileNameMeteo.c_str());
	double val;
	std::vector<double> meteos;
	while (fileMeteo >> val)
	{
		meteos.push_back(val);
	}
	int pos=0;
	for (unsigned int i=0;i<(meteos.size()/3);i++)
	{	
		Days.push_back(meteos[pos]);
		PAR0s.push_back(meteos[pos+1]);
		Tres.push_back(meteos[pos+2]);
		pos+=3;
	}
	
}

ParamEnv::~ParamEnv()
{
}


////////////////////////////// Private Methods ////////////////////////////////
double ParamEnv::GetFirstDouble(string line)
{
	stringstream oss;
	oss << line;
	double tempValue;
	oss >> tempValue;
	return tempValue;
}

//ADD ELMER
/*
	int ParamEnv::GetNbRow()
	{
	return nbrow;
	}	
*/
//ADD ELMER



/////////////////////////////// Public Methods ////////////////////////////////

//Individual parameter values
double ParamEnv::GetTransmission()
{
	return transmission;
}
double ParamEnv::GetCoef_Ext()
{
	return coef_ext;
}


// Vectors of parameter values and parameter names:
std::vector<double> ParamEnv::GetParams()
{
	 // voir si je dois mettre en attribut ce vecteur de paramètres
    return params;
}
std::vector<string> ParamEnv::GetParamNames()
{
	return paramNames;
}

//Environmental variables values through time
std::vector<double> ParamEnv::GetNsoils()
{
	return Nsoils;
}
std::vector<double> ParamEnv::GetTimeSoils()
{
	return TimeSoils;
}
std::vector<double> ParamEnv::GetDays() // pour la météo: changer de nom!!
{
	return Days;
}
std::vector<double> ParamEnv::GetPAR0s()
{
	return PAR0s;
}
std::vector<double> ParamEnv::GetTres()
{
	return Tres;
}

//WQL-SA
void ParamEnv::SetTransmission(double value)
{
    transmission=value;
}
void ParamEnv::SetCoef_Ext(double value)
{
	coef_ext=value;
}
void ParamEnv::SetParamEnv(double Transmission,double Coef_Ext)
{
	SetTransmission(Transmission);
	SetCoef_Ext(Coef_Ext);
}
//NNP
void ParamEnv::SetNNP(ParamEnv* Sigmas)
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
void ParamEnv::Random(int randUG,ParamEnv* orig,ParamEnv* Sigmas)
{
	int unsigned i;
	switch(randUG)
  {
     case 1://uniform distribution
		 for(i=0;i<this->params.size();i++)
		 {
			 this->params[i]=orig->params[i]+((Sigmas->params[i]==0)?(0):(Sigmas->params[i]*2*(genrand_real3()-0.5)));
		 }
		 this->SetParamEnv(params[0],params[1]);
		 break;
	 case 0://Guauss distribution
		 for(i=0;i<this->params.size();i++)
		 {
			 this->params[i]=orig->params[i]+((Sigmas->params[i]==0)?(0):(Sigmas->params[i]*Moro_NormSInv(genrand_real3())));
		 }
		 this->SetParamEnv(params[0],params[1]);
		 break;
	 
  }
}

void ParamEnv::Copy(ParamEnv* source)
{
	//Individual parameter values
	transmission=source->transmission;
	coef_ext=source->coef_ext;


	//Environmental variables values through time
	Nsoils.assign(source->Nsoils.begin(),source->Nsoils.end());
	TimeSoils.assign(source->TimeSoils.begin(),source->TimeSoils.end());
	Days=source->Days;
	PAR0s=source->PAR0s;
	Tres=source->Tres;
	/*Days.assign(source->Days.begin(),source->Days.end());
	PAR0s.assign(source->PAR0s.begin(),source->PAR0s.end());
	Tres.assign(source->Tres.begin(),source->Tres.end());*/

	//vector of parameter values
	params=source->params;

	//Vector of parameter names
	paramNames=source->paramNames;

	//WQL-SA index vetor of None_zero_sigmas parameters
	params_NNP=source->params_NNP;
	//std::vector<int> paramNames_NNP;

}

//WQL-SA

// END OF File : ParamEnv.cpp

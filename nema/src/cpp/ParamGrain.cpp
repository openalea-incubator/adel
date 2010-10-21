/**
 * @file ParamGrain.cpp
 * @brief Implementation of class ParamGrain
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamGrain.h"
#include "mt_generator.h"

// System includes ////////////////////////////////////////// System include //
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
ParamGrain::ParamGrain()
{
}
ParamGrain::ParamGrain(std::string fileNameParam, std::string fileNameVar0)
{
	//Individual parameter values
	ifstream fileParam;
	fileParam.open(fileNameParam.c_str());
	string lineParam;
	TTinit = 0;
	getline(fileParam, lineParam);
	p = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	tt = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	g = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	sink = GetFirstDouble(lineParam); 
	getline(fileParam, lineParam);
	alpha = GetFirstDouble(lineParam); 
	getline(fileParam, lineParam);
	beta = GetFirstDouble(lineParam); 
	getline(fileParam, lineParam);
	TTexp = GetFirstDouble(lineParam); 
	fileParam.close();

	//vector of parameter values
	params.push_back(p);
	params.push_back(tt);
	params.push_back(g);
	params.push_back(sink);
	params.push_back(alpha);
	params.push_back(beta);
	params.push_back(TTexp);
	params.push_back(TTinit);

	//Vector of parameter names
	paramNames.push_back("p");
	paramNames.push_back("tgrain");
	paramNames.push_back("ggrain");
	paramNames.push_back("sink");
	paramNames.push_back("alpha");
	paramNames.push_back("beta");
	paramNames.push_back("TTexp");
	paramNames.push_back("TTinit");

	// Variable Initial state
	ifstream fileVar0;
	fileVar0.open(fileNameVar0.c_str());
	string lineVar0;
	getline(fileVar0, lineVar0);
	NGrainInit = GetFirstDouble(lineVar0);
	//added Elmer 12/08/2010
	getline(fileVar0, lineVar0);
	DMGrainInit = GetFirstDouble(lineVar0);
	//added Elmer 12/08/2010
	fileVar0.close();
	
//	DMGrainInit = 0.12;//voir si pertinent: Cf les grains sont d�j?pr�sents ?floraison: les faire appara�tre bien plus t�t au stade double rides ?
	//cout << t << " " << a1 << " " << a2 << endl;
}

ParamGrain::~ParamGrain()
{
}



////////////////////////////// Private Methods ////////////////////////////////
double ParamGrain::GetFirstDouble(string line)
{
	stringstream oss;
	oss << line;
	double tempValue;
	oss >> tempValue;
	return tempValue;
}


/////////////////////////////// Public Methods ////////////////////////////////

//Individual parameter values
double ParamGrain::GetTTinit()
{
	return TTinit;
}

double ParamGrain::GetP()
{
	return p;
}
double ParamGrain::GetTT()
{
	return tt;
}
double ParamGrain::GetG()
{
	return g;
}
double ParamGrain::GetSink()
{
	return sink;
}
double ParamGrain::GetAlpha()
{
	return alpha;
}
double ParamGrain::GetBeta()
{
	return beta;
}
double ParamGrain::GetTTExp()
{
	return TTexp;
}

// Vectors of parameter values and parameter names:
std::vector<double> ParamGrain::GetParams()
{
	return params;
}
std::vector<string> ParamGrain::GetParamNames()
{
	return paramNames;
}

// Variable Initial state
double ParamGrain::GetNGrainInit()
{
	return NGrainInit;
}
void ParamGrain::SetNGrainInit(double value)
{
	NGrainInit=value;
}
double ParamGrain::GetDMGrainInit()
{
	return DMGrainInit;
}
void ParamGrain::SetDMGrainInit(double value)
{
	DMGrainInit=value;
}


//WQL-SA

void ParamGrain::SetP(double value)
{
	p=value;
        params[0] = p;
        assert(paramNames[ 0 ] == "p");
}
void ParamGrain::SetTT(double value)
{
	tt=value;
        params[1] = value;
        assert(paramNames[ 1 ] == "tgrain");
}
void ParamGrain::SetG(double value)
{
	g=value;
        params[2] = value;
        assert(paramNames[2] == "ggrain");
}
void ParamGrain::SetSink(double value)
{
	sink=value;
        params[3] = value;
        assert(paramNames[3] == "sink");
}
void ParamGrain::SetAlpha(double value)
{
	alpha=value;
        params[4] = value;
        assert(paramNames[4] == "alpha");
}
void ParamGrain::SetBeta(double value)
{
	beta=value;
        params[5] = value;
        assert(paramNames[5] == "beta");
}
void ParamGrain::SetTTExp(double value)
{
	TTexp=value;
        params[6] = value;
        assert(paramNames[6] == "TTexp");
}
void ParamGrain::SetTTinit(double value)
{
	p=value;
        params[7] = p;
        assert(paramNames[ 7 ] == "TTinit");
}

void ParamGrain::SetParamGrain(double p,double tt,double g,double sink,double alpha,double beta,double TTexp,double TTinit)
{
	SetP(p);
	SetTT(tt);
	SetG(g);
	SetSink(sink);
	SetAlpha(alpha);
	SetBeta(beta);
	SetTTExp(TTexp);
	SetTTinit(TTinit);
}
//NNP
void ParamGrain::SetNNP(ParamGrain* Sigmas)
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
void ParamGrain::Random(int randUG,ParamGrain* orig,ParamGrain* Sigmas)
{
	int unsigned i;
	switch(randUG)
  {
     case 1://uniform distribution
		for(i=0;i<this->params.size();i++)
		 {
			 this->params[i]=orig->params[i]+((Sigmas->params[i]==0)?(0):(Sigmas->params[i]*2*(genrand_real3()-0.5)));
		 }
		this->SetParamGrain(params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7]);
		 break;
	 case 0://Guauss distribution
		 for(i=0;i<this->params.size();i++)
		 {
			 this->params[i]=orig->params[i]+((Sigmas->params[i]==0)?(0):(Sigmas->params[i]*Moro_NormSInv(genrand_real3())));
		 }
		 this->SetParamGrain(params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7]);
		 break;
	 
  }
}

void ParamGrain::Copy(ParamGrain* source)
{
	//Individual parameter values
	NGrainInit=source->NGrainInit;
	DMGrainInit=source->DMGrainInit;
	p=source->p;
	tt=source->tt;
	g=source->g;
	sink=source->sink;
	alpha=source->alpha;
	beta=source->beta;
	TTexp=source->TTexp;
	TTinit=source->TTinit;


	//vector of parameter values
	params=source->params;
	//Vector of parameter names
	paramNames=source->paramNames;

	//WQL-SA index vetor of None_zero_sigmas parameters
	params_NNP=source->params_NNP;
	//std::vector<int> paramNames_NNP;
	
}
//WQL-SA
// END OF File : ParamGrain.cpp

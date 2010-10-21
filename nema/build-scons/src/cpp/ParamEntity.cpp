/**
 * @file ParamLamina.cpp
 * @brief Implementation of class ParamLamina
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamEntity.h"
#include "mt_generator.h"

// System includes ////////////////////////////////////////// System include //
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>

// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
ParamEntity::ParamEntity()
{
}
ParamEntity::ParamEntity(std::string fileNameParam, std::string fileNameVar0, int nrow)
{
        nbrow = nrow;
	//Individual parameter values
	ifstream fileParam;
	fileParam.open(fileNameParam.c_str());
	string lineParam;
	getline(fileParam, lineParam);
	p = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	syntcoef = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	degcoef = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	k1 = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	k2 = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	peff= GetFirstDouble(lineParam);
	/*getline(fileParam, lineParam);
	pe2= GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	pe3= GetFirstDouble(lineParam);*/
	getline(fileParam, lineParam);
	pm1 = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	pm2 = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	pdeath = GetFirstDouble(lineParam);
	getline(fileParam, lineParam);
	sink = GetFirstDouble(lineParam); // ex 1
	getline(fileParam, lineParam);
	alpha = GetFirstDouble(lineParam); // ex 2
	getline(fileParam, lineParam);
	beta = GetFirstDouble(lineParam); // ex 2
	getline(fileParam, lineParam);
	TTexp = GetFirstDouble(lineParam); // ex 600
	getline(fileParam, lineParam);
	InsertionAngle = GetFirstDouble(lineParam); 
	getline(fileParam, lineParam);
	degDMcoef = GetFirstDouble(lineParam); 
	fileParam.close();


	//vector of parameter values
	params.push_back(p);
	params.push_back(syntcoef);
	params.push_back(degcoef);
	params.push_back(k1);
	params.push_back(k2);
	params.push_back(peff);
	//params.push_back(pe2);
	//params.push_back(pe3);
	params.push_back(pm1);
	params.push_back(pm2);
	params.push_back(pdeath);
	params.push_back(sink);
	params.push_back(alpha);
	params.push_back(beta);
	params.push_back(TTexp);
	params.push_back(InsertionAngle);
	params.push_back(degDMcoef);

	//Vector of parameter names
	paramNames.push_back("p");
	paramNames.push_back("SynthRate");
	paramNames.push_back("DegRate");
	paramNames.push_back("k1");
	paramNames.push_back("k2");
	paramNames.push_back("peff");
	//paramNames.push_back("peff2");
	//paramNames.push_back("peff3");
	paramNames.push_back("pm1");
	paramNames.push_back("pm2");
	paramNames.push_back("pdeath");
	paramNames.push_back("sink");
	paramNames.push_back("alpha");	
	paramNames.push_back("beta");
	paramNames.push_back("TTexp");
	paramNames.push_back("InsertionAngle");
	paramNames.push_back("DegDMRate");

	//Variable Initial values
	ifstream fileVar0;
	fileVar0.open(fileNameVar0.c_str());
	string lineVar0;
	for (int r=0;r<nrow;r++)
	{
		getline(fileVar0, lineVar0);
		Nstructs0.push_back(GetFirstDouble(lineVar0));
		getline(fileVar0, lineVar0);
		Nphs0.push_back(GetFirstDouble(lineVar0));
		getline(fileVar0, lineVar0);
		Lengths0.push_back(GetFirstDouble(lineVar0));
		getline(fileVar0, lineVar0);
		Areas0.push_back(GetFirstDouble(lineVar0));
		/*getline(fileVar0, lineVar0);
		Agreens0.push_back(GetFirstDouble(lineVar0));*/
		getline(fileVar0, lineVar0); 
		TTinits0.push_back(GetFirstDouble(lineVar0));
		getline(fileVar0, lineVar0);
		DMrems0.push_back(GetFirstDouble(lineVar0));
		getline(fileVar0, lineVar0);
		DMstructs0.push_back(GetFirstDouble(lineVar0));
	}
	fileVar0.close();
}

ParamEntity::~ParamEntity()
{
}


////////////////////////////// Private Methods ////////////////////////////////
double ParamEntity::GetFirstDouble(std::string line)
{
	stringstream oss;
	oss << line;
	double tempValue;
	oss >> tempValue;
	return tempValue;
}

//ADD ELMER
int ParamEntity::GetNbRow()
	{
	return nbrow;
	}	
/*
std::vector<std::vector<double> > ParamEntity::GetNstructs() {
    return Nstructs;
}
*/
//ADD ELMER

//ADD CPL
void ParamEntity::AddEntity(double Nstruct, double Nph, double Length, 
                       double Area, double TTinit, double DMrem, double DMstruct)
{
    nbrow += 1;
    Nstructs0.push_back(Nstruct);
    Nphs0.push_back(Nph);
    Lengths0.push_back(Length);
    Areas0.push_back(Area);
    TTinits0.push_back(TTinit);
    DMrems0.push_back(DMrem);
    DMstructs0.push_back(DMstruct);
    
}
//ADD CPL


/////////////////////////////// Public Methods ////////////////////////////////

//Individual parameter values
double ParamEntity::GetP()
{
	return p;
}
double ParamEntity::GetSYNTCOEF()
{
	return syntcoef;
}
double ParamEntity::GetDEGCOEF()
{
	return degcoef;
}
double ParamEntity::GetK1()
{
	return k1;
}
double ParamEntity::GetK2()
{
	return k2;
}
double ParamEntity::GetPeff()
{
	return peff;
}
/*double ParamEntity::GetPe2()
{
	return pe2;
}
double ParamEntity::GetPe3()
{
	return pe3;
}*/
double ParamEntity::GetPm1()
{
	return pm1;
}
double ParamEntity::GetPm2()
{
	return pm2;
}
double ParamEntity::GetPdeath()
{
	return pdeath;
}
double ParamEntity::GetSink()
{
	return sink;
}
double ParamEntity::GetAlpha()
{
	return alpha;
}
double ParamEntity::GetBeta()
{
	return beta;
}
double ParamEntity::GetTTExp()
{
	return TTexp;
}
double ParamEntity::GetInsertionAngle()
{
	return InsertionAngle;
}
double ParamEntity::GetDEGDMCOEF()
{
	return degDMcoef;
}

// Vectors of parameter values and parameter names:
std::vector<double> ParamEntity::GetParams()
{
	 // voir si je dois mettre en attribut ce vecteur de paramètres
	return params;
}
std::vector<string> ParamEntity::GetParamNames()
{
	return paramNames;
}

//Variable Initial values
std::vector<double> ParamEntity::GetNstructs0()
{
	return Nstructs0;
}
std::vector<double> ParamEntity::GetNphs0()
{
	return Nphs0;
}
std::vector<double> ParamEntity::GetLengths0()
{
	return Lengths0;
}
std::vector<double> ParamEntity::GetAreas0()
{
	return Areas0;
}
/*std::vector<double> ParamEntity::GetAgreens0()
{
	return Agreens0;
}*/
std::vector<double> ParamEntity::GetTTinits0()
{
	return TTinits0;
}
std::vector<double> ParamEntity::GetDMrems0()
{
	return DMrems0;
}
std::vector<double> ParamEntity::GetDMstructs0()
{
	return DMstructs0;
}

void ParamEntity::SetNstructs0(const std::vector<double>& values)
{
    Nstructs0= values;
    nbrow = values.size();
    
}
void ParamEntity::SetNphs0(const std::vector<double>& values)
{
    Nphs0= values;
    nbrow = values.size();
}
void ParamEntity::SetLengths0(const std::vector<double>& values)
{
    Lengths0= values;
    nbrow = values.size();
}
void ParamEntity::SetAreas0(const std::vector<double>& values)
{
    Areas0= values;
    nbrow = values.size();
}
//std::vector<double> GetAgreens0()
void ParamEntity::SetTTinits0(const std::vector<double>& values)
{
    TTinits0= values;
    nbrow = values.size();
}
void ParamEntity::SetDMrems0(const std::vector<double>& values)
{
    DMrems0= values;
    nbrow = values.size();
}
void ParamEntity::SetDMstructs0(const std::vector<double>& values)
{
    DMstructs0= values;
    nbrow = values.size();
}
//WQL-SA

#define DEFINE_SET(_METH,_PAR,_INDEX)		\
void ParamEntity::_METH(double value)\
{\
	_PAR=value;\
        params[_INDEX] = _PAR;\
        assert(paramNames[ _INDEX ] == #_PAR);\
}

DEFINE_SET(SetP,p,0)
DEFINE_SET(SetSYNTCOEF,syntcoef,1)
DEFINE_SET(SetDEGCOEF,degcoef,2)
DEFINE_SET(SetK1,k1,3)
DEFINE_SET(SetK2,k2,4)
DEFINE_SET(SetPeff,peff,5)
DEFINE_SET(SetPm1,pm1,6)
DEFINE_SET(SetPm2,pm2,7)
DEFINE_SET(SetPdeath,pdeath,8)
DEFINE_SET(SetSink,sink,9)
DEFINE_SET(SetAlpha,alpha,10)
DEFINE_SET(SetBeta,beta,11)
DEFINE_SET(SetTTExp,TTexp,12)
DEFINE_SET(SetInsertionAngle,InsertionAngle,13)
DEFINE_SET(SetDEGDMCOEF,degDMcoef,14)

       


void ParamEntity::SetParamEntity(double p,double syntcoef,double degcoef,double k1,double k2,double peff,double pm1,double pm2,double pdeath,double sink,double alpha,double beta,double TTexp,double InsertionAngle, double degDMcoef)
{
	SetP(p);
	SetSYNTCOEF(syntcoef);
	SetDEGCOEF(degcoef);
	SetK1(k1);
	SetK2(k2);
	SetPeff(peff);
	//SetPe2(pe2);
	//SetPe3(pe3);
	SetPm1(pm1);
	SetPm2(pm2);
	SetPdeath(pdeath);
	SetSink(sink);
	SetAlpha(alpha);
	SetBeta(beta);
	SetTTExp(TTexp);
	SetInsertionAngle(InsertionAngle);
	SetDEGDMCOEF(degDMcoef);
}
//NNP
void ParamEntity::SetNNP(ParamEntity* Sigmas)
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
void ParamEntity::Random(int randUG,ParamEntity* orig,ParamEntity* Sigmas)
{
	int unsigned i;
	switch(randUG)
  {
     case 1://uniform distribution
		 for(i=0;i<this->params.size();i++)
		 {
			 this->params[i]=orig->params[i]+((Sigmas->params[i]==0)?(0):(Sigmas->params[i]*2*(genrand_real3()-0.5)));
		 }
		 this->SetParamEntity(params[0],params[1],params[2],params[3],params[4],params[5],params[6],
			 params[7],params[8],params[9],params[10],params[11],params[12],params[13],params[14]);
		 break;
	 case 0://Guauss distribution
		 for(i=0;i<this->params.size();i++)
		 {
			 this->params[i]=orig->params[i]+((Sigmas->params[i]==0)?(0):(Sigmas->params[i]*Moro_NormSInv(genrand_real3())));
		 }
		 this->SetParamEntity(params[0],params[1],params[2],params[3],params[4],params[5],params[6],
			 params[7],params[8],params[9],params[10],params[11],params[12],params[13],params[14]);
		 break;
	 
  }
}

void ParamEntity::Copy(ParamEntity* source)
{
	//Individual parameter values
	p=source->p;
	syntcoef=source->syntcoef;
	degcoef=source->degcoef;
	k1=source->k1;
	k2=source->k2;
	peff=source->peff;
	//pe2=source->pe2;
	//pe3=source->pe3;
	pm1=source->pm1;
	pm2=source->pm2;
	pdeath=source->pdeath;
	sink=source->sink;
	alpha=source->alpha;
	beta=source->beta;
	TTexp=source->TTexp;
	InsertionAngle=source->InsertionAngle;
	degDMcoef=source->degDMcoef;

	//Variable Initial values
	Nstructs0=source->Nstructs0;
	Nphs0=source->Nphs0;
	Lengths0=source->Lengths0;
	Areas0=source->Areas0;
	//Agreens0=source->Agreens0;
	TTinits0=source->TTinits0;
	DMrems0=source->DMrems0;
	DMstructs0=source->DMstructs0;

	//vector of parameter values
	params=source->params;
	//Vector of parameter names
	paramNames=source->paramNames;

	//WQL-SA index vetor of None_zero_sigmas parameters
	params_NNP=source->params_NNP;
	//std::vector<int> paramNames_NNP;
	
}
//WQL-SA

// END OF File : ParamLamina.cpp

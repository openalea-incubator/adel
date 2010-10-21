/**
 * @file ParamRoot.cpp
 * @brief Implementation of class ParamRoot
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamRoot.h"
#include "mt_generator.h"

// System includes ////////////////////////////////////////// System include //
#include <fstream>
#include <sstream>
#include <iostream>

#include <cassert>
// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
ParamRoot::ParamRoot()
{
}

ParamRoot::ParamRoot(std::string fileName)
//ParamRoot::ParamRoot(std::string fileName, std::string fileNameVar0) //ELMER

{
	//Individual parameter values
	ifstream file;
	file.open(fileName.c_str());
	string line;
	TTinit = - 1000;// en supposant que la floraison a lieu ?0 dd
	getline(file, line);
	Umax = GetFirstDouble(line);
	getline(file, line);
	kroot1 = GetFirstDouble(line);
	getline(file, line);
	kroot2 = GetFirstDouble(line);
	getline(file, line);
	QonDmin = GetFirstDouble(line);
	getline(file, line);
	beta1 = GetFirstDouble(line);
	getline(file, line);
	beta2 = GetFirstDouble(line);
	getline(file, line);
	Nsoilmin = GetFirstDouble(line);
	getline(file, line);
	p = GetFirstDouble(line);
	getline(file, line);
	sink = GetFirstDouble(line); 
	getline(file, line);
	alpha = GetFirstDouble(line); 
	getline(file, line);
	beta = GetFirstDouble(line); 
	getline(file, line);
	TTexp = GetFirstDouble(line); 
	getline(file, line);
	degDMcoef = GetFirstDouble(line); 
	getline(file, line);
	degcoef = GetFirstDouble(line); 
	file.close();


	//Vector of parameter values
	params.push_back(Umax);
	params.push_back(kroot1);
	params.push_back(kroot2);
	params.push_back(QonDmin);
	params.push_back(beta1);
	params.push_back(beta2);
	params.push_back(Nsoilmin);
	params.push_back(p);
	params.push_back(sink);
	params.push_back(alpha);
	params.push_back(beta);
	params.push_back(TTexp);
	params.push_back(degDMcoef);
	params.push_back(degcoef);
	params.push_back(TTinit);

	//Vector of parameter names
	paramNames.push_back("Umax");
	paramNames.push_back("kroot1");
	paramNames.push_back("kroot2");
	paramNames.push_back("QonDmin");
	paramNames.push_back("beta1");
	paramNames.push_back("beta2");
	paramNames.push_back("Nsoilmin");
	paramNames.push_back("p");
	paramNames.push_back("sink");
	paramNames.push_back("alpha");
	paramNames.push_back("beta");
	paramNames.push_back("TTexp");
	paramNames.push_back("degDMcoef");
	paramNames.push_back("degcoef");
	paramNames.push_back("TTinit");

	// Variable initial state
	// 3.75= MS totale= Cf donnï¿½es PM 1994 (traitement H2) * 0.15 = 0.56 g pour fraction allouÃ©e aux racines: Cf Siddique et al., 1990
	//(3g pour MS shoot) / 0.8 (de JL Drouet: 0.2= fraction MS totale allouï¿½e aux racines)
	
	DMremInit = 0.056;
	DMstructInit = 0.504;
	NremInit = 0.0011;
	NstructInit = 0.01;// still to be checked
	
//DMremInit=0.0616; //Trial 1
//DMstructInit=0.5544; //Trial 1
//NstructInit=0.009; //Trial 1
//NremInit=0.00099; //Trial 1


//  DMremInit=0.0504; //Trial 2
//  DMstructInit=0.4536; //Trial 2
// NstructInit=0.011; //Trial 2
// NremInit=0.00121; //Trial 2

//DMremInit=0.0616; //Trial 3
//DMstructInit=0.5544; //Trial 3
//NstructInit=0.011; //Trial 3
//NremInit=0.00121; //Trial 3

//DMremInit=0.0616; //Trial 4
//DMstructInit=0.5544; //Trial 4
//NstructInit=0.009; //Trial 4
//NremInit=0.00099; //Trial 4

//DMremInit=0.0504; //Trial 5
//DMstructInit=0.4536; //Trial 5
//NstructInit=0.009; //Trial 5
//NremInit=0.00099; //Trial 5

//DMremInit=0.0616; //Trial 6
//DMstructInit=0.5544; //Trial 6
//NstructInit=0.011; //Trial 6
//NremInit=0.00121; //Trial 6

//DMremInit=0.0504; //Trial 7
//DMstructInit=0.4536; //Trial 7
//NstructInit=0.011; //Trial 7
//NremInit=0.00121; //Trial 7

//DMremInit=0.0504; //Trial 8
//DMstructInit=0.4536; //Trial 8
//NstructInit=0.011; //Trial 8
//NremInit=0.00121; //Trial 8

//DMremInit=0.0616; //Trial 9
//DMstructInit=0.5544; //Trial 9
//NstructInit=0.011; //Trial 9
//NremInit=0.00121; //Trial 9

//DMremInit=0.0504; //Trial 10
//DMstructInit=0.4536; //Trial 10
//NstructInit=0.009; //Trial 10
//NremInit=0.00099; //Trial 10

//DMremInit=0.0616; //Trial 11
//DMstructInit=0.5544; //Trial 11
//NstructInit=0.009; //Trial 11
//NremInit=0.00099; //Trial 11

//DMremInit=0.0616; //Trial 12
//DMstructInit=0.5544; //Trial 12
//NstructInit=0.009; //Trial 12
//NremInit=0.00099; //Trial 12

//DMremInit=0.0504; //Trial 13
//DMstructInit=0.4536; //Trial 13
//NstructInit=0.009; //Trial 13
//NremInit=0.00099; //Trial 13

//DMremInit=0.0616; //Trial 14
//DMstructInit=0.5544; //Trial 14
//NstructInit=0.011; //Trial 14
//NremInit=0.00121; //Trial 14

//DMremInit=0.0504; //Trial 15
//DMstructInit=0.4536; //Trial 15
//NstructInit=0.011; //Trial 15
//NremInit=0.00121; //Trial 15

//DMremInit=0.0616; //Trial 16
//DMstructInit=0.5544; //Trial 16
//NstructInit=0.011; //Trial 16
//NremInit=0.00121; //Trial 16

//DMremInit=0.0504; //Trial 17
//DMstructInit=0.4536; //Trial 17
//NstructInit=0.011; //Trial 17
//NremInit=0.00121; //Trial 17

//DMremInit=0.0616; //Trial 18
//DMstructInit=0.5544; //Trial 18
//NstructInit=0.011; //Trial 18
//NremInit=0.00121; //Trial 18

//DMremInit=0.0504; //Trial 19
//DMstructInit=0.4536; //Trial 19
//NstructInit=0.011; //Trial 19
//NremInit=0.00121; //Trial 19

//DMremInit=0.0504; //Trial 20
//DMstructInit=0.4536; //Trial 20
//NstructInit=0.011; //Trial 20
//NremInit=0.00121; //Trial 20

//DMremInit=0.0616; //Trial 21
//DMstructInit=0.5544; //Trial 21
//NstructInit=0.011; //Trial 21
//NremInit=0.00121; //Trial 21

//DMremInit=0.0616; //Trial 22
//DMstructInit=0.5544; //Trial 22
//NstructInit=0.009; //Trial 22
//NremInit=0.00099; //Trial 22

//DMremInit=0.0504; //Trial 23
//DMstructInit=0.4536; //Trial 23
//NstructInit=0.009; //Trial 23
//NremInit=0.00099; //Trial 23

//DMremInit=0.0504; //Trial 24
//DMstructInit=0.4536; //Trial 24
//NstructInit=0.011; //Trial 24
//NremInit=0.00121; //Trial 24

//DMremInit=0.0616; //Trial 25
//DMstructInit=0.5544; //Trial 25
//NstructInit=0.011; //Trial 25
//NremInit=0.00121; //Trial 25

//DMremInit=0.0616; //Trial 26
//DMstructInit=0.5544; //Trial 26
//NstructInit=0.009; //Trial 26
//NremInit=0.00099; //Trial 26

//DMremInit=0.0504; //Trial 27
//DMstructInit=0.4536; //Trial 27
//NstructInit=0.009; //Trial 27
//NremInit=0.00099; //Trial 27

//DMremInit=0.0504; //Trial 28
//DMstructInit=0.4536; //Trial 28
//NstructInit=0.011; //Trial 28
//NremInit=0.00121; //Trial 28

//DMremInit=0.0616; //Trial 29
//DMstructInit=0.5544; //Trial 29
//NstructInit=0.011; //Trial 29
//NremInit=0.00121; //Trial 29

//DMremInit=0.0504; //Trial 30
//DMstructInit=0.4536; //Trial 30
//NstructInit=0.009; //Trial 30
//NremInit=0.00099; //Trial 30

//DMremInit=0.0616; //Trial 31
//DMstructInit=0.5544; //Trial 31
//NstructInit=0.009; //Trial 31
//NremInit=0.00099; //Trial 31

//DMremInit=0.0504; //Trial 32
//DMstructInit=0.4536; //Trial 32
//NstructInit=0.009; //Trial 32
//NremInit=0.00099; //Trial 32

//DMremInit=0.0616; //Trial 33
//DMstructInit=0.5544; //Trial 33
//NstructInit=0.009; //Trial 33
//NremInit=0.00099; //Trial 33

//DMremInit=0.0504; //Trial 34
//DMstructInit=0.4536; //Trial 34
//NstructInit=0.009; //Trial 34
//NremInit=0.00099; //Trial 34

//DMremInit=0.0616; //Trial 35
//DMstructInit=0.5544; //Trial 35
//NstructInit=0.009; //Trial 35
//NremInit=0.00099; //Trial 35

//DMremInit=0.0616; //Trial 36
//DMstructInit=0.5544; //Trial 36
//NstructInit=0.009; //Trial 36
//NremInit=0.00099; //Trial 36

//DMremInit=0.0504; //Trial 37
//DMstructInit=0.4536; //Trial 37
//NstructInit=0.009; //Trial 37
//NremInit=0.00099; //Trial 37

//DMremInit=0.0616; //Trial 38
//DMstructInit=0.5544; //Trial 38
//NstructInit=0.011; //Trial 38
//NremInit=0.00121; //Trial 38

//DMremInit=0.0504; //Trial 39
//DMstructInit=0.4536; //Trial 39
//NstructInit=0.011; //Trial 39
//NremInit=0.00121; //Trial 39

//DMremInit=0.0504; //Trial 40
//DMstructInit=0.4536; //Trial 40
//NstructInit=0.009; //Trial 40
//NremInit=0.00099; //Trial 40


	

	
	//ELMER
	//ifstream fileVar0;
	//fileVar0.open(fileNameVar0.c_str());
	//string lineVar0;
	//getline(fileVar0, lineVar0);
	//DMremInit = GetFirstDouble(lineVar0);
	//getline(fileVar0, lineVar0);
	//DMstructInit = GetFirstDouble(lineVar0);
	//getline(fileVar0, lineVar0);
	//NremInit = GetFirstDouble(lineVar0);
	//getline(fileVar0, lineVar0);
	//NstructInit = GetFirstDouble(lineVar0);
	//fileVar0.close();
	//ELMER
	
}

ParamRoot::~ParamRoot()
{
}


////////////////////////////// Private Methods ////////////////////////////////
double ParamRoot::GetFirstDouble(string line)
{
	stringstream oss;
	oss << line;
	double tempValue;
	oss >> tempValue;
	return tempValue;
}

/////////////////////////////// Public Methods ////////////////////////////////

//Individual parameter values
double ParamRoot::GetUMAX()
{
	return Umax;
}
double ParamRoot::GetKroot1()
{
	return kroot1;
}
double ParamRoot::GetKroot2()
{
	return kroot2;
}
double ParamRoot::GetQONDMIN()
{
	return QonDmin;
}
double ParamRoot::GetBETA1()
{
	return beta1;
}
double ParamRoot::GetBETA2()
{
	return beta2;
}
double ParamRoot::GetNsoilMin()
{
	return Nsoilmin;
}
double ParamRoot::GetP()
{
	return p;
}

double ParamRoot::GetSink()
{
	return sink;
}
double ParamRoot::GetAlpha()
{
	return alpha;
}
double ParamRoot::GetBeta()
{
	return beta;
}
double ParamRoot::GetTTExp()
{
	return TTexp;
}
double ParamRoot::GetDEGDMCOEF()
{
	return degDMcoef;
}
double ParamRoot::GetDEGCOEF()
{
	return degcoef;
}
double ParamRoot::GetTTinit()
{
	return TTinit;
}


// Vectors of parameter values and parameter names:
std::vector<double> ParamRoot::GetParams()
{
	// voir si je dois mettre en attribut ce vecteur de paramètres
	return params;
}
std::vector<string> ParamRoot::GetParamNames()
{
	return paramNames;
}
double ParamRoot::GetDMremInit()
{
	return DMremInit;
}
double ParamRoot::GetDMstructInit()
{
	return DMstructInit;
}
double ParamRoot::GetNremInit()
{
	return NremInit;
}
double ParamRoot::GetNstructInit()
{
	return NstructInit;
}
void ParamRoot::SetDMremInit(double value)
{
	DMremInit=value;
}
void ParamRoot::SetNremInit(double value)
{
	NremInit=value;
}
void ParamRoot::SetDMstructInit(double value)
{
	DMstructInit=value;
}
void ParamRoot::SetNstructInit(double value)
{
	NstructInit=value;
}

//WQL-SA
#define DEFINE_SET(_METH,_PAR,_INDEX)		\
void ParamRoot::_METH(double value)\
{\
	_PAR=value;\
        params[_INDEX] = _PAR;\
        assert(paramNames[ _INDEX ] == #_PAR);\
}

DEFINE_SET(SetUMAX,Umax,0)
DEFINE_SET(SetKroot1,kroot1,1)
DEFINE_SET(SetKroot2,kroot2,2)
DEFINE_SET(SetQONDMIN,QonDmin,3)
DEFINE_SET(SetBETA1,beta1,4)
DEFINE_SET(SetBETA2,beta2,5)
DEFINE_SET(SetNsoilMin,Nsoilmin,6)
DEFINE_SET(SetP,p,7)
DEFINE_SET(SetSink,sink,8)
DEFINE_SET(SetAlpha,alpha,9)
DEFINE_SET(SetBeta,beta,10)
DEFINE_SET(SetTTExp,TTexp,11)
DEFINE_SET(SetDEGDMCOEF,degDMcoef,12)
DEFINE_SET(SetDEGCOEF,degcoef,13)
DEFINE_SET(SetTTinit,TTinit,14)

void ParamRoot::SetParamRoot(double Umax,double kroot1,double kroot2,double QonDmin,double beta1,double beta2,
			     double Nsoilmin,double p,double sink,double alpha,double beta,double TTexp, double degDMcoef, double degcoef,double TTinit)
{
	SetUMAX(Umax);
	SetKroot1(kroot1);
	SetKroot2(kroot2);
	SetQONDMIN(QonDmin);
	SetBETA1(beta1);
	SetBETA2(beta2);
	SetNsoilMin(Nsoilmin);
	SetP(p);
	SetSink(sink);
	SetAlpha(alpha);
	SetBeta(beta);
	SetTTExp(TTexp);
	SetDEGDMCOEF(degDMcoef);
	SetDEGCOEF(degcoef);
	SetTTinit(TTinit);
}
//NNP
void ParamRoot::SetNNP(ParamRoot* Sigmas)
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
void ParamRoot::Random(int randUG,ParamRoot* orig,ParamRoot* Sigmas)
{
	int unsigned i;
	switch(randUG)
  {
     case 1://uniform distribution
		 for(i=0;i<this->params.size();i++)
		 {
			 this->params[i]=orig->params[i]+((Sigmas->params[i]==0)?(0):(Sigmas->params[i]*2*(genrand_real3()-0.5)));
		 }
		 this->SetParamRoot(params[0],params[1],params[2],params[3],params[4],params[5],params[6],
			 params[7],params[8],params[9],params[10],params[11],params[12],params[13],params[14]);
		 break;
	 case 0://Guauss distribution
		 for(i=0;i<this->params.size();i++)
		 {
			 this->params[i]=orig->params[i]+((Sigmas->params[i]==0)?(0):(Sigmas->params[i]*Moro_NormSInv(genrand_real3())));
		 }
		 this->SetParamRoot(params[0],params[1],params[2],params[3],params[4],params[5],params[6],
			 params[7],params[8],params[9],params[10],params[11],params[12],params[13],params[14]);
		 break;
	 
  }
}

void ParamRoot::Copy(ParamRoot* source)
{
	//Individual parameter values
	Umax=source->Umax;
	kroot1=source->kroot1;
	kroot2=source->kroot2;
	QonDmin=source->QonDmin;
	beta1=source->beta1;
	beta2=source->beta2;
	Nsoilmin=source->Nsoilmin;
	p=source->p;
	sink=source->sink;
	alpha=source->alpha;
	beta=source->beta;
	TTexp=source->TTexp;
	degDMcoef=source->degDMcoef;
	degcoef=source->degcoef;
	TTinit = source->TTinit;

	//Vector of parameter values
	params=source->params; 
	//Vector of parameter names
	paramNames=source->paramNames;

	//WQL-SA index vetor of None_zero_sigmas parameters
	params_NNP=source->params_NNP;
	//std::vector<int> paramNames_NNP;
	
}
//WQL-SA
// END OF File : ParamRoot.cpp

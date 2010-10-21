/**
 * @file ParamEnv.h
 * @brief Declaration of class ParamEnv
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifndef PARAMENV_H_
#define PARAMENV_H_

// Local includes /////////////////////////////////////////// Local Includes //

// System includes ////////////////////////////////////////// System include //
#include <string>
#include <vector>

class ParamEnv
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	ParamEnv();
	ParamEnv(std::string fileNameParam,std::string fileNameN, std::string fileNameT, std::string fileNameMeteo);
    
	virtual ~ParamEnv();

// Public Methods /////////////////////////////////////////// Public Methods //
public:
	//Individual parameter values
	double GetTransmission();
	double GetCoef_Ext();// for vertical entities

	//Vector of parameter values and names
	std::vector<double> GetParams();
	std::vector<std::string> GetParamNames();

	//Environmental variables values through time
	std::vector<double> GetNsoils();
	std::vector<double> GetTimeSoils();
	std::vector<double> GetDays();
	std::vector<double> GetPAR0s();
	std::vector<double> GetTres();

	//WQL-SA
	void SetTransmission(double value);
	void SetCoef_Ext(double value);
	void SetParamEnv(double Transmission,double Coef_Ext);
	//NNP
	void SetNNP(ParamEnv* Sigmas);
	//random
	void Random(int randUG,ParamEnv* orig,ParamEnv* Sigmas);
	//Copy
	void Copy(ParamEnv* source);
	//WQL-SA

//ADD ELMER
	//int GetNbRow();
//ADD ELMER



// Private Methods ///////////////////////////////////////// Private Methods //
private:
	double GetFirstDouble(std::string line);
	

// Attributes /////////////////////////////////////////////////// Attributes //
private:
	//Individual parameter values
	double transmission;
	double coef_ext;


	//Environmental variables values through time
	std::vector<double> Nsoils;
	std::vector<double> TimeSoils;
	std::vector<double> Days;
	std::vector<double> PAR0s;
	std::vector<double> Tres;

public:
	//vector of parameter values
	std::vector<double> params;

	//Vector of parameter names
	std::vector<std::string> paramNames;

	//WQL-SA index vetor of None_zero_sigmas parameters
    std::vector<int> params_NNP;
	//std::vector<int> paramNames_NNP;
	//WQL-SA



};

#endif /* PARAMENV_H_ */

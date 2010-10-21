/**
 * @file ParamGrain.h
 * @brief Declaration of class ParamGrain
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifndef PARAMGRAIN_H_
#define PARAMGRAIN_H_

// Local includes /////////////////////////////////////////// Local Includes //

// System includes ////////////////////////////////////////// System include //
#include <string>
#include <vector>

class ParamGrain
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	ParamGrain();
	ParamGrain(std::string fileNameParam, std::string fileNameVar0);
    virtual ~ParamGrain();

// Public Methods /////////////////////////////////////////// Public Methods //
public:
	//Individual parameter values
	double GetNGrainInit();
	double GetDMGrainInit();
	double GetP();
	double GetTT();
	double GetG();
	double GetSink();
	double GetAlpha();
	double GetBeta();
	double GetTTExp();
	double GetTTinit();

	//Vector of parameter values and names
	std::vector<double> GetParams();
	std::vector<std::string> GetParamNames();

	//WQL-SA
        void SetNGrainInit(double value);
        void SetDMGrainInit(double value);
	void SetP(double value);
	void SetTT(double value);
	void SetG(double value);
	void SetSink(double value);
	void SetAlpha(double value);
	void SetBeta(double value);
	void SetTTExp(double value);
	void SetTTinit(double value);
	void SetParamGrain(double p,double tt,double g,double sink,double alpha,double beta,double TTexp,double TTinit);

	//NNP
	void SetNNP(ParamGrain* Sigmas);
	//random
	void Random(int randUG,ParamGrain* orig,ParamGrain* Sigmas);
	//copy
	void Copy(ParamGrain* source);
	//WQL-SA

// Private Methods ///////////////////////////////////////// Private Methods //
private:
	double GetFirstDouble(std::string line);
	

// Attributes /////////////////////////////////////////////////// Attributes //
private:
	// parametres d'initialisation
	double NGrainInit;
	double DMGrainInit;
	//Individual parameter values
	double TTinit;
	double p;
	double tt;
	double g;
	double sink;
	double alpha;
	double beta;
	double TTexp;

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

#endif /* PARAMGRAIN_H_ */

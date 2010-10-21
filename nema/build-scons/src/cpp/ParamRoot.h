/**
 * @file ParamRoot.h
 * @brief Declaration of class ParamRoot
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifndef PARAMROOT_H_
#define PARAMROOT_H_

// Local includes /////////////////////////////////////////// Local Includes //

// System includes ////////////////////////////////////////// System include //
#include <string>
#include <vector>

class ParamRoot
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	ParamRoot();
    
	
	ParamRoot(std::string fileName);
	
	//ParamRoot(std::string fileName, std::string fileNameVar0); //ELMER
	
    virtual ~ParamRoot();

// Public Methods /////////////////////////////////////////// Public Methods //
public:
	//Individual parameter values
	double GetUMAX();
	double GetKroot1();
	double GetKroot2();
	double GetQONDMIN();
	double GetBETA1();
	double GetBETA2();
	double GetNsoilMin();
	double GetP();
	double GetSink();
	double GetAlpha();
	double GetBeta();
	double GetTTExp();
	double GetDEGDMCOEF();
	double GetDEGCOEF();
	double GetTTinit();

	//initial values
	double GetDMremInit();
	double GetDMstructInit();
	double GetNremInit();
	double GetNstructInit();
	//Vector of parameter values and names
	std::vector<double> GetParams();
	std::vector<std::string> GetParamNames();

	//WQL-SA
	void SetUMAX(double value);
	void SetKroot1(double value);
	void SetKroot2(double value);
	void SetQONDMIN(double value);
	void SetBETA1(double value);
	void SetBETA2(double value);
	void SetNsoilMin(double value);
	void SetP(double value);
	void SetSink(double value);
	void SetAlpha(double value);
	void SetBeta(double value);
	void SetTTExp(double value);
	void SetDEGDMCOEF(double value);
	void SetDEGCOEF(double value);
	void SetTTinit(double value);
	void SetParamRoot(double Umax,double kroot1,double kroot2,double QonDmin,double beta1,double beta2,
			  double Nsoilmin,double p,double sink,double alpha,double beta,double TTexp,double degDMcoef,double degcoef,double TTinit);
	
	//Set initials
	void SetDMremInit(double value);
	void SetDMstructInit(double value);
	void SetNremInit(double value);
	void SetNstructInit(double value);

	//NNP
	void SetNNP(ParamRoot* Sigmas);
	//random
	void Random(int randUG,ParamRoot* orig,ParamRoot* Sigmas);
	//copy
	void Copy(ParamRoot* source);
//WQL-SA

// Private Methods ///////////////////////////////////////// Private Methods //
private:
	double GetFirstDouble(std::string line); //?ne mettre qu'une fois dans tout le programme par la suite

// Attributes /////////////////////////////////////////////////// Attributes //
private:
	//parametre d'initialisation
	double DMremInit;
	double DMstructInit;
	double NstructInit;
	double NremInit;
	//Individual parameter values
	double Umax;
	double kroot1;
	double kroot2;
	double QonDmin;
	double beta1;
	double beta2;
	double Nsoilmin;
	double p;
	double sink;
	double alpha;
	double beta;
	double TTexp;
	double degDMcoef;
	double degcoef;
	double TTinit;

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

#endif /* PARAMROOT_H_ */

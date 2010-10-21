/**
 * @file ParamLamina.h
 * @brief Declaration of class ParamLamina
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifndef PARAMENTITY_H_
#define PARAMENTITY_H_

// Local includes /////////////////////////////////////////// Local Includes //

// System includes ////////////////////////////////////////// System include //
#include <string>
#include <vector>

class ParamEntity
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	ParamEntity();
	ParamEntity(std::string fileNameParam, std::string fileNameVar0, int nrow);
    virtual ~ParamEntity();

// Public Methods /////////////////////////////////////////// Public Methods //
public:
	// Add a new Entity
	void AddEntity(double Nstruct, double Nph, double Length, 
			double Area, double TTinit, double DMrem, double DMstruct);
	//Individual parameter values
	double GetP();
	double GetSYNTCOEF();
	double GetDEGCOEF();
	double GetK1();
	double GetK2();
	double GetPeff();
	//double GetPe2();
	//double GetPe3();
	double GetPm1();
	double GetPm2();
	double GetPdeath();
	double GetSink();
	double GetAlpha();
	double GetBeta();
	double GetTTExp();
	double GetInsertionAngle();
	double GetDEGDMCOEF();

	//Vector of parameter values and names
	std::vector<double> GetParams();
	std::vector<std::string> GetParamNames();

	//Variable Initial values
	std::vector<double> GetNstructs0();
	std::vector<double> GetNphs0();
	std::vector<double> GetLengths0();
	std::vector<double> GetAreas0();
	//std::vector<double> GetAgreens0();
	std::vector<double> GetTTinits0();
	std::vector<double> GetDMrems0();
	std::vector<double> GetDMstructs0();

	void SetNstructs0(const std::vector<double>& values);
	void SetNphs0(const std::vector<double>& values);
	void SetLengths0(const std::vector<double>& values);
	void SetAreas0(const std::vector<double>& values);
	//std::vector<double> GetAgreens0();
	void SetTTinits0(const std::vector<double>& values);
	void SetDMrems0(const std::vector<double>& values);
	void SetDMstructs0(const std::vector<double>& values);
//WQL-SA
	void SetP(double value);
	void SetSYNTCOEF(double value);
	void SetDEGCOEF(double value);
	void SetK1(double value);
	void SetK2(double value);
	void SetPeff(double value);
	//void SetPe2(double value);
	//void SetPe3(double value);
	void SetPm1(double value);
	void SetPm2(double value);
	void SetPdeath(double value);
	void SetSink(double value);
	void SetAlpha(double value);
	void SetBeta(double value);
	void SetTTExp(double value);
	void SetInsertionAngle(double value);
	void SetDEGDMCOEF(double value);
	void SetParamEntity(double p,double syntcoef,double degcoef,double k1,double k2,double peff,double pm1,double pm2,double pdeath,double sink,double alpha,double beta,double TTexp,double InsertionAngle, double degDMcoef);
	//NNP
	void SetNNP(ParamEntity* Sigmas);
	//random
	void Random(int randUG,ParamEntity* orig,ParamEntity* Sigmas);
	//copy
	void Copy(ParamEntity* source);
//WQL-SA

//ADD ELMER
	int GetNbRow();
	//std::vector<std::vector<double> > GetNstructs();
//ADD ELMER



// Private Methods ///////////////////////////////////////// Private Methods //
private:
	double GetFirstDouble(std::string line);


// Attributes /////////////////////////////////////////////////// Attributes //
private:
	//Individual parameter values
	double p;
	double syntcoef;
	double degcoef;
	double k1;
	double k2;
	double peff;
	//double pe2;
	//double pe3;
	double pm1;
	double pm2;
	double pdeath;
	double sink;
	double alpha;
	double beta;
	double TTexp;
	double InsertionAngle;
	double degDMcoef;

	//Variable Initial values
	std::vector<double> Nstructs0;
	std::vector<double> Nphs0;
	std::vector<double> Lengths0;
	std::vector<double> Areas0;
	//std::vector<double> Agreens0;
	std::vector<double> TTinits0;
	std::vector<double> DMrems0;
	std::vector<double> DMstructs0;
protected:
	// Number of row
	int nbrow;

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

#endif /* PARAMENTITY_H_ */

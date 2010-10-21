/**
 * @file Environment.h
 * @brief Declaration of class Environment
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamEnv.h"

// System includes ////////////////////////////////////////// System include //
#include <string>
#include <vector>

class Environment:public ParamEnv
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	Environment();
    Environment(std::string fileNameParam, std::string fileNameN, std::string fileNameT, std::string fileNameMeteo);
    virtual ~Environment();

// Public Methods /////////////////////////////////////////// Public Methods //
public:
	double GetNsoil(int Time);
	double ThermalTime(int time, double tre_base, int step, int nbsubsteps);
	std::vector<double> PARs(int time, int nbsubsteps, int type, int nrow, std::vector<int> rows_exposed, std::vector<std::vector<double> > gai_cums, double InsertionAngleLamina, std::vector<double> HeightVis);
	//void ActuPARs(int ntype, int nrow, double gai_cum);
	//std::vector<std::vector<double> > GetPARss();
	//std::vector<double> GetPARs(int type);

// Private Methods ///////////////////////////////////////// Private Methods //
private:
	double InterpolLin(int t, std::vector<double> fileTps, std::vector<double> fileVal);
	double PAR0(int time,int nbsubsteps);
	double Tre(int time);

// Attributes /////////////////////////////////////////////////// Attributes //
//private:
	//std::vector<double> GAICUMs;
	//double PAR0;
	//ParamEnv* parameters;
	//std::vector<double> PARs;
	//std::vector<std::vector<double> > PARss;
};

#endif /* ENVIRONMENT_H_ */

/**
 * @file Grain.h
 * @brief Declaration of class Grain
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifndef GRAIN_H_
#define GRAIN_H_

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamGrain.h"
// System includes ////////////////////////////////////////// System include //
#include <vector>
#include <string>

class Grain:public ParamGrain
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	Grain();
    Grain(std::string fileNameParam, std::string fileNameVar0,std::string outputdir);
    virtual ~Grain();

// Public Methods /////////////////////////////////////////// Public Methods //
public:
	void InitialGrain();
	void ResetGrain();
	void ActuGrain(double tt, double dtt, double dt, double Ncmob, double NmobMin, double Q, double D, double ConcNmob);
	double ExportGrain(double tt, double dtt, double Ncmob, double NmobMin, double Q, double D, double ConcNmob);
	//double ExportGrain_RK(double dt, double tt, double dtt, double Ncmob, double NmobMin);
	double ComputeDemand(double TT);
	void ActuDryMass(double TT, double Q, double D);
	double GetDMgrain();
	double GetNgrain();
	void Affich(std::vector<double> TTimes,std::string outputdir);

//ADD ELMER
	//int GetNbRow();
	std::vector<double> GetNgrains();
	std::vector<double> GetNPots();
	std::vector<double> GetDMs();
//ADD ELMER


// Private Methods ///////////////////////////////////////// Private Methods //
private:
	double ExportGrainPot(double tt, double dtt, double Npot);
	double NSynthRate(double ConcNmob,double dDM);
	//double ExportGrainPot_RK(double dt, double tt, double dtt);

// Attributes /////////////////////////////////////////////////// Attributes //
private:
	double Ngrain;
	std::vector <double> Ngrains;
	double NPot;
	std::vector <double> NPots;
	double TTinit; // pourquoi la mettre en attribut ? Est-ce vraiment nécessaire ?
	double NormSink;
	double DM;
	std::vector<double> DMs;
};

#endif /* GRAIN_H_ */

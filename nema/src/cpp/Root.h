/**
 * @file Root.h
 * @brief Declaration of class Root
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifndef ROOT_H_
#define ROOT_H_

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamRoot.h"

// System includes ////////////////////////////////////////// System include //
#include <vector>
#include<string>

class Root:public ParamRoot
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	Root();
	
	
	Root(std::string fileName,std::string outputdir);
	//Root(std::string fileName, std::string fileNameVar0, std::string outputdir); //ELMER
	
	
    virtual ~Root();

// Public Methods /////////////////////////////////////////// Public Methods //
public:
	void InitialRoot();
	void ResetRoot();
	double Ceffect(double Cstatus, int nbsubsteps);
	double Neffect(double Nstatus);
	double ImportRootPot(double dTT, double Nsoil);
	double ImportRoot(double dTT, int nbsubsteps, double Nstatus,double Cstatus, double Nsoil);
	double Ninflux(double ConcNmob,double ConcNmobMin,double TT, double dTT, double Q, double D);
	void ActuNrem(double ConcNmob, double ConcNmobMin, double TT, double dTT, double Q, double D);
	double ComputeDemand(double TT);
	void ActuDMrem(double TT, double dTT, double Q, double D);
	double ComputeRemobilizedDM(double dTT, double TT);
	void Affich(std::vector<double> TTimes,std::vector<double> dTTs, std::string outputdir);

//ADD ELMER
	std::vector<double> GetNstructs();
	std::vector<double> GetNrems();
	std::vector<double> GetDMstructs();
	std::vector<double> GetDMrems();
//ADD ELMER



// Private Methods ///////////////////////////////////////// Private Methods //
private:
	double NremSynthRate(double ConcNmob, double dDM); //voir si je dois mettre ConcNmobMin
	double DMremDegRate(double dTT, double dm);
	double NremDegRate(double dTT, double N);

// Attributes /////////////////////////////////////////////////// Attributes //
private:
	//ParamRoot* parameters;
	double Nstruct;
	double Nrem;
	std::vector<double> Nstructs;
	std::vector<double> Nrems;
	double TTinit; // pourquoi la mettre en attribut ? Est-ce vraiment nécessaire ?
	double NormSink;
	double DMstruct;
	double DMrem;
	std::vector<double> DMstructs;
	std::vector<double> DMrems;
};

#endif /* ROOT_H_ */

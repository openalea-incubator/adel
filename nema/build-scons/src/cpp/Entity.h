/**
 * @file Entity.h
 * @brief Declaration of class Entity
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifndef ENTITY_H_
#define ENTITY_H_

// Local includes /////////////////////////////////////////// Local Includes //
#include "ParamEntity.h"

// System includes ////////////////////////////////////////// System include //
#include <vector>
#include <string>
//using namespace std;

class Entity:public ParamEntity
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	Entity();
	Entity(std::string fileNameParam, std::string fileNameVar0,int nrow,std::string outputdir);
    virtual ~Entity();

// Public Methods /////////////////////////////////////////// Public Methods //
public:
	void InitialEntity();
	void ResetEntity();
	double DMGreen(int row);
	double VisibleLength(int type, int row, double Height, double HeightPreviousEntity);
	double VisibleArea(int row, double VisibleLength);
	void InitializeAreaGreensNm2Max(int nrow, std::vector<double> visibleAreas);
	//void ActuVisibleAreas(int row, double height, double heightPreviousOrgan);
	/*double NphSynthRate(int row,double dTT,double ConcNmob, double PAR);
	double NphDegRate(int row,double dTT);*/
	double NphInput(int row, double TT, double dTT, double ConcNmob, double ConcNmobMin, double PAR, double Q, double D);
	void ActuNstructs(int nrow, double ConcNmob, std::vector<double> PARs, double TT, double Q, double D);
	void ActuNphs(int nrow, double TT, double dTT, double ConcNmob, double ConcNmobMin, std::vector<double> PARs, double Q, double D);
	void ActuAreaGreens(int nrow, double TT, double dTT, double ConcNmob, double ConcNmobMin, std::vector<double> PARs, double Q, double D);
	void ActuNphM2Max(int nrow, double TT, double dTT, double ConcNmob, double ConcNmobMin, std::vector<double> PARs, double Q, double D);
	void ActuDryMass(int nrow, double TT, double dTT, double Q, double D);
	double ComputeProduction(int nrow, std::vector<double> PARs, int T, int nbsubsteps, double TT);
	double ComputeRemobilizedDM(int nrow, double dTT, double TT);
	double ComputeDemand(int nrow, double TT);
	//double GetAreaGreen0(int row);
	double GetAreaGreenTot(int nrow);
	void Affich(std::vector<double> dTTs,std::vector<double> TT,std::string type,int nrow,std::string outputdir);

	int GetNbRow();
	//Add Elmer
	std::vector<std::vector<double> > GetNstructss();
	std::vector<std::vector<double> > GetNphss();
	std::vector<std::vector<double> > GetTTinitss();
	std::vector<std::vector<double> > GetDMremss();
	std::vector<std::vector<double> > GetDMstructss();
	std::vector<std::vector<double> > GetAreaGreenss();
	std::vector<std::vector<double> > GetNphM2Maxss();
	std::vector<std::vector<double> > GetNphSynthRatess();
	std::vector<std::vector<double> > GetNphloemSynthRatess();
	std::vector<std::vector<double> > GetDMremSynthRatess();
	std::vector<std::vector<double> > GetPhotoss();
	//Add Elmer

// Private Methods ///////////////////////////////////////// Private Methods //
private:
	double NstructSynthRate(double ConcNmob,double photo, double dDM); //voir si je dois mettre ConcNmobMin
	double NphSynthRate(double dTT,double ConcNmob, double DMgreen, double PAR);
	double NphDegRate(double dTT,double Nph);
	double NphM2Min(double Nphm2max, int row);
	double DeltaAreaGreen(double Areagreen, double Nph, double Nphinput, double Nphm2max, int row);
	//double Peff(double Nph, double Areagreen);
	double Pmax(double Nph, double Areagreen);
	double Photo(double Nph, double Areagreen, double PAR);
	double ComputeSink(double TT, double TTin);
	double DMDegRate(double dTT, double dm);
	

// Attributes /////////////////////////////////////////////////// Attributes //
#define LAMINA 1 // ?mettre ici ?
#define SHEATH 2
#define INTERNODE 3 
#define PEDUNCLE 4
#define CHAFF 5

private:
	//int nbrow;
	std::vector<double> Nstructs;
	std::vector<std::vector<double> > Nstructss;
	std::vector<double> Nphs;
	std::vector<std::vector<double> > Nphss;
	std::vector<double> Area0s;// Area at entity maturity
	std::vector<double> TTinits;
	std::vector<std::vector<double> > TTinitss;
	double NormSink;
	std::vector<double> DMrems;
	std::vector<std::vector<double> > DMremss;
	std::vector<double> DMstructs;
	std::vector<std::vector<double> > DMstructss;

	std::vector<double> Lengths;
	std::vector<double> AreaGreens;
	std::vector<std::vector<double> > AreaGreenss;
	std::vector<double> NphM2Maxs;
	std::vector<std::vector<double> > NphM2Maxss;

	std::vector<std::vector<double> > NphSynthRatess;
	std::vector<std::vector<double> > NphloemSynthRatess;
	std::vector<std::vector<double> > DMremSynthRatess;

	std::vector<std::vector<double> > Photoss;

};

#endif /* ENTITY_H_ */

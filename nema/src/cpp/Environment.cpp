/**
 * @file Environment.cpp
 * @brief Implementation of class Environment
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "Environment.h"
#include "math.h"

// System includes ////////////////////////////////////////// System include //
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
Environment::Environment()
{
}
Environment::Environment(std::string fileNameParam, std::string fileNameN, std::string fileNameT, std::string fileNameMeteo):ParamEnv(fileNameParam, fileNameN, fileNameT, fileNameMeteo)
{
	//this->parameters = parameters;
//	PAR0=10368000;

	//std::vector<double> PARs1,PARs2;
//	PARs.push_back(237600);// trouver une autre solution par la suite que de donner ces valeurs par d�fault ! + pour les gaines; ce ne sont pas les bonnes valeurs
//	PARs.push_back(95040);
//	PARs.push_back(21600);
//	PARs.push_back(12960);
//	PARs2.push_back(100000); 
//	PARs2.push_back(30000);
//	PARs2.push_back(15000);
//	PARs2.push_back(10000);
//	PARss.push_back(PARs1);
//	PARss.push_back(PARs2);
}

Environment::~Environment()
{
}

////////////////////////////// Private Methods ////////////////////////////////

double Environment::InterpolLin(int t, std::vector<double> fileTps, std::vector<double> fileVal)
{
  double r=0;
  int ind=-1;
  int NBPT=fileTps.size();
  for (int k=0;k<NBPT;k++)
  {
    if (t<=fileTps[k] && ind==-1)
	{
      ind=k;
	}
  }
  if (t==fileTps[ind])
  {
    r=fileVal[ind];
  }
  else 
  {
	double slope=(fileVal[ind]-fileVal[ind-1])/(fileTps[ind]-fileTps[ind-1]);
	double intercept=fileVal[ind-1]-slope*fileTps[ind-1];
	r=slope*t+intercept;
  }
  return r;
}



double Environment::PAR0(int time, int nbsubsteps) // pas besoin de mettre nbsubsteps en réalité car (i) la production de chaque entité est déjà divisée par 10 à chaque pas de temps; (ii) d'exprimer le PAR en jour ou jour/10 n'affecte pas le taux de synthèse de Nph.
{
	int ind=-1;
	for (unsigned int i=0;i<(GetDays()).size();i++) 
	{
		if (GetDays()[i]==time)
		{
			ind=i;
		}
	}
	if (ind >=0)
	{
		return GetPAR0s()[ind];
		cout << "The days of the meteo file do not correspond to the model days" << endl;
	}
	return 0.;
}

double Environment::Tre(int time)
{
	int ind=-1;
	for (unsigned int i=0;i<(GetDays()).size();i++) 
	{
		//cout << GetDays()[i] << "\t" << time << endl;
		if (GetDays()[i]==time)
		{
			ind=i;
		}
		//cout << i << endl;
	}
	if (ind >= 0)
	{
		return GetTres()[ind];
	}else
	{
		cout << "The days of the meteo file do not correspond to the model days" << endl;
	}
	return 0.;
}

/////////////////////////////// Public Methods ////////////////////////////////
double Environment::GetNsoil(int Time)
{
	return InterpolLin(Time,GetTimeSoils(),GetNsoils());
}

double Environment::ThermalTime(int time, double tre_base, int step, int nbsubsteps)
{
	if (Tre(time)>tre_base)
	{
		return (Tre(time)- tre_base)*step/nbsubsteps;
	}else
	{
		return 0;
	}
}


std::vector<double> Environment::PARs(int time, int nbsubsteps, int type, int nrow, std::vector<int> rows_exposed, std::vector<std::vector<double> > gai_cumss, double InsertionAngleLamina, std::vector<double> HeightVis) // pour l'instant, vu qu'il n'y a que des feuilles, c'est du lai: changer la valeur de k pour avoir du gai
{
	std::vector<double> PARs;
	for (int r=0;r<nrow;r++)
	{
		int OrganPositionInGAIs=-1;
		for (unsigned int ri=0;ri<rows_exposed.size();ri++)
		{
			if ( (r+1)==rows_exposed[ri] )
			{
				OrganPositionInGAIs=ri;
			}
		}
		if (OrganPositionInGAIs>(-1))
		{
			// modifier l'équation suivante pour la rendre plus robuste (CF je suis obligée de mettre les numéros de type d'organes)
			double Coef_ExtLamina=GetCoef_Ext()+ (1-GetCoef_Ext())*cos(InsertionAngleLamina); // arbitrary equation for the variation of the extinction coefficient of the laminae in function of their angle;
			// GetCoef_Ext() est le coefficient d'extinction pour les entités verticales
			double Ext_Canopy = gai_cumss[0][OrganPositionInGAIs]* Coef_ExtLamina + gai_cumss[1][OrganPositionInGAIs]*GetCoef_Ext() + gai_cumss[2][OrganPositionInGAIs]*GetCoef_Ext() + gai_cumss[3][OrganPositionInGAIs]*GetCoef_Ext() + gai_cumss[4][OrganPositionInGAIs]*GetCoef_Ext();
			//cout << Ext_Canopy << endl;
			double GAICUM=gai_cumss[0][OrganPositionInGAIs]+gai_cumss[1][OrganPositionInGAIs]+gai_cumss[2][OrganPositionInGAIs]+gai_cumss[3][OrganPositionInGAIs]+gai_cumss[4][OrganPositionInGAIs];
			//GAICUMs.push_back(GAICUM);
			
			double Coef_Ext;
			if (type==1)
			{
				Coef_Ext=Coef_ExtLamina;
			}else
			{
				Coef_Ext=GetCoef_Ext(); // arbitrary value for vertical organs
			}
			//cout << time << "\t" << type <<  "\t" << Coef_Ext << endl;
			double PAR;
			//if (type!=5)
			//{
			//	PAR= PAR0(time,nbsubsteps) * Coef_Ext * exp( (-1) * Ext_Canopy );//il faudrait faire quelque chose de plus robuste !!
			//}else
			//{
			PAR= PAR0(time,nbsubsteps) * Coef_Ext * exp( (-1) * Ext_Canopy ); // à voir, ce n'est pas propre ...
			//cout << PAR << "\t" << type << "\t" << Ext_Canopy << "\t" << GAICUM << endl;
			//}
			PARs.push_back(PAR);
			
		}else
		{
			PARs.push_back(0);
		}

	}
	return PARs;

}

//double Environment::GetPAR(int row)
//{
//	return PARs[row-1];
//}


//std::vector<std::vector<double> > Environment::GetPARss()
//{
//	return PARss;
//}

//std::vector<double> Environment::GetPARs(int type)
//{
//	return PARss[type-1];
//}


// END OF File : Environment.cpp

/**
 * @file Grain.cpp
 * @brief Implementation of class Grain
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "Grain.h"

// System includes ////////////////////////////////////////// System include //
#include <iostream>
#include <fstream>
#include "math.h"
#include <string>

// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
Grain::Grain()
{
//vector of parameter values
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);

	//Vector of parameter names
	paramNames.push_back("p");
	paramNames.push_back("tgrain");
	paramNames.push_back("ggrain");
	paramNames.push_back("sink");
	paramNames.push_back("alpha");
	paramNames.push_back("beta");
	paramNames.push_back("TTexp");
	paramNames.push_back("TTinit");
}
Grain::Grain(std::string fileNameParam, std::string fileNameVar0,std::string outputdir):ParamGrain(fileNameParam, fileNameVar0)
{
}

Grain::~Grain()
{
}

////////////////////////////// Private Methods ////////////////////////////////

double Grain::ExportGrainPot(double tt, double dtt, double Npot)
{
	double r=0;
	if (tt <= (GetTT())){
		r = GetG()* Npot * dtt;
		//cout << r << endl;
	}else{
		r = GetG() * dtt * GetNGrainInit() * exp(GetG() * (GetTT())); 
		//cout << r << endl;
	}
	return r;
}

double Grain::NSynthRate(double ConcNmob,double dDM)
{
	return ConcNmob*GetP()*dDM; 
}


/*double Grain::ExportGrainPot_RK(double dt,double tt, double dtt)
{
	double dNPotdt1 = ExportGrainPot(tt,dtt,NPot);

	double NPot2 = NPot + dNPotdt1*dt/2;
	double dNPotdt2 = ExportGrainPot(tt+dtt/2,dtt,NPot2); 
	
	double NPot3 = NPot + dNPotdt2*dt/2;
	double dNPotdt3 = ExportGrainPot(tt+dtt/2,dtt,NPot3); 

	double NPot4 = NPot + dNPotdt3*dt;
	double dNPotdt4 = ExportGrainPot(tt+dtt,dtt,NPot4); 

	return (dNPotdt1/6 + dNPotdt2/3 + dNPotdt3/3 + dNPotdt4/6);
}*/



/////////////////////////////// Public Methods ////////////////////////////////
void Grain::InitialGrain()
{
	Ngrain = GetNGrainInit();
	Ngrains.push_back(Ngrain);
	
	NPot = GetNGrainInit();
	NPots.push_back(NPot);
	//this->parameters = parameters;

	TTinit=GetTTinit();

	double a=GetAlpha();
	double b=GetBeta();
	if(a>1 && b>1)
	{
		NormSink=pow((a-1)/(a+b-2),a-1)*pow((b-1)/(a+b-2),b-1);
	}
	else
	{
		NormSink=1.;
		cout << "Attention: coefficient fonction beta < 1"<<endl;
	}

	//DM=0.12; //voir si pertinent: Cf les grains sont d�j?pr�sents ?floraison: les faire appara�tre bien plus t�t au stade double rides ?
	
	//DM=0.108;
	DM = GetDMGrainInit(); 
	DMs.push_back(DM);
}

void Grain::ResetGrain()
{
	Ngrain=0;
	Ngrains.clear();
	NPot=0;
	NPots.clear();
	TTinit=0; // pourquoi la mettre en attribut ? Est-ce vraiment n�cessaire ?
	NormSink=0;
	DM=0;
	DMs.clear();
}

//ADD ELMER	
std::vector<double> Grain::GetNgrains() {
    return Ngrains;
}
std::vector<double> Grain::GetNPots() {
    return NPots;
}
std::vector<double> Grain::GetDMs() {
    return DMs;
}
//ADD ELMER



/*double Grain::ExportGrain_RK(double dt, double tt, double dtt, double Ncmob, double NmobMin)
{
	double dNPotdt1 = ExportGrainPot(tt,dtt,NPot);
	double dNdt1 = ExportGrain(tt,dtt,Ncmob,NmobMin,NPot);

	double NPot2 = NPot + dNPotdt1*dt/2;
	double dNPotdt2 = ExportGrainPot(tt+dtt/2,dtt,NPot2); 
	double dNdt2 = ExportGrain(tt+dtt/2,dtt,Ncmob,NmobMin,NPot2); 
	
	double NPot3 = NPot + dNPotdt2*dt/2;
	double dNPotdt3 = ExportGrainPot(tt+dtt/2,dtt,NPot3); 
	double dNdt3 = ExportGrain(tt+dtt/2,dtt,Ncmob,NmobMin,NPot3);

	double NPot4 = NPot + dNPotdt3*dt;
	double dNPotdt4 = ExportGrainPot(tt+dtt,dtt,NPot4); 
	double dNdt4 = ExportGrain(tt+dtt,dtt,Ncmob,NmobMin,NPot4); 

	return (dNdt1/6 + dNdt2/3 + dNdt3/3 + dNdt4/6);
}*/


double Grain::ExportGrain(double tt, double dtt, double Ncmob, double NmobMin, double Q, double D, double ConcNmob) // mettre le terme N quelque part pour ne pas confondre avec DM
{
	//if (tt <= GetTT())
	//{
		if (ExportGrainPot(tt,dtt,NPot)<(Ncmob-NmobMin))
		{
			return ExportGrainPot(tt,dtt,NPot);
		}else{
			if ((Ncmob-NmobMin)>0)
			{
				return (Ncmob-NmobMin);
			}else
			{
				return 0;
			}
		}
	//}else
	//{
	//	return NSynthRate(ConcNmob,ComputeDemand(tt)*Q/D);
	//}
}


void Grain::ActuGrain(double tt, double dtt, double dt, double Ncmob, double NmobMin, double Q, double D, double ConcNmob)// mettre le terme N quelque part pour ne pas confondre avec DM (et pas besoin de mettre Grain)
{
	//Ngrain += ExportGrain(tt,dtt,Ncmob,NmobMin)*dt;
	Ngrain += ExportGrain(tt,dtt,Ncmob,NmobMin,Q,D,ConcNmob);
	//Ngrain += NSynthRate(ConcNmob,ComputeDemand(tt)*Q/D);
	NPot +=  ExportGrainPot(tt,dtt,NPot);

	Ngrains.push_back(Ngrain);
	NPots.push_back(NPot);

	//cout<< endl;
	//cout << endl;
	//cout << tt << " " << NPot << endl;
}

double Grain::ComputeDemand(double TT) // mettre TT ? Mettre toutes ces m�thodes greenlab en commun pour Entity, Root, Grain ? 
{
	double u=(TT-TTinit)/GetTTExp();
	if(u<=0||u>=1)
	{
		return 0;
	}
	else
	{
		return GetSink()*pow(u,GetAlpha()-1)*pow(1-u,GetBeta()-1)/NormSink; // voir ?quoi correspond le param�tre Sink
	}
}

void Grain::ActuDryMass(double TT, double Q, double D)
{
	if (D==0) DM+=0;
	else DM+=ComputeDemand(TT)*Q/D;
	DMs.push_back(DM);
}

double Grain::GetDMgrain()
{
	return DM;
}

double Grain::GetNgrain()
{
	return Ngrain;
}

void Grain::Affich(std::vector<double> TTimes,std::string outputdir)
{
	string buffer;

	buffer=outputdir+"outputs/NgrainPred.txt";
	ofstream file1;
	file1.open(buffer.c_str());
	for (unsigned int i=0;i<Ngrains.size();i++)
	{
		//cout<<Ngrains[i]<<endl;
		file1 << TTimes[i] << "\t" << Ngrains[i] << "\t" << NPots[i] << "\n";
	}
	file1.close();

	buffer=outputdir+"outputs/MSgrainPred.txt";
	ofstream file2;
	file2.open(buffer.c_str());
	for (unsigned int i=0;i<DMs.size();i++)
	{
		file2 << TTimes[i]<< "\t" << DMs[i] << "\n";
	}
	file2.close();
}


// END OF File : Grain.cpp

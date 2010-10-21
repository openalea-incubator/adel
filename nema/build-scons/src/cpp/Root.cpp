/**
 * @file Root.cpp
 * @brief Implementation of class Root
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "Root.h"
#include "math.h"

// System includes ////////////////////////////////////////// System include //
#include <iostream>
#include <fstream>
#include <string>

// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
Root::Root()
{
	//Vector of parameter values
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);
	params.push_back(0);

	//Vector of parameter names
	paramNames.push_back("Umax");
	paramNames.push_back("kroot1");
	paramNames.push_back("kroot2");
	paramNames.push_back("QonDmin");
	paramNames.push_back("beta1");
	paramNames.push_back("beta2");
	paramNames.push_back("Nsoilmin");
	paramNames.push_back("p");
	paramNames.push_back("sink");
	paramNames.push_back("alpha");
	paramNames.push_back("beta");
	paramNames.push_back("TTexp");
	paramNames.push_back("DegDMrootRate");
	paramNames.push_back("DegrootRate");
	paramNames.push_back("TTinit");

}
Root::Root(std::string fileName,std::string outputdir):ParamRoot(fileName)
//Root::Root(std::string fileName, std::string fileNameVar0, std::string outputdir):ParamRoot(fileName, fileNameVar0) //ELMER


{
}

Root::~Root()
{
}

////////////////////////////// Private Methods ////////////////////////////////

double Root::NremSynthRate(double ConcNmob,double dDM)
{
	return ConcNmob*GetP()*dDM; 
}

double Root::NremDegRate(double dTT, double N)
{
	return dTT * GetDEGCOEF()*N;
}

double Root::DMremDegRate(double dTT, double dm)
{
	return dTT * GetDEGDMCOEF()*dm;
}

//ADD ELMER
std::vector<double> Root::GetNstructs() {
    return Nstructs;
}
std::vector<double> Root::GetNrems() {
    return Nrems;
}
std::vector<double> Root::GetDMstructs() {
    return DMstructs;
}
std::vector<double> Root::GetDMrems() {
    return DMrems;
}
//ADD ELMER


/////////////////////////////// Public Methods ////////////////////////////////

void Root::InitialRoot()
{
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

	
	DMrem = GetDMremInit();
	DMstruct = GetDMstructInit();
	
//	DMrem=0.0616;
//	DMstruct=0.5544;
	
	DMrems.push_back(DMrem);
	DMstructs.push_back(DMstruct);

//	Nstruct=0.009;
//	Nrem=0.00099;
	
	Nstruct = GetNstructInit();
	Nrem = GetNremInit();
	//Nrem=0;
	Nstructs.push_back(Nstruct);
	Nrems.push_back(Nrem);
}

void Root::ResetRoot()
{
  TTinit=GetTTinit(); // pourquoi la mettre en attribut ? Est-ce vraiment n�cessaire ?
	NormSink=0;
	DMrem=0;
	DMstruct=0;
	DMrems.clear();
	DMstructs.clear();
	Nstruct=0;
	Nrem=0;
	Nstructs.clear();
	Nrems.clear();

}
double Root::Neffect(double Nstatus)
{
	if (Nstatus <=0)
	{
		return 0;
	}else
	{
		return exp(-1*(GetBETA1())*(Nstatus));
	}
}

double Root::Ceffect(double Cstatus, int nbsubsteps)
{
	if (Cstatus<=0)
	{
		//cout << Cstatus << "\t" << 0 << endl;
		return 0;
	}
	else
	{
		//cout << Cstatus << "\t" << 1-exp(-1*(GetBETA2())*(Cstatus-(GetQONDMIN()))) << endl;
		return 1-exp(-1*(GetBETA2())*Cstatus);
	//return 1-exp(-1*(GetBETA2())*(Cstatus));
	}
}

double Root::ImportRoot(double dTT, int nbsubsteps, double Nstatus, double Cstatus, double Nsoil)// RK: cela doit être résolu à l'échelle plante entière
{
	if (Nsoil > GetNsoilMin())
	{
	//cout<< ((GetUMAX() * Nsoil)/((GetK1()) + Nsoil)) + (GetK2() * Nsoil) << "\t" << Neffect(Nstatus) << "\t" << Ceffect(Cstatus) << endl;
		if (Neffect(Nstatus)<=0 || Ceffect(Cstatus,nbsubsteps)<=0 || Cstatus<GetQONDMIN())
		{
			return 0;
		}else
		{// I put 1 to say that I need a minimum N concentration in the soil for the transpoters to be active
			return Neffect(Nstatus) * Ceffect(Cstatus,nbsubsteps)* dTT * (((GetUMAX() * (Nsoil-GetNsoilMin()) )/((GetKroot1()) + (Nsoil-GetNsoilMin())) + (GetKroot2() * (Nsoil-GetNsoilMin()) )));
		}
	}else
	{
		return 0;
	}
}

double Root::ImportRootPot(double dTT, double Nsoil)
{
	if (Nsoil > GetNsoilMin())
	{
		return dTT * (((GetUMAX() * (Nsoil-GetNsoilMin()) )/((GetKroot1()) + (Nsoil-GetNsoilMin())) + (GetKroot2() * (Nsoil-GetNsoilMin()) )));
	}else
	{
		return 0;
	}
}

double Root::Ninflux(double ConcNmob,double ConcNmobMin,double TT, double dTT, double Q, double D)
{
	if (NremSynthRate(ConcNmob-ConcNmobMin,ComputeDemand(TT)*Q/D)<=0 || D==0)
	{
		return (-1)*NremDegRate(dTT,Nrem);
	}else
	{
		return NremSynthRate(ConcNmob-ConcNmobMin,ComputeDemand(TT)*Q/D)-NremDegRate(dTT,Nrem);
	}
}

void Root::ActuNrem(double ConcNmob,double ConcNmobMin,double TT, double dTT, double Q, double D) // à fusionner avec ActuNph pour donner ActuN: attention: à bien actualiser avant Nph (Cf c'est Nph au pas de temps précédent qui compte)
{
	Nrem+=Ninflux(ConcNmob,ConcNmobMin,TT,dTT,Q,D);
	Nrems.push_back(Nrem);
	Nstructs.push_back(Nstruct);
}

double Root::ComputeDemand(double TT) // mettre TT ? Mettre toutes ces m�thodes greenlab en commun pour Entity, Root, Grain ? 
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

void Root::ActuDMrem(double TT, double dTT, double Q, double D)
{
	if (D==0) DMrem+=(-1)*DMremDegRate(dTT,DMrem);
	else DMrem+=ComputeDemand(TT)*Q/D-DMremDegRate(dTT,DMrem);
	DMrems.push_back(DMrem);
	DMstructs.push_back(DMstruct);
}

double Root::ComputeRemobilizedDM(double dTT, double TT)
{
	double remobilizedDM=DMremDegRate(dTT,DMrem);
	if (TT >= 800) remobilizedDM=0; // pas très catholique !!
	return remobilizedDM;
}

void Root::Affich(std::vector<double> TTimes,std::vector<double> dTTs, std::string outputdir)
{
	string buffer;

	buffer=outputdir+"outputs/NstructPredRoot.txt";
	ofstream file0;
	file0.open(buffer.c_str());
	for (unsigned int j=0;j<Nstructs.size();j++) 
	{
		file0 << TTimes[j] << "\t" << Nstructs[j] << "\n";
	}
	file0.close();

	buffer=outputdir+"outputs/NremPredRoot.txt";
	ofstream file1;
	file1.open(buffer.c_str());
	for (unsigned int j=0;j<Nrems.size();j++) 
	{
		file1 << TTimes[j] << "\t" << Nrems[j] << "\n";
	}
	file1.close();

	buffer=outputdir+"outputs/DMstructRootPred.txt";
	ofstream file2;
	file2.open(buffer.c_str());
	for (unsigned int i=0;i<DMstructs.size();i++)
	{
		file2 << TTimes[i] << "\t" << DMstructs[i] << "\n";
	}
	file2.close();

	buffer=outputdir+"outputs/DMremRootPred.txt";
	ofstream file3;
	file3.open(buffer.c_str());
	for (unsigned int i=0;i<DMrems.size();i++)
	{
		file3 << TTimes[i] << "\t" << DMrems[i] << "\n";
	}
	file3.close();

	buffer=outputdir+"outputs/NremDegRatePred.txt";
	ofstream file4;
	file4.open(buffer.c_str());
	for (unsigned int i=0;i<Nrems.size();i++)
	{
		if (i==0) file4 << TTimes[i] << "\t" << 0 << "\n";
		else file4 << TTimes[i] << "\t" << NremDegRate(dTTs[i],Nrems[i]) << "\n";
	}
	file4.close();

}



// END OF File : Root.cpp

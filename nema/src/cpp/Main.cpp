/**
 * @file Main.cpp
 * @brief Contain the main function
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

#ifdef WIN32
#include <cstdlib>
#endif
 
// Local includes /////////////////////////////////////////// Local Includes //
#include "Plant.h"
// CPL, ECR
//#include "sa.h"

// System includes ////////////////////////////////////////// System include //
#include <iostream>
#include <string>

// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

int main(int argc, char *argv[])
{
	//cout<<"Hello World !!"<<endl;
	string inputdir,outputdir,outf;
	//directory for parameters files
	string DParamEnv,DParamGrain,DParamRoot,DParamPlant;
	string DParamLamina,DParamSheath,DParamInternode,DParamPeduncle,DParamChaff;
	//directory for DrivingVariables
	string DNsoilH0,DTimeSoil,Dmeteo;
	//directory for State0dd
	
	
	//string DState0Grain,DState0Lamina,DState0Sheath,DState0Internode,DState0Peduncle,DState0Chaff; //ELMER
	string DState0Grain,DState0Lamina,DState0Sheath,DState0Internode,DState0Peduncle,DState0Chaff,DState0Root;
	
	
	//directory for sigmas files
	string DSigmasParamEnv,DSigmasParamGrain,DSigmasParamRoot,DSigmasParamPlant;
	string DSigmasParamLamina,DSigmasParamSheath,DSigmasParamInternode,DSigmasParamPeduncle,DSigmasParamChaff;
	//directory for simulation control flags
	//string DSim_base;
	
	if (argc==4)
	{
		inputdir=argv[1];
		outputdir=argv[2];
		outf=argv[3];
	}
	else
	{
		inputdir="";
		outputdir="";
		outf="";
	}

        /*	
        cout << "inputdir "<<inputdir<<endl;
        cout << "outputdir "<<outputdir<<endl;
        cout << "outf "<<outf<<endl;
        */

	//path of parameters files
	DParamPlant=inputdir+"parameters/ParamPlant.txt";
	DParamEnv=inputdir+"parameters/ParamEnv.txt";
	DParamGrain=inputdir+"parameters/ParamGrain.txt";
	DParamLamina=inputdir+"parameters/ParamLamina.txt";
	DParamSheath=inputdir+"parameters/ParamSheath.txt";
	DParamInternode=inputdir+"parameters/ParamInternode.txt";
	DParamPeduncle=inputdir+"parameters/ParamPeduncle.txt";
	DParamChaff=inputdir+"parameters/ParamChaff.txt";
	DParamRoot=inputdir+"parameters/ParamRoot.txt";

	//std::cout<<DParamPlant<<std::endl;

	//path of DrivingVariables 
	DNsoilH0=inputdir+"DrivingVariables/NsoilH2.txt"; // H0, H1, H2 refers to treatments H0, H3 and H15 (paper definition) 
	DTimeSoil=inputdir+"DrivingVariables/TimeSoil.txt";
	Dmeteo=inputdir+"DrivingVariables/meteo(old2).txt";

	//path of State0dd
	DState0Grain=inputdir+"State0dd/State0Grain.txt";
	DState0Lamina=inputdir+"State0dd/State0Lamina.txt";
	DState0Sheath=inputdir+"State0dd/State0Sheath.txt";
	DState0Internode=inputdir+"State0dd/State0Internode.txt";
	DState0Peduncle=inputdir+"State0dd/State0Peduncle.txt";
	DState0Chaff=inputdir+"State0dd/State0Chaff.txt";
	
	DState0Chaff=inputdir+"State0dd/State0Root.txt"; //ELMER
	
    //WQL For the path of Sigmas of every block of parameters
	DSigmasParamEnv=inputdir+"Sigmas/Sigmas_ParamEnv.txt";
	DSigmasParamGrain=inputdir+"Sigmas/Sigmas_ParamGrain.txt";
	DSigmasParamRoot=inputdir+"Sigmas/Sigmas_ParamRoot.txt";
	DSigmasParamPlant=inputdir+"Sigmas/Sigmas_ParamPlant.txt";
	DSigmasParamLamina=inputdir+"Sigmas/Sigmas_ParamLamina.txt";
	DSigmasParamSheath=inputdir+"Sigmas/Sigmas_ParamSheath.txt";
	DSigmasParamInternode=inputdir+"Sigmas/Sigmas_ParamInternode.txt";
	DSigmasParamPeduncle=inputdir+"Sigmas/Sigmas_ParamPeduncle.txt";
	DSigmasParamChaff=inputdir+"Sigmas/Sigmas_ParamChaff.txt";

	//WQL For path of simulation control flags
	//DSim_base=inputdir+"SA/simulations.txt";

	
	Plant maPlant(DParamPlant,
		DParamEnv,DNsoilH0,DTimeSoil,Dmeteo,
		DParamGrain,DState0Grain,
		DParamLamina,DState0Lamina,4,
		DParamSheath,DState0Sheath,4,
		DParamInternode,DState0Internode,4,
		DParamPeduncle,DState0Peduncle,1,
		DParamChaff,DState0Chaff,1,
		
		DParamRoot, 
		//DParamRoot,//ELMER
			
		outputdir);

        /*
	cout << DParamPlant << endl;
	cout << DParamEnv << "\t" << DNsoilH0<< "\t" <<DTimeSoil<< "\t" <<Dmeteo << endl;
	cout << DParamGrain << "\t" << DState0Grain << endl;
	cout << DParamLamina << "\t" << DState0Lamina << endl;
	cout << DParamSheath << "\t" << DState0Sheath << endl;
	cout <<	DParamInternode << "\t" << DState0Internode << endl;
	cout <<	DParamPeduncle << "\t" << DState0Peduncle << endl;
	cout <<	DParamChaff << "\t" << DState0Chaff << endl;
	cout <<	DParamRoot << endl;
	cout <<	outputdir << endl;
        */

	//all the sigmas
	ParamPlant SigmasParamPlant(DSigmasParamPlant,
		DSigmasParamEnv,DNsoilH0,DTimeSoil,Dmeteo,
		DSigmasParamGrain,DState0Grain,
		

		
		DSigmasParamLamina,DState0Lamina,4,
		DSigmasParamSheath,DState0Sheath,4,
		DSigmasParamInternode,DState0Internode,4,
		DSigmasParamPeduncle,DState0Peduncle,1,
		DSigmasParamChaff,DState0Chaff,1,
		
		//DSigmasParamRoot,DState0Root); //ELMER
		
		DSigmasParamRoot);
		
	
	/*ParamEnv SigmasParamEnv(DSigmasParamEnv,DNsoilH0,DTimeSoil,Dmeteo);
	ParamGrain SigmasParamGrain(DSigmasParamGrain,DState0Grain);
	ParamRoot SigmasParamRoot(DSigmasParamRoot);
	ParamEntity SigmasParamLamina(DSigmasParamLamina,DState0Lamina,4);
	ParamEntity SigmasParamSheath(DSigmasParamSheath,DState0Sheath,4);
	ParamEntity SigmasParamInternode(DSigmasParamInternode,DState0Internode,4);
	ParamEntity SigmasParamPeduncle(DSigmasParamPeduncle,DState0Peduncle,1);
	ParamEntity SigmasParamChaff(DSigmasParamChaff,DState0Chaff,10);*/
	
	/*paramLamina = new ParamEntity();
	paramSheath = new ParamEntity();
	paramInternode = new ParamEntity();
	paramPeduncle = new ParamEntity();
	paramChaff = new ParamEntity();*/
	//double GrowthDuration=800.; // en degr?jour ici
	
	// double GrowthDuration=40.; // en jour ici //  Variable not in use !
	int nSTEPS=49;// ?mettre en argv dans un premier temps
	int step=1;
	int nbsubsteps=4; //subdivision of the step nb in order to correctly solve differential equations with the euler method  
	maPlant.InitialPlant(nbsubsteps);//need to check parameter's value every time before model running
        
        cout << "Simulation of "<< nSTEPS << " steps "<<endl;
	for (int i=0;i<nSTEPS;i++)
	{
                //cout << "Simulation step "<<i<<endl;
		maPlant.ActuPlant(step,nbsubsteps);
	}
	maPlant.Affich(outputdir,nbsubsteps);

	//sa SA_Plant(DSim_base,outputdir,(Plant*)&maPlant,&SigmasParamPlant);

	//SA_Plant.SRC("DMgrain",SA_Plant.Getmaplant()-> DMgrains); // nécessité de mettre InitialPlant + écoulement du temps avant sinon le modèle ne tourne pas !!
	
	//SA_Plant.Sobol("RootAbsorb",SA_Plant.Getmaplant()->Imports);

        cout << "Simulation done."<<endl;
#ifdef WIN32
	system("pause");
#endif
	return 0;
}

// END OF File : Main.cpp

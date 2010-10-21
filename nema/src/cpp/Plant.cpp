/**
 * @file Plant.cpp
 * @brief Implementation of class Plant
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "Plant.h"

// System includes ////////////////////////////////////////// System include //
#include <iostream>
#include <fstream>
#include<string>
#include <cmath>
#include <cassert>

// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
Plant::Plant()
{
	grain = new Grain();
	root = new Root();
	lamina = new Entity();
	sheath = new Entity();
	internode = new Entity();
	peduncle = new Entity();
	chaff = new Entity();
	environment = new Environment();

	//vector of parameter values
	params.push_back(0);
	//Vector of parameter names
	paramNames.push_back("ConcNmobMin");
}
Plant::Plant(std::string fileNamePlant,
		std::string fileNameParamEnv,std::string fileNameNEnv, std::string fileNameTEnv, std::string fileNameMeteoEnv,
		std::string fileNameParamGrain, std::string fileNameVar0Grain,
		std::string fileNameParamLamina, std::string fileNameVar0Lamina, int nrowLamina,
		std::string fileNameParamSheath, std::string fileNameVar0Sheath, int nrowSheath,
		std::string fileNameParamInternode, std::string fileNameVar0Internode, int nrowInternode,
		std::string fileNameParamPeduncle, std::string fileNameVar0Peduncle, int nrowPeduncle,
		std::string fileNameParamChaff, std::string fileNameVar0Chaff, int nrowChaff,
		
		std::string fileNameRoot,
		//std::string fileNameRoot, std::string fileNameVar0Root,//ELMER
		
		std::string outputdir)
		:ParamPlant(fileNamePlant,
		 fileNameParamEnv, fileNameNEnv,  fileNameTEnv,  fileNameMeteoEnv,
		 fileNameParamGrain,  fileNameVar0Grain,
		 fileNameParamLamina,  fileNameVar0Lamina, nrowLamina,
		 fileNameParamSheath,  fileNameVar0Sheath, nrowSheath,
		 fileNameParamInternode,  fileNameVar0Internode, nrowInternode,
		 fileNameParamPeduncle,  fileNameVar0Peduncle, nrowPeduncle,
		 fileNameParamChaff,  fileNameVar0Chaff, nrowChaff,
		 
		 fileNameRoot)
		 //fileNameRoot, fileNameVar0Root)//ELMER
{
	//parameters = paramPlant;
	grain = new Grain(fileNameParamGrain, fileNameVar0Grain,outputdir);
	
	
	root = new Root(fileNameRoot,outputdir);
	//root = new Root(fileNameRoot,fileNameVar0Root,outputdir); //ELMER
	
	lamina = new Entity(fileNameParamLamina, fileNameVar0Lamina, nrowLamina, outputdir);
	sheath = new Entity(fileNameParamSheath, fileNameVar0Sheath, nrowSheath, outputdir);
	internode = new Entity(fileNameParamInternode, fileNameVar0Internode, nrowInternode, outputdir);
	peduncle = new Entity(fileNameParamPeduncle, fileNameVar0Peduncle, nrowPeduncle, outputdir);
	chaff = new Entity(fileNameParamChaff, fileNameVar0Chaff, nrowChaff, outputdir);
	environment = new Environment(fileNameParamEnv,fileNameNEnv, fileNameTEnv, fileNameMeteoEnv);
}

//virere environment in this constructor ??

Plant::Plant( const Grain& _grain, 
              const Root& _root, 
              const Environment& _environment, 
              const Entity& _lamina, 
              const Entity& _sheath, 
              const Entity& _internode, 
              const Entity& _peduncle, 
              const Entity& _chaff)
{
    grain = new Grain( _grain );
    root = new Root( _root );
    environment = new Environment( _environment );
    lamina = new Entity( _lamina );
    sheath = new Entity(_sheath);
    internode = new Entity( _internode );
    peduncle = new Entity(_peduncle);
    chaff = new Entity( _chaff );
}

Plant::~Plant()
{
	delete grain;
	delete root;
	delete lamina;
	delete sheath;
	delete internode;
	delete peduncle;
	delete chaff;
	delete environment;
}

////////////////////////////// Public Methods ////////////////////////////////



void Plant::InitialPlant(int nbsubsteps) // Attributes initialization for each class
{
	grain->InitialGrain();
	root->InitialRoot();
	//environment->InitialEnv();

	//cout << "lamina" << endl;
	lamina->InitialEntity();
	//cout << "sheath" << endl;
	sheath->InitialEntity();
	//cout << "internode" << endl;
	internode->InitialEntity();
	//cout << "peduncle" << endl;
	peduncle->InitialEntity();
	//cout << "chaff" << endl;
	chaff->InitialEntity();

	Time = 1;
	TTime = 1.;
	TTimes.push_back(TTime);
	dTT=0.;
	dTTs.push_back(dTT);
	plt_density = 410; // mettre cela dans arg de main

	types.push_back(LAMINA);
	types.push_back(SHEATH);
	types.push_back(INTERNODE);
	types.push_back(PEDUNCLE);
	types.push_back(CHAFF);
	
	//Nmob=0.0023;
	Nmob=0.0024; 
	Nmobs.push_back(Nmob);
	DMgrains.push_back(grain->GetDMgrain());// pour SA: à virer ensuite
	Ngrains.push_back(grain->GetNgrain());// pour SA: à virer ensuite

	ComputeOrganHeightMaxFromSoil(types); 
	
	ComputeVisibleAreas(types);
	ActuRowExposed(types);
	
	ComputeAreaGreenTot(types);
	ComputeDMGreenTot(types);

	//std::vector<std::vector<double> > PARss;
	for (int unsigned ty=0;ty<types.size();ty++)
	{
		//cout << types[ty] << "\t";
		//calcul de la matrice de GAI correspondant au GAI cumulé pour les différents types d'organes 
		std::vector<std::vector<double> > gai_cumss;
		for (unsigned int ty2=0;ty2<types.size();ty2++)
		{
			gai_cumss.push_back(GAI_cums(types[ty],types[ty2]));
			
		}
		gai_cumsss.push_back(gai_cumss); // dimension 1= type de l'organe considéré; dimension 2=type2(=type pour lequel on cumul le GAI);dimension 3= rang de l'organe considéré
		PARss.push_back(environment->PARs(Time,nbsubsteps,types[ty],nrow(types[ty]),rowss_exposed[types[ty]-1],gai_cumss,lamina->GetInsertionAngle(),HeightViss(types[ty])));
		/*for (int r=0;r<nrow(types[ty]);r++)
		{
			cout << PARss[types[ty]-1][r] << "\t";
		}
		//cout << "-------------Nouveau type d'organe--------------------------" << endl;
		cout << endl;*/
	}
	PARsss.push_back(PARss);

	Imports.push_back(0);
	ImportPots.push_back(0);
	Nsoil=environment->GetNsoil(Time);
	//cout << Nsoil << endl;
	ComputeProduction(types,PARss,Time,nbsubsteps,TTime);
	ComputeRemobilizedDM(types,dTT,TTime);
	ComputeDemand(types,TTime); 

	/*std::vector<std::vector<double> > gai_cumss;
	for (unsigned int ty2=0;ty2<types.size();ty2++)
	{
		gai_cumss.push_back(GAI_cums(SHEATH,types[ty2]));
		cout << types[ty2] << "\t";
		for (int r=0;r<nrow(types[ty2]);r++)
		{
			cout << gai_cumss[ty2][r] << "\t";
		}
		cout << endl;
	}*/


	//GAI_cums(LAMINA,SHEATH);
	/*for (unsigned int r=0;r<GAI_cums(CHAFF,LAMINA).size();r++)
	{
		cout << GAI_cums(CHAFF,LAMINA)[r] << endl;
	}*/
	
}

//version that don't need environment

void Plant::InitialPlantExt(int nbsubsteps,double _Nmob,std::vector<std::vector<double> > _PARss,double _Nsoil) // Attributes initialization for each class
{
	grain->InitialGrain();
	root->InitialRoot();
	//environment->InitialEnv();

	//cout << "lamina" << endl;
	lamina->InitialEntity();
	//cout << "sheath" << endl;
	sheath->InitialEntity();
	//cout << "internode" << endl;
	internode->InitialEntity();
	//cout << "peduncle" << endl;
	peduncle->InitialEntity();
	//cout << "chaff" << endl;
	chaff->InitialEntity();

	Time = 1;
	TTime = 1.;
	TTimes.push_back(TTime);
	dTT=0.;
	dTTs.push_back(dTT);
	plt_density = 410; // mettre cela dans arg de main

	types.push_back(LAMINA);
	types.push_back(SHEATH);
	types.push_back(INTERNODE);
	types.push_back(PEDUNCLE);
	types.push_back(CHAFF);
	
	
	Nmob = _Nmob; 
	Nmobs.push_back(Nmob);
	DMgrains.push_back(grain->GetDMgrain());// pour SA: à virer ensuite
	Ngrains.push_back(grain->GetNgrain());// pour SA: à virer ensuite

	ComputeOrganHeightMaxFromSoil(types); 
	
	ComputeVisibleAreas(types);
	ActuRowExposed(types);
	
	ComputeAreaGreenTot(types);
	ComputeDMGreenTot(types);

	//std::vector<std::vector<double> > PARss;
	for (int unsigned ty=0;ty<types.size();ty++)
	{
		//cout << types[ty] << "\t";
		//calcul de la matrice de GAI correspondant au GAI cumulé pour les différents types d'organes 
		std::vector<std::vector<double> > gai_cumss;
		for (unsigned int ty2=0;ty2<types.size();ty2++)
		{
			gai_cumss.push_back(GAI_cums(types[ty],types[ty2]));
			
		}
		gai_cumsss.push_back(gai_cumss); // dimension 1= type de l'organe considéré; dimension 2=type2(=type pour lequel on cumul le GAI);dimension 3= rang de l'organe considéré
		//PARss.push_back(environment->PARs(Time,nbsubsteps,types[ty],nrow(types[ty]),rowss_exposed[types[ty]-1],gai_cumss,lamina->GetInsertionAngle(),HeightViss(types[ty])));
		/*for (int r=0;r<nrow(types[ty]);r++)
		{
			cout << PARss[types[ty]-1][r] << "\t";
		}
		//cout << "-------------Nouveau type d'organe--------------------------" << endl;
		cout << endl;*/
	}
	assert( _PARss.size() == types.size()); 
	PARss = _PARss;
	PARsss.push_back(_PARss);

	Imports.push_back(0);
	ImportPots.push_back(0);
	//Nsoil=environment->GetNsoil(Time);
	Nsoil = _Nsoil;
	//cout << Nsoil << endl;
	ComputeProduction(types,PARss,Time,nbsubsteps,TTime);
	ComputeRemobilizedDM(types,dTT,TTime);
	ComputeDemand(types,TTime); 

	/*std::vector<std::vector<double> > gai_cumss;
	for (unsigned int ty2=0;ty2<types.size();ty2++)
	{
		gai_cumss.push_back(GAI_cums(SHEATH,types[ty2]));
		cout << types[ty2] << "\t";
		for (int r=0;r<nrow(types[ty2]);r++)
		{
			cout << gai_cumss[ty2][r] << "\t";
		}
		cout << endl;
	}*/


	//GAI_cums(LAMINA,SHEATH);
	/*for (unsigned int r=0;r<GAI_cums(CHAFF,LAMINA).size();r++)
	{
		cout << GAI_cums(CHAFF,LAMINA)[r] << endl;
	}*/
	
}


void Plant::ResetPlant()
{
	grain->ResetGrain();
	root->ResetRoot();
	//environment->ResetEnv();
	lamina->ResetEntity();
	sheath->ResetEntity();
	internode->ResetEntity();
	peduncle->ResetEntity();
	chaff->ResetEntity();

	Time=1;
	TTime=1.;
	TTimes.clear();
	dTTs.clear();
	dTT=0.;
	plt_density=410;
	types.clear();
	Nsoil=0;
	Nmob=0;
	Nmobs.clear();
	DMgrains.clear();
	Ngrains.clear();
	DMGreenTot=0;
	DMGreenTots.clear();
	Imports.clear();
	ImportPots.clear();
	gai_cumsss.clear();
	PARss.clear();
	PARsss.clear();
	AreaGreenTots.clear();
	HeightMaxss.clear();
	visibleAreass.clear();
	rowss_exposed.clear();

	Production=0; 
	Productions.clear();
	Demand=0;
	Demands.clear();

}

int Plant::nrow(int type)
{
  Entity* e=  TypeToClass(type);
  if( e != 0)
    return e->GetNbRow();
  else
    return 0;
}


Entity* Plant::TypeToClass(int type)
{
  switch(type)
  {
     case LAMINA:
		 return lamina;
		 break;
	 case SHEATH:
		 return sheath;
		 break;
	 case INTERNODE:
		 return internode;
		 break;
	 case PEDUNCLE:
		 return peduncle;
		 break;
	 case CHAFF:
		 return chaff;
		 break;
	 default:
	     return NULL;
  }

}

void Plant::ComputeDMGreenTot(std::vector<int> typs)
{
	DMGreenTot=0;
	for (int unsigned ty=0;ty<typs.size();ty++)
	{
		for (int r=0;r<nrow(typs[ty]);r++)
		{
			DMGreenTot+=(TypeToClass(typs[ty]))->DMGreen(r+1);
		}
	}
	DMGreenTots.push_back(DMGreenTot);
}


void Plant::ComputeNmob(std::vector<int> typs,double ConcNmob,std::vector<std::vector<double> > PARss, double Absorb, double NInfluxGrain, double NInfluxRoot)
{
	double NphInfluxs=0;
	for (int unsigned ty=0;ty<typs.size();ty++)
	{
		for (int r=0;r<nrow(typs[ty]);r++)
		{
			NphInfluxs+=(TypeToClass(typs[ty]))->NphInput(r+1,TTime,dTT,ConcNmob, GetConcNmobMin(),PARss[typs[ty]-1][r],Production,Demand)*(-1);
		}
	}
	//cout << NphInfluxs << "\t" << Absorb << "\t" << NInfluxGrain <<"\t" << NInfluxRoot <<"\t" << endl;
	Nmob+=NphInfluxs+Absorb-NInfluxGrain-NInfluxRoot; // voir si je mets pas directement la valeur de absorb et Ninfluxgrain
	//if (Nmob<=0) Nmob=0;
	Nmobs.push_back(Nmob);
}

double Plant::ConcNmob(double Nmob,double DMgreen)
{
	if (DMgreen==0) return 0;
	else return Nmob/DMgreen;
}

void Plant::ComputeProduction(std::vector<int> typs, std::vector<std::vector<double> > PARss, int time, int nbsubsteps, double ttime)
{
	Production=0;
	//cout << "PAR lamina rang 1" << PARss[0][0] << endl;
	for (int unsigned ty=0;ty<typs.size();ty++)
	{
		Production+=(TypeToClass(typs[ty]))->ComputeProduction(nrow(typs[ty]),PARss[typs[ty]-1],time,nbsubsteps,ttime);// voir si je dois mettre ttime
	}
	//if (ttime>=800) Production=0; 
	Productions.push_back(Production);
}

void Plant::ComputeRemobilizedDM(std::vector<int> typs, double dTT, double TT)
{
	remobilizedDM=0;
	for (int unsigned ty=0;ty<typs.size();ty++)
	{
		remobilizedDM+=(TypeToClass(typs[ty]))->ComputeRemobilizedDM(nrow(typs[ty]),dTT,TT);
	}
	remobilizedDM+=root->ComputeRemobilizedDM(dTT,TT);
	remobilizedDMs.push_back(remobilizedDM);
}

void Plant::ComputeDemand(std::vector<int> typs,double ttime)
{
	Demand=0;
	for (int unsigned ty=0;ty<typs.size();ty++)
	{
		Demand+=(TypeToClass(typs[ty]))->ComputeDemand(nrow(typs[ty]),ttime);// voir pourquoi besoin de mettre TTime
		//cout << Demand << endl;
	}
	Demand+=grain->ComputeDemand(ttime);
	Demand+=root->ComputeDemand(ttime);
	Demands.push_back(Demand);
}


void Plant::ComputeOrganHeightMaxFromSoil(std::vector<int> typs) // height of the upper part of the organ from the soil
{
	//cout << "type" << "\t" << "row" << "\t" << "length" << "\t" << endl;
	for (int unsigned ty=0;ty<typs.size();ty++)
	{
		std::vector<double> vect;
		HeightMaxss.push_back(vect);
		for (int r=1;r<(nrow(typs[ty])+1);r++)
		{
			//cout << typs[ty] << "\t" << r << "\t";
			if (typs[ty]==LAMINA)
			{
				double height_internodes=0.;
				for (int i=1;i<(r+1);i++) // sum of the length of the internodes below lamina r (comprising internode r)
				{
					//cout << i << "-";
					height_internodes += internode->GetLengths0()[i-1];
				}
				HeightMaxss[ty].push_back( height_internodes + sheath->GetLengths0()[r-1] + cos(1.57-lamina->GetInsertionAngle()) * lamina->GetLengths0()[r-1] );
				//cout << HeightMaxss[ty][r-1] << "\t" << height_internodes << "\t" << sheath->GetLengths0()[r-1] << "\t" << cos(lamina->GetInsertionAngle()) * lamina->GetLengths0()[r-1]  << endl;
			}else
			{
				if (typs[ty]==SHEATH)
				{
					double height_internodes=0.;
					for (int i=1;i<(r+1);i++)
					{
						//cout << i << "-";
						height_internodes += internode->GetLengths0()[i-1];
					}
					HeightMaxss[ty].push_back( height_internodes + sheath->GetLengths0()[r-1] );
					//cout << sheath->GetLengths0()[r-1] << endl;
				}else
				{
					if (typs[ty]==INTERNODE)
					{
						double height_internodes=0.;
						for (int i=1;i<(r+1);i++)
						{
							//cout << i << "-";
							height_internodes += internode->GetLengths0()[i-1];
						}
						HeightMaxss[ty].push_back(height_internodes);
						//cout << internode->GetLengths0()[r-1] << endl;
					}else
					{
						if (typs[ty]==PEDUNCLE)
						{
							double height_internodes=0.;
							for (int i=1;i<(internode->GetNbRow()+1);i++)
							{
								//cout << i << "-";
								height_internodes += internode->GetLengths0()[i-1];
							}
							HeightMaxss[ty].push_back( height_internodes + peduncle->GetLengths0()[0] );
							//cout << peduncle->GetLengths0()[r-1] << endl;
						}else
						{
							if (typs[ty]==CHAFF)
							{
								double height_internodes=0.;
								for (int i=1;i<(internode->GetNbRow()+1);i++)
								{
									//cout << i << "-";
									height_internodes += internode->GetLengths0()[i-1];
								}
								HeightMaxss[ty].push_back( height_internodes + peduncle->GetLengths0()[0] + chaff->GetLengths0()[0] );
								//cout << chaff->GetLengths0()[r-1] << endl;
							}
						}
					}
				}
			}
			//cout << "\t" << HeightMaxss[ty][r-1] << endl;
		}
	}

	/*cout << "type" << "\t" << "row" << "\t" << "hauteur" << endl;
	for (unsigned int ty=0;ty<typs.size();ty++)
	{
		for (int r=0;r<nrow(typs[ty]);r++)
		{
			cout << typs[ty] << "\t" << r << "\t" << HeightMaxss[ty][r] << endl;
		}
	}*/
}

double Plant::HeightPreviousOrgan(int type, int row)
{
	if (type==LAMINA)
	{
		return HeightMaxss[SHEATH-1][row-1];
	}else
	{
		if(type==SHEATH || type==INTERNODE) // attention: gaine 1 et EN1: mettre des valeurs par défault vu que les autres organes ne sont pas là mais qu'ils les cahe dans la réalité
		{
			if (row>1)
			{
				return HeightMaxss[SHEATH-1][row-2];
			}else
			{
				if (type==SHEATH) // I suppose that half of the first sheath is exposed: this has to be improve later
				{
					return HeightMaxss[SHEATH-1][row-1]/2;
				}else
				{
					return HeightMaxss[INTERNODE-1][row-1]+1; // I suppose that the first internode is not exposed
				}
			}
		}else
		{
			if (type==PEDUNCLE)
			{
			  return HeightMaxss[SHEATH-1][nrow(SHEATH)-1];
			}else
			{
				if (type==CHAFF)
				{
					return HeightMaxss[PEDUNCLE-1][0];
				}
			}
		}
	}
}

std::vector<double> Plant::HeightViss(int type) // pour prendre en compte la hauteur de l'organe dans le calcul du rayonnement intercepté
{ // en fait, cela ne fonctionne pas: donc, cette fonction n'est pas utilisée en réalité
	std::vector<double> HeightViss;
	for (int r=0;r<nrow(type);r++)
	{
		double HeightVis;
		if (type==LAMINA)
		{
			HeightVis=TypeToClass(type)->VisibleLength(type,r+1,HeightMaxss[type-1][r],HeightPreviousOrgan(type,r+1)*cos(1.57-lamina->GetInsertionAngle()));
		}else
		{
			HeightVis=TypeToClass(type)->VisibleLength(type,r+1,HeightMaxss[type-1][r],HeightPreviousOrgan(type,r+1));
		}
		HeightViss.push_back(HeightVis);
	}
	return HeightViss;
}

void Plant::ComputeVisibleAreas(std::vector<int> typs) // valable si uniquement utilisé dans l'initialisation (Cf dans .push_back)
{ 
	//cout << endl;
	//cout << "Stot" << "\t" << "Lvis" << "\t" << "Ltot" << "\t" << "type" << "\t" << "row" << "\t" << "Svis" << endl;
	for (unsigned int ty=0;ty<typs.size();ty++)
	{
	//int ty=1;
		std::vector<double> visibleAreas;
		for (int r=0;r<nrow(typs[ty]);r++)
		{
			//cout << HeightMaxss[ty][r] << "\t" << endl;
			//cout << (TypeToClass(typs[ty]))->VisibleLength(HeightMaxss[ty][r],HeightPreviousOrgan(typs[ty],r+1)) << endl;
			double visibleArea=(TypeToClass(typs[ty]))->VisibleArea(r+1,(TypeToClass(typs[ty]))->VisibleLength(typs[ty],r+1,HeightMaxss[ty][r],HeightPreviousOrgan(typs[ty],r+1)));
			if ( visibleArea>0 )
			{
				visibleAreas.push_back(visibleArea);
			}else
			{
				visibleAreas.push_back(0);
			}
			//cout << typs[ty] << "\t" << r+1 << "\t" << visibleAreas[r] << endl;
		}
		//cout << typs[ty] << "\t";
		TypeToClass(typs[ty])->InitializeAreaGreensNm2Max(nrow(typs[ty]),visibleAreas);
		visibleAreass.push_back(visibleAreas);
	}
}


void Plant::ActuRowExposed(std::vector<int> typs)
{
	for (unsigned int ty=0;ty<typs.size();ty++)
	{
		std::vector<int> vect;
		for (int r=0;r<nrow(typs[ty]);r++)
		{
			if ( (TypeToClass(typs[ty])->VisibleLength(typs[ty],r+1,HeightMaxss[typs[ty]-1][r],HeightPreviousOrgan(typs[ty],r+1))) >0)
			{
				vect.push_back(r+1);
			}
		}
		rowss_exposed.push_back(vect);
		/*cout << typs[ty] << "\t";
		for (int r=0;r<rowss_exposed[ty].size();r++)
		{
			cout << rowss_exposed[ty][r] << "\t";
		}
		cout << endl;*/
	}
}

double Plant::MeanHeightExposedPart(int type, int row)
{
	double visible_length=TypeToClass(type)->VisibleLength(type,row,HeightMaxss[type-1][row-1],HeightPreviousOrgan(type,row));
	double InsertionAngle_Factor;
	if (type==LAMINA)
	{
		InsertionAngle_Factor=cos(1.57-lamina->GetInsertionAngle());
	}else
	{
		InsertionAngle_Factor=1;
	}
	if (visible_length!=0.)
	{
		return HeightMaxss[type-1][row-1]-(visible_length/2*InsertionAngle_Factor);
	}else
	{
		cout << "error in the program: the entity considered is not exposed" << endl;// Cf i do not compute the mean height of the organ if it is not exposed 
	}
}



std::vector<double> Plant::GAI_cums(int type, int type2)// this is only computed for organs that are visible (see the PAR calculation fonction)
{ // the first type of the function corresponds to the type of organ for which cumulated gai above is computed
	// the second type corresponds to the type of organ for which we sum the GAIs (this is necessary for the considereration of different extinction coefficients in the PAR calculation)
	std::vector<double> gai_cums; 
	//cout << "nombres de rangs exposes=" << rowss_exposed[type-1].size() <<"\t";
	if (rowss_exposed[type-1].size()>0)
	{
		// sum of the areas that are above the organ
		for (unsigned int r=0;r<rowss_exposed[type-1].size();r++)// I have to make the loop for gai calculation only on the exposed organs since the calculation of the mean height is only for exposed parts of organs
		{
			//cout << "-------------rang expose" << r+1 << " -------------------" << endl;
			//I scan all the organs (types and rows) to check which one is above
			double area_above=0; // for the organ considered
			//for (unsigned int ty=0;ty<types.size();ty++)
			//{
			int ty2=type2-1;
			//cout << "type2=" << ty2 << "\t" << "nrow_ty2= " << nrow(types[ty2]) << endl;

			for (int ri=0;ri<nrow(types[ty2]);ri++)
			{
				//cout << "HtMax type2=" << HeightMaxss[types[ty2]-1][ri] << "\t" << "HtMoy type=" << MeanHeightExposedPart(type,rowss_exposed[type-1][r]) << "\t" << "HtMax PreviousType2=" << HeightPreviousOrgan(types[ty2],ri+1) << endl;
				if ( HeightMaxss[types[ty2]-1][ri] >= MeanHeightExposedPart(type,rowss_exposed[type-1][r]) )
				{
					if ( HeightPreviousOrgan(types[ty2],ri+1) >= MeanHeightExposedPart(type,rowss_exposed[type-1][r]) ) // if the base of the organ is above the average height of the organ considered (the one for which we compute GAI), we add the entire visible area of the organ
					{
						//cout << "l'organe entier est au dessus" << "\t" << "Lvisible=" << (TypeToClass(types[ty2]))->VisibleLength(types[ty2],ri+1,HeightMaxss[types[ty2]-1][ri],HeightPreviousOrgan(types[ty2],ri+1)) << "\t" << "Surface visible=" << (TypeToClass(types[ty2]))->VisibleArea(ri+1,(TypeToClass(types[ty2]))->VisibleLength(types[ty2],ri+1,HeightMaxss[types[ty2]-1][ri],HeightPreviousOrgan(types[ty2],ri+1))) << endl;
						area_above += (TypeToClass(types[ty2]))->VisibleArea(ri+1,(TypeToClass(types[ty2]))->VisibleLength(types[ty2],ri+1,HeightMaxss[types[ty2]-1][ri],HeightPreviousOrgan(types[ty2],ri+1)));
					}else // else we add the fraction of area above
					{
						double visible_area = (TypeToClass(types[ty2]))->VisibleArea(ri+1,(TypeToClass(types[ty2]))->VisibleLength(types[ty2],ri+1,HeightMaxss[types[ty2]-1][ri],HeightPreviousOrgan(types[ty2],ri+1)));
						double visible_height = HeightMaxss[types[ty2]-1][ri]-HeightPreviousOrgan(types[ty2],ri+1);
						double height_above = HeightMaxss[types[ty2]-1][ri]-MeanHeightExposedPart(type,rowss_exposed[type-1][r]);
						area_above += visible_area*height_above/visible_height;
						//cout << "surface visible de l'organe=" << visible_area << endl;
						//cout << "une fraction de l'organe est au dessus" << "\t" << "Ht au-dessus=" << height_above << "\t" << "ht visible=" << visible_height << "\t" << "S au-dessus=" << visible_area*height_above/visible_height << endl;
					}
				}
			}
			//}
			/*if (type2==type)
			{
				area_above += ( (TypeToClass(type))->VisibleArea(rowss_exposed[type-1][r],(TypeToClass(type))->VisibleLength(type,rowss_exposed[type-1][r],HeightMaxss[type-1][rowss_exposed[type-1][r]-1],HeightPreviousOrgan(type,rowss_exposed[type-1][r]))) )/2; // half of the visible area of the organ considered
			}*/

			gai_cums.push_back(area_above*plt_density);
			//cout << endl;
			//cout << "Stotale au dessus=" << area_above << endl;
		}
	}
	//cout << endl;
	return gai_cums;

}


/*std::vector<double> Plant::GAI_cums(int type)// mettre en fonction des longueurs et de l'angle des feuilles
{	std::vector<double> area_cums;
	std::vector<double> gai_cums;
	for (int row=1;row<(nrow(type)+1);row++)
	{
		area_cums.push_back(0);
		if (type==LAMINA)
		{
			if (row<nrow(type)) 
			{
				for (int r=row+1;r<(nrow(type)+1);r++)
				{
					area_cums[row-1]+=lamina->GetAreaGreen0(r)+sheath->GetAreaGreen0(r)+internode->GetAreaGreen0(r);
				}
				area_cums[row-1]+= lamina->GetAreaGreen0(row)/2 + peduncle->GetAreaGreen0(1) + chaff->GetAreaGreen0(1);
			}else
			{
				area_cums[row-1]+=peduncle->GetAreaGreen0(1)+chaff->GetAreaGreen0(1);
			}
		}else
		{
			if (type==SHEATH)// ici, je consid�re que le limbe qui est parall�le ?la gaine ne lui fait pas d'ombre: en v�rifier la pertinence
			{
				area_cums[row-1]+=sheath->GetAreaGreen0(row)/2;
				if (row==nrow(type))
				{
					area_cums[row-1]+=lamina->GetAreaGreen0(row) + peduncle->GetAreaGreen0(1) + chaff->GetAreaGreen0(1);
				}
				if (row<nrow(type)) 
				{
					area_cums[row-1]+=lamina->GetAreaGreen0(row);
					for	(int r=row+1;r<(nrow(type)+1);r++)
					{
						area_cums[row-1]+=lamina->GetAreaGreen0(r);
						area_cums[row-1]+=sheath->GetAreaGreen0(r);
						area_cums[row-1]+=internode->GetAreaGreen0(r);
					}
					area_cums[row-1]+=peduncle->GetAreaGreen0(1)+chaff->GetAreaGreen0(1);
				}
			}else
			{
				if (type==PEDUNCLE)
				{
					area_cums[row-1]+=peduncle->GetAreaGreen0(1)/2;
					area_cums[row-1]+=chaff->GetAreaGreen0(1);
				}else
				{
					if (type==CHAFF)
					{
						area_cums[row-1]+=chaff->GetAreaGreen0(1)/2;
					}else
					{
						area_cums[row-1]+=internode->GetAreaGreen0(1)/2;
						if (row==nrow(type))
						{
							area_cums[row-1]+=lamina->GetAreaGreen0(row)+sheath->GetAreaGreen0(row)+peduncle->GetAreaGreen0(1)+chaff->GetAreaGreen0(1);
						}
						if (row<nrow(type)) 
						{
							area_cums[row-1]+=lamina->GetAreaGreen0(row)+sheath->GetAreaGreen0(row);
							for	(int r=row+1;r<(nrow(type)+1);r++)
							{
								area_cums[row-1]+=lamina->GetAreaGreen0(r);
								area_cums[row-1]+=sheath->GetAreaGreen0(r);
								area_cums[row-1]+=internode->GetAreaGreen0(r);
							}
							area_cums[row-1]+=peduncle->GetAreaGreen0(1)+chaff->GetAreaGreen0(1);
						}
					}
				}
			}
		}
		gai_cums.push_back(area_cums[row-1]*plt_density);
		//printf("%g\n",area_cums[row-1]);
	}
	//cout<<Area0s[row-1]<<endl;
	return gai_cums;
}*/

void Plant::ComputeAreaGreenTot(std::vector<int> typs)
{
	double AreaGreenTot=0;
	for (int unsigned ty=0;ty<typs.size();ty++)
	{
		AreaGreenTot+=(TypeToClass(typs[ty]))->GetAreaGreenTot(nrow(typs[ty]));
	}
	AreaGreenTots.push_back(AreaGreenTot);
}

std::vector<double> Plant::CumulGAI(int typ, std::vector<std::vector<double> > gaiss)
{
	std::vector<double> gais;
	//cout << rowss_exposed[typ-1].size() << endl;
	for (int r=0;r<nrow(typ);r++)
	{
		int OrganPositionInGAIs=-1;
		for (unsigned int ri=0;ri<rowss_exposed[typ-1].size();ri++)
		{
			if ( (r+1)==rowss_exposed[typ-1][ri] )
			{
				OrganPositionInGAIs=ri;
			}
		}
		if (OrganPositionInGAIs>(-1))
		{
			double gai=0;
			for (int unsigned ty2=0;ty2<gaiss.size();ty2++)
			{
				gai+=gaiss[ty2][OrganPositionInGAIs];
			}
			//cout << typ << "\t" << r <<  "\t" << gai << endl;  
			gais.push_back(gai);
		}
	}
	return gais;
}

/////////////////////////////// Public Methods ////////////////////////////////
/*void Plant::ComputeFluxes(int step) // voir si cela sert de mettre step en argument
{
	Time=step;
	double FluxGrain=grain->ExportGrain(Time,Nmob,GetConcNmobMin()*DMGreenTot);

}*/


void Plant::ActuPlant(int step, int nbsubsteps) 
{
	Time+=step;
	for (int j=0;j<nbsubsteps;j++)
	{
		dTT=environment->ThermalTime(Time,0.,step,nbsubsteps);
		dTTs.push_back(dTT);
		TTime+=dTT;
		TTimes.push_back(TTime);
	//cout << Time << "\t" << dTT << endl;
	//cout << PARss.size() << endl;

	//cout << chaff->DMGreen(1) << endl;

	//actualisation de l'architecture de la plante fait une fois pour toute à l'initialisation
	// Par contre actualisation du rayonnement sur cette architecture à chaque pas de temps

		ComputeProduction(types,PARss,Time,nbsubsteps,TTime); // du coup c'est la production et la demande au temps t-1 qui s'affiche (ainsi que le PAR) (pas obligatoire de mettre TT dans la production)
		ComputeRemobilizedDM(types,dTT,TTime);
		ComputeDemand(types,TTime); // Mettre quelque chose qui fait qu'on sort du programme quand ?l'�chelle de la plante, la demande est nulle

		grain->ActuGrain(TTime,dTT,step,Nmob,(GetConcNmobMin()/nbsubsteps)*DMGreenTot,Production+remobilizedDM,Demand,ConcNmob(Nmob,DMGreenTot));

	//cout << ConcNmob(Nmob,DMGreenTot)-GetConcNmobMin()<< endl;
		ImportPots.push_back(root->ImportRootPot(dTT,Nsoil));

		double ImpTP;
		if (Demand<=0){
			ImpTP=root->ImportRoot(dTT,nbsubsteps,ConcNmob(Nmob,DMGreenTot)-(GetConcNmobMin()/nbsubsteps),0,Nsoil);
			Imports.push_back(ImpTP);
			
			//cout << TTime << "\t" << Nmob << "\t" << ImpTP << "\t" << "\t";
		}
		else{
			ImpTP=root->ImportRoot(dTT,nbsubsteps,ConcNmob(Nmob,DMGreenTot)-(GetConcNmobMin()/nbsubsteps),root->ComputeDemand(TTime)*(Production+remobilizedDM)/Demand,Nsoil);
			Imports.push_back(ImpTP);
			//cout << Nmob << "\t" << root->ComputeDemand(TTime)*(Production+remobilizedDM)/Demand << "\t" << Nsoil << "\t" << endl;
			//cout << Nmob << "\t" << DMGreenTot << endl;
		}
			//cout << root->ImportRoot(dTT,ConcNmob(Nmob,DMGreenTot)-GetConcNmobMin(),root->ComputeDemand(TTime)*Production/Demand,environment->GetNsoil(Time))<<endl; 
		if(Nsoil>root->GetNsoilMin())
		{
			if (Demand<=0) Nsoil += (-1)* root->ImportRoot(dTT,nbsubsteps,ConcNmob(Nmob,DMGreenTot)-(GetConcNmobMin()/nbsubsteps),0,Nsoil);
			else Nsoil += (-1)* root->ImportRoot(dTT,nbsubsteps,ConcNmob(Nmob,DMGreenTot)-(GetConcNmobMin()/nbsubsteps),root->ComputeDemand(TTime)*(Production+remobilizedDM)/Demand,Nsoil);
		}
		//cout << ConcNmob(Nmob,DMGreenTot) << "\t" << (Production+remobilizedDM)/Demand << endl;

	//cout << lamina->NphInput(1,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin(),PARss[LAMINA-1][0]) << endl;
		lamina->ActuNphM2Max(nrow(LAMINA),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[LAMINA-1],Production+remobilizedDM,Demand);
		sheath->ActuNphM2Max(nrow(SHEATH),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[SHEATH-1],Production+remobilizedDM,Demand); 
		internode->ActuNphM2Max(nrow(INTERNODE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[INTERNODE-1],Production+remobilizedDM,Demand);
		peduncle->ActuNphM2Max(nrow(PEDUNCLE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[PEDUNCLE-1],Production+remobilizedDM,Demand);
		chaff->ActuNphM2Max(nrow(CHAFF),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[CHAFF-1],Production+remobilizedDM,Demand);

	//cout << lamina->NphSynthRate(1,dTT,ConcNmob(Nmob,DMGreenTot),PARss[LAMINA-1][0]) << "\t"  << lamina->NphDegRate(1,dTT) << "\t" << (lamina->NphSynthRate(1,dTT,ConcNmob(Nmob,DMGreenTot),PARss[LAMINA-1][0])- lamina->NphDegRate(1,dTT) ) << "\t" << lamina->NphInput(1,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin(),PARss[LAMINA-1][0]) << endl;

		lamina->ActuNphs(nrow(LAMINA),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[LAMINA-1],Production+remobilizedDM,Demand);
		sheath->ActuNphs(nrow(SHEATH),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[SHEATH-1],Production+remobilizedDM,Demand);
		internode->ActuNphs(nrow(INTERNODE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[INTERNODE-1],Production+remobilizedDM,Demand);
		peduncle->ActuNphs(nrow(PEDUNCLE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[PEDUNCLE-1],Production+remobilizedDM,Demand);
		chaff->ActuNphs(nrow(CHAFF),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[CHAFF-1],Production+remobilizedDM,Demand);

		lamina->ActuNstructs(nrow(LAMINA),ConcNmob(Nmob,DMGreenTot),PARss[LAMINA-1],TTime,Production+remobilizedDM,Demand);//voir si je mets ConcNmobMin
		sheath->ActuNstructs(nrow(SHEATH),ConcNmob(Nmob,DMGreenTot),PARss[SHEATH-1],TTime,Production+remobilizedDM,Demand);
		internode->ActuNstructs(nrow(INTERNODE),ConcNmob(Nmob,DMGreenTot),PARss[INTERNODE-1],TTime,Production+remobilizedDM,Demand);
		peduncle->ActuNstructs(nrow(PEDUNCLE),ConcNmob(Nmob,DMGreenTot),PARss[PEDUNCLE-1],TTime,Production+remobilizedDM,Demand);
		chaff->ActuNstructs(nrow(CHAFF),ConcNmob(Nmob,DMGreenTot),PARss[CHAFF-1],TTime,Production+remobilizedDM,Demand);
		root->ActuNrem(ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,TTime,dTT,Production+remobilizedDM,Demand);

		
		ComputeNmob(types,ConcNmob(Nmob,DMGreenTot),PARss,ImpTP,grain->ExportGrain(TTime,dTT,Nmob,GetConcNmobMin()/nbsubsteps*DMGreenTot,Production+remobilizedDM,Demand,ConcNmob(Nmob,DMGreenTot)),root->Ninflux(ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,TTime,dTT,Production+remobilizedDM,Demand));
		//cout << ConcNmob(Nmob,DMGreenTot) << "\t" << ImpTP << endl;
		//ComputeNmob(types,ConcNmob(Nmob,DMGreenTot),PARss,root->ImportRoot(dTT,nbsubsteps,ConcNmob(Nmob,DMGreenTot)-(GetConcNmobMin()/nbsubsteps),root->ComputeDemand(TTime)*(Production+remobilizedDM)/Demand,environment->GetNsoil(Time)),grain->NSynthRate(ConcNmob(Nmob,DMGreenTot),grain->ComputeDemand(TTime)*(Production+remobilizedDM)/Demand));
		// en fait, il faudrait faire l'actualisation de Nmob avant d'actualiser les valeurs de Nph et Nrem !!!

	//temporaire:
	//Forcs.push_back(NphInfluxBalles(parameters->GetInfluxStem(TTime),ConcNmob(Nmob,DMGreenTot),PARss));
	// fin temporaire

		lamina->ActuAreaGreens(nrow(LAMINA),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[LAMINA-1],Production+remobilizedDM,Demand);
		sheath->ActuAreaGreens(nrow(SHEATH),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[SHEATH-1],Production+remobilizedDM,Demand);
		internode->ActuAreaGreens(nrow(INTERNODE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[INTERNODE-1],Production+remobilizedDM,Demand);
		peduncle->ActuAreaGreens(nrow(PEDUNCLE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[PEDUNCLE-1],Production+remobilizedDM,Demand);
		chaff->ActuAreaGreens(nrow(CHAFF),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[CHAFF-1],Production+remobilizedDM,Demand);

	//cout << GAI_cums(1,nrow)[0] << "," << GAI_cums(2,nrow)[0] << "\t" << GAI_cums(1,nrow)[1] << "," << GAI_cums(2,nrow)[1] << "\t" << GAI_cums(1,nrow)[2] << "," << GAI_cums(2,nrow)[2] << "\t" << GAI_cums(1,nrow)[3] << "," << GAI_cums(2,nrow)[3] <<endl;

		lamina->ActuDryMass(nrow(LAMINA),TTime,dTT,Production+remobilizedDM,Demand);
		sheath->ActuDryMass(nrow(SHEATH),TTime,dTT,Production+remobilizedDM,Demand);
		internode->ActuDryMass(nrow(INTERNODE),TTime,dTT,Production+remobilizedDM,Demand);
		peduncle->ActuDryMass(nrow(PEDUNCLE),TTime,dTT,Production+remobilizedDM,Demand);
		chaff->ActuDryMass(nrow(CHAFF),TTime,dTT,Production+remobilizedDM,Demand);
		grain->ActuDryMass(TTime,Production+remobilizedDM,Demand);
		root->ActuDMrem(TTime,dTT,Production+remobilizedDM,Demand); //regarder s'il n'y a pas un artefact par le fait que GreenLab n'autorise pas de r�allocation d'assimilats carbon�s.

		DMgrains.push_back(grain->GetDMgrain());
		Ngrains.push_back(grain->GetNgrain());

	// actualisation du rayonnement, de la production et de la demande sur la nouvelle structure
		for (int unsigned ty=0;ty<types.size();ty++)
		{
	// dimension 1= type de l'organe considéré; dimension 2=type2(=type pour lequel on cumul le GAI);dimension 3= rang de l'organe considéré
			PARss.push_back(environment->PARs(Time,nbsubsteps,types[ty],nrow(types[ty]),rowss_exposed[types[ty]-1],gai_cumsss[ty],lamina->GetInsertionAngle(),HeightViss(types[ty])));
		}
		//cout << (environment->PARs(Time,nbsubsteps,LAMINA,4,rowss_exposed[LAMINA-1],gai_cumsss[0],lamina->GetInsertionAngle()))[0]<< endl;
		PARsss.push_back(PARss);
		ComputeAreaGreenTot(types);
		ComputeDMGreenTot(types);
	}

}


//version avec PAR externe

void Plant::ActuPlantExt(int step, int nbsubsteps,double dTTstep,std::vector<std::vector<double> > _PARss) 
{
  Time+=step;
  for (int j=0;j<nbsubsteps;j++)
    {
      //dTT=environment->ThermalTime(Time,0.,step,nbsubsteps);
      dTT = dTTstep;
      dTTs.push_back(dTT);
      TTime+=dTT;
      TTimes.push_back(TTime);
      
      ComputeProduction(types,PARss,Time,nbsubsteps,TTime); // du coup c'est la production et la demande au temps t-1 qui s'affiche (ainsi que le PAR) (pas obligatoire de mettre TT dans la production)
      ComputeRemobilizedDM(types,dTT,TTime);
      ComputeDemand(types,TTime); // Mettre quelque chose qui fait qu'on sort du programme quand ?l'�chelle de la plante, la demande est nulle
      
      grain->ActuGrain(TTime,dTT,step,Nmob,(GetConcNmobMin()/nbsubsteps)*DMGreenTot,Production+remobilizedDM,Demand,ConcNmob(Nmob,DMGreenTot));
      
      //cout << ConcNmob(Nmob,DMGreenTot)-GetConcNmobMin()<< endl;
      ImportPots.push_back(root->ImportRootPot(dTT,Nsoil));

      double ImpTP;
      if (Demand<=0){
	ImpTP=root->ImportRoot(dTT,nbsubsteps,ConcNmob(Nmob,DMGreenTot)-(GetConcNmobMin()/nbsubsteps),0,Nsoil);
	Imports.push_back(ImpTP);
      }
      else {
	ImpTP=root->ImportRoot(dTT,nbsubsteps,ConcNmob(Nmob,DMGreenTot)-(GetConcNmobMin()/nbsubsteps),root->ComputeDemand(TTime)*(Production+remobilizedDM)/Demand,Nsoil);
	Imports.push_back(ImpTP);
      }
      
      if(Nsoil>root->GetNsoilMin()) {
	if (Demand<=0) 
	  Nsoil += (-1)* root->ImportRoot(dTT,nbsubsteps,ConcNmob(Nmob,DMGreenTot)-(GetConcNmobMin()/nbsubsteps),0,Nsoil);
	else 
	  Nsoil += (-1)* root->ImportRoot(dTT,nbsubsteps,ConcNmob(Nmob,DMGreenTot)-(GetConcNmobMin()/nbsubsteps),root->ComputeDemand(TTime)*(Production+remobilizedDM)/Demand,Nsoil);
      }
      
      lamina->ActuNphM2Max(nrow(LAMINA),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[LAMINA-1],Production+remobilizedDM,Demand);
      sheath->ActuNphM2Max(nrow(SHEATH),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[SHEATH-1],Production+remobilizedDM,Demand); 
      internode->ActuNphM2Max(nrow(INTERNODE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[INTERNODE-1],Production+remobilizedDM,Demand);
      peduncle->ActuNphM2Max(nrow(PEDUNCLE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[PEDUNCLE-1],Production+remobilizedDM,Demand);
      chaff->ActuNphM2Max(nrow(CHAFF),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[CHAFF-1],Production+remobilizedDM,Demand);

	//cout << lamina->NphSynthRate(1,dTT,ConcNmob(Nmob,DMGreenTot),PARss[LAMINA-1][0]) << "\t"  << lamina->NphDegRate(1,dTT) << "\t" << (lamina->NphSynthRate(1,dTT,ConcNmob(Nmob,DMGreenTot),PARss[LAMINA-1][0])- lamina->NphDegRate(1,dTT) ) << "\t" << lamina->NphInput(1,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin(),PARss[LAMINA-1][0]) << endl;

      lamina->ActuNphs(nrow(LAMINA),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[LAMINA-1],Production+remobilizedDM,Demand);
      sheath->ActuNphs(nrow(SHEATH),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[SHEATH-1],Production+remobilizedDM,Demand);
      internode->ActuNphs(nrow(INTERNODE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[INTERNODE-1],Production+remobilizedDM,Demand);
      peduncle->ActuNphs(nrow(PEDUNCLE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[PEDUNCLE-1],Production+remobilizedDM,Demand);
      chaff->ActuNphs(nrow(CHAFF),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[CHAFF-1],Production+remobilizedDM,Demand);

      lamina->ActuNstructs(nrow(LAMINA),ConcNmob(Nmob,DMGreenTot),PARss[LAMINA-1],TTime,Production+remobilizedDM,Demand);//voir si je mets ConcNmobMin
      sheath->ActuNstructs(nrow(SHEATH),ConcNmob(Nmob,DMGreenTot),PARss[SHEATH-1],TTime,Production+remobilizedDM,Demand);
      internode->ActuNstructs(nrow(INTERNODE),ConcNmob(Nmob,DMGreenTot),PARss[INTERNODE-1],TTime,Production+remobilizedDM,Demand);
      peduncle->ActuNstructs(nrow(PEDUNCLE),ConcNmob(Nmob,DMGreenTot),PARss[PEDUNCLE-1],TTime,Production+remobilizedDM,Demand);
      chaff->ActuNstructs(nrow(CHAFF),ConcNmob(Nmob,DMGreenTot),PARss[CHAFF-1],TTime,Production+remobilizedDM,Demand);
      root->ActuNrem(ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,TTime,dTT,Production+remobilizedDM,Demand);

		
      ComputeNmob(types,ConcNmob(Nmob,DMGreenTot),PARss,ImpTP,grain->ExportGrain(TTime,dTT,Nmob,GetConcNmobMin()/nbsubsteps*DMGreenTot,Production+remobilizedDM,Demand,ConcNmob(Nmob,DMGreenTot)),root->Ninflux(ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,TTime,dTT,Production+remobilizedDM,Demand));
		//cout << ConcNmob(Nmob,DMGreenTot) << "\t" << ImpTP << endl;
		//ComputeNmob(types,ConcNmob(Nmob,DMGreenTot),PARss,root->ImportRoot(dTT,nbsubsteps,ConcNmob(Nmob,DMGreenTot)-(GetConcNmobMin()/nbsubsteps),root->ComputeDemand(TTime)*(Production+remobilizedDM)/Demand,environment->GetNsoil(Time)),grain->NSynthRate(ConcNmob(Nmob,DMGreenTot),grain->ComputeDemand(TTime)*(Production+remobilizedDM)/Demand));
		// en fait, il faudrait faire l'actualisation de Nmob avant d'actualiser les valeurs de Nph et Nrem !!!

	//temporaire:
	//Forcs.push_back(NphInfluxBalles(parameters->GetInfluxStem(TTime),ConcNmob(Nmob,DMGreenTot),PARss));
	// fin temporaire

      lamina->ActuAreaGreens(nrow(LAMINA),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[LAMINA-1],Production+remobilizedDM,Demand);
      sheath->ActuAreaGreens(nrow(SHEATH),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[SHEATH-1],Production+remobilizedDM,Demand);
      internode->ActuAreaGreens(nrow(INTERNODE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[INTERNODE-1],Production+remobilizedDM,Demand);
      peduncle->ActuAreaGreens(nrow(PEDUNCLE),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[PEDUNCLE-1],Production+remobilizedDM,Demand);
      chaff->ActuAreaGreens(nrow(CHAFF),TTime,dTT,ConcNmob(Nmob,DMGreenTot),GetConcNmobMin()/nbsubsteps,PARss[CHAFF-1],Production+remobilizedDM,Demand);

	//cout << GAI_cums(1,nrow)[0] << "," << GAI_cums(2,nrow)[0] << "\t" << GAI_cums(1,nrow)[1] << "," << GAI_cums(2,nrow)[1] << "\t" << GAI_cums(1,nrow)[2] << "," << GAI_cums(2,nrow)[2] << "\t" << GAI_cums(1,nrow)[3] << "," << GAI_cums(2,nrow)[3] <<endl;

      lamina->ActuDryMass(nrow(LAMINA),TTime,dTT,Production+remobilizedDM,Demand);
      sheath->ActuDryMass(nrow(SHEATH),TTime,dTT,Production+remobilizedDM,Demand);
      internode->ActuDryMass(nrow(INTERNODE),TTime,dTT,Production+remobilizedDM,Demand);
      peduncle->ActuDryMass(nrow(PEDUNCLE),TTime,dTT,Production+remobilizedDM,Demand);
      chaff->ActuDryMass(nrow(CHAFF),TTime,dTT,Production+remobilizedDM,Demand);
      grain->ActuDryMass(TTime,Production+remobilizedDM,Demand);
      root->ActuDMrem(TTime,dTT,Production+remobilizedDM,Demand); //regarder s'il n'y a pas un artefact par le fait que GreenLab n'autorise pas de r�allocation d'assimilats carbon�s.

      DMgrains.push_back(grain->GetDMgrain());
      Ngrains.push_back(grain->GetNgrain());

	// actualisation du rayonnement, de la production et de la demande sur la nouvelle structure
		//for (int unsigned ty=0;ty<types.size();ty++)
		//{
	// dimension 1= type de l'organe considéré; dimension 2=type2(=type pour lequel on cumul le GAI);dimension 3= rang de l'organe considéré
		//PARss.push_back(environment->PARs(Time,nbsubsteps,types[ty],nrow(types[ty]),rowss_exposed[types[ty]-1],gai_cumsss[ty],lamina->GetInsertionAngle(),HeightViss(types[ty])));
		//}
		//cout << (environment->PARs(Time,nbsubsteps,LAMINA,4,rowss_exposed[LAMINA-1],gai_cumsss[0],lamina->GetInsertionAngle()))[0]<< endl;
      PARsss.push_back(_PARss);
      ComputeAreaGreenTot(types);
      ComputeDMGreenTot(types);
    }
}

void Plant::Affich(std::string outputdir, int nbsubsteps)
{
	lamina->Affich(dTTs,TTimes,"Lamina",nrow(LAMINA),outputdir);
	sheath->Affich(dTTs,TTimes,"Sheath",nrow(SHEATH),outputdir);
	internode->Affich(dTTs,TTimes,"Internode",nrow(INTERNODE),outputdir);
	peduncle->Affich(dTTs,TTimes,"Peduncle",nrow(PEDUNCLE),outputdir);
	chaff->Affich(dTTs,TTimes,"Chaff",nrow(CHAFF),outputdir);
	grain->Affich(TTimes,outputdir);
	root->Affich(TTimes,dTTs,outputdir);

	string buffer;

    buffer=outputdir+"outputs/NmobPred.txt";
	ofstream file1;
	file1.open(buffer.c_str());
	for (unsigned int j=0;j<Nmobs.size();j++)
	{
		file1 << TTimes[j] << "\t" <<Nmobs[j] << "\t" << ConcNmob(Nmobs[j],DMGreenTots[j]) << "\n";
	}
	file1.close();

	buffer=outputdir+"outputs/ProdDemPred.txt";
	ofstream file2;
	file2.open(buffer.c_str());
	for (unsigned int j=0;j<Productions.size();j++)
	{
		file2 << TTimes[j] << "\t" << Productions[j] << "\t" << Demands[j];
		if (Demands[j]<=0) file2 << "\t" << 0 <<"\n";
		else file2 << "\t" << (Productions[j]+remobilizedDMs[j])/Demands[j] <<"\n";
	}
	file2.close();

	buffer=outputdir+"outputs/RootSynthPred.txt";
	ofstream file14;
	file14.open(buffer.c_str());
	for (unsigned int j=0;j<Productions.size();j++)
	{
		file14 << TTimes[j] << "\t" << root->ComputeDemand(TTimes[j])*(Productions[j]+remobilizedDMs[j])/Demands[j] <<"\n";
	}
	file14.close();

	buffer=outputdir+"outputs/AbsorbPred.txt";
	ofstream file3;
	file3.open(buffer.c_str());
	for (unsigned int i=0;i<Imports.size();i++)
	{
		if (Demands[i]<=0) file3 << TTimes[i] << "\t" << Imports[i] << "\t" << root->Neffect(ConcNmob(Nmobs[i],DMGreenTots[i])) << "\t" << root->Ceffect(0,nbsubsteps) << "\t" << ImportPots[i] << "\n";
		else file3 << TTimes[i] << "\t" << Imports[i] << "\t" << root->Neffect(ConcNmob(Nmobs[i],DMGreenTots[i])) << "\t" << root->Ceffect(root->ComputeDemand(TTimes[i])*(Productions[i]+remobilizedDMs[i])/Demands[i],nbsubsteps) << "\t" << ImportPots[i] << "\n";
	}
	file3.close();

	// .csv outputs:
	buffer=outputdir+"outputs/AgreenTot.csv";
	ofstream file4;
	file4.open(buffer.c_str());
	file4 << "dd" << "," << "Total Green Area" << "\n";
	for (unsigned int i=0;i<AreaGreenTots.size();i++)
	{
		file4 << TTimes[i] << "," << AreaGreenTots[i] << "\n";
	}
	file4.close();

	buffer=outputdir+"outputs/Ngrains.csv";
	ofstream file12;
	file12.open(buffer.c_str());
	file12 << "dd" << "," << "Ngrains" << "\n";
	for (unsigned int i=0;i<Ngrains.size();i++)
	{
		file12 << TTimes[i] << "," << Ngrains[i] << "\n";
	}
	file12.close();

	buffer=outputdir+"outputs/RootAbsorbs.csv";
	ofstream file13;
	file13.open(buffer.c_str());
	file13 << "dd" << "," << "RootAbsorbs" << "\n";
	for (unsigned int i=0;i<Imports.size();i++)
	{
		file13 << TTimes[i] << "," << Imports[i] << "\n";
	}
	file13.close();

	buffer=outputdir+"outputs/Production.csv";
	ofstream file11;
	file11.open(buffer.c_str());
	file11 << "dd" << "," << "Total Production" << "\n";
	for (unsigned int i=0;i<Productions.size();i++)
	{
		file11 << TTimes[i] << "," << Productions[i] << "\n";
	}
	file11.close();
	//////////

	buffer=outputdir+"outputs/AgreenTot.txt";
	ofstream file10;
	file10.open(buffer.c_str());
	for (unsigned int i=0;i<AreaGreenTots.size();i++)
	{
		file10 << TTimes[i] << "\t" << AreaGreenTots[i] << "\n";
	}
	file10.close();

	buffer=outputdir+"outputs/PARPredLamina.txt";
	ofstream file5;
	file5.open(buffer.c_str());
	for (unsigned int j=0;j<PARsss.size();j++)
	{
		file5 << TTimes[j];
		for (int r=0;r<nrow(LAMINA);r++) file5 << "\t" << PARsss[j][LAMINA-1][r];
		file5 << "\n";
	}
	file5.close();

	/*cout << CumulGAI(LAMINA,gai_cumsss[LAMINA-1])[0]<< endl;
	cout << CumulGAI(SHEATH,gai_cumsss[SHEATH-1])[0]<< endl;
	cout << CumulGAI(INTERNODE,gai_cumsss[INTERNODE-1])[0]<< endl;
	cout << CumulGAI(PEDUNCLE,gai_cumsss[PEDUNCLE-1])[0]<< endl;
	cout << CumulGAI(CHAFF,gai_cumsss[CHAFF-1])[0]<< endl;*/

	buffer=outputdir+"outputs/GAIPredLamina.txt";
	ofstream file5b;
	file5b.open(buffer.c_str());
	for (unsigned int r=0;r<rowss_exposed[LAMINA-1].size();r++) file5b << "\t" << CumulGAI(LAMINA,gai_cumsss[LAMINA-1])[r];
	file5b << "\n";
	file5b.close();

	buffer=outputdir+"outputs/PARPredSheath.txt";
	ofstream file6;
	file6.open(buffer.c_str());
	for (unsigned int j=0;j<PARsss.size();j++)
	{
		file6 << TTimes[j];
		for (int r=0;r<nrow(SHEATH);r++) file6 << "\t" << PARsss[j][SHEATH-1][r];
		file6 << "\n";
	}
	file6.close();

	buffer=outputdir+"outputs/GAIPredSheath.txt";
	ofstream file6b;
	file6b.open(buffer.c_str());
	for (unsigned int r=0;r<rowss_exposed[SHEATH-1].size();r++) file6b << "\t" << CumulGAI(SHEATH,gai_cumsss[SHEATH-1])[r];
	file6b << "\n";
	file6b.close();

	buffer=outputdir+"outputs/PARPredInternode.txt";
	ofstream file7;
	file7.open(buffer.c_str());
	for (unsigned int j=0;j<PARsss.size();j++)
	{
		file7 << TTimes[j];
		for (int r=0;r<nrow(INTERNODE);r++) file7 << "\t" << PARsss[j][INTERNODE-1][r];
		file7 << "\n";
	}
	file7.close();

	buffer=outputdir+"outputs/GAIPredInternode.txt";
	ofstream file7b;
	file7b.open(buffer.c_str());
	for (unsigned int r=0;r<rowss_exposed[INTERNODE-1].size();r++) file7b << "\t" << CumulGAI(INTERNODE,gai_cumsss[INTERNODE-1])[r];
	file7b << "\n";
	file7b.close();

	buffer=outputdir+"outputs/PARPredPeduncle.txt";
	ofstream file8;
	file8.open(buffer.c_str());
	for (unsigned int j=0;j<PARsss.size();j++)
	{
		file8 << TTimes[j];
		for (int r=0;r<nrow(PEDUNCLE);r++) file8 << "\t" << PARsss[j][PEDUNCLE-1][r];
		file8 << "\n";
	}
	file8.close();

	buffer=outputdir+"outputs/GAIPredPeduncle.txt";
	ofstream file8b;
	file8b.open(buffer.c_str());
	for (unsigned int r=0;r<rowss_exposed[PEDUNCLE-1].size();r++) file8b << "\t" << CumulGAI(PEDUNCLE,gai_cumsss[PEDUNCLE-1])[r];
	file8b << "\n";
	file8b.close();

	buffer=outputdir+"outputs/PARPredChaff.txt";
	ofstream file9;
	file9.open(buffer.c_str());
	for (unsigned int j=0;j<PARsss.size();j++)
	{
		file9 << TTimes[j];
		for (int r=0;r<nrow(CHAFF);r++) file9 << "\t" << PARsss[j][CHAFF-1][r];
		file9 << "\n";
	}
	file9.close();

	buffer=outputdir+"outputs/GAIPredChaff.txt";
	ofstream file9b;
	file9b.open(buffer.c_str());
	for (unsigned int r=0;r<rowss_exposed[CHAFF-1].size();r++) file9b << "\t" << CumulGAI(CHAFF,gai_cumsss[CHAFF-1])[r];
	file9b << "\n";
	file9b.close();

}
//WQL-SA
void Plant::Random(int randUG,ParamPlant* orig_Plant,ParamPlant* Sigmas_Plant,ParamEnv* orig_Env,ParamEnv* Sigmas_Env,ParamGrain* orig_Grain,ParamGrain* Sigmas_Grain,
		ParamRoot* orig_Root,ParamRoot* Sigmas_Root,ParamEntity* orig_lamina,ParamEntity* Sigmas_lamina,ParamEntity* orig_sheath,ParamEntity* Sigmas_sheath,
		ParamEntity* orig_internode,ParamEntity* Sigmas_internode,ParamEntity* orig_peduncle,ParamEntity* Sigmas_peduncle,
		ParamEntity* orig_chaff,ParamEntity* Sigmas_chaff)
{
	this->ParamPlant::Random(randUG,orig_Plant,Sigmas_Plant);
	environment->ParamEnv::Random(randUG,orig_Env,Sigmas_Env);
	grain->ParamGrain::Random(randUG,orig_Grain,Sigmas_Grain);
	root->ParamRoot::Random(randUG,orig_Root,Sigmas_Root);
	
	lamina->ParamEntity::Random(randUG,orig_lamina,Sigmas_lamina);
	sheath->ParamEntity::Random(randUG,orig_sheath,Sigmas_sheath);
	internode->ParamEntity::Random(randUG,orig_internode,Sigmas_internode);
	peduncle->ParamEntity::Random(randUG,orig_peduncle,Sigmas_peduncle);
	chaff->ParamEntity::Random(randUG,orig_chaff,Sigmas_chaff);
}
int Plant::SetNNP(ParamPlant* Sigmas_Plant,ParamEnv* Sigmas_Env,ParamGrain* Sigmas_Grain,ParamRoot* Sigmas_Root,
		ParamEntity* Sigmas_lamina,ParamEntity* Sigmas_sheath,ParamEntity* Sigmas_internode,ParamEntity* Sigmas_peduncle,ParamEntity* Sigmas_chaff)
{
	this->ParamPlant::SetNNP(Sigmas_Plant);
	environment->ParamEnv::SetNNP(Sigmas_Env);
	grain->ParamGrain::SetNNP(Sigmas_Grain);
	root->ParamRoot::SetNNP(Sigmas_Root);

	lamina->ParamEntity::SetNNP(Sigmas_lamina);
	sheath->ParamEntity::SetNNP(Sigmas_sheath);
	internode->ParamEntity::SetNNP(Sigmas_internode);
	peduncle->ParamEntity::SetNNP(Sigmas_peduncle);
	chaff->ParamEntity::SetNNP(Sigmas_chaff);

	NbNNP=this->params_NNP.size()+environment->params_NNP.size()+grain->params_NNP.size()+root->params_NNP.size()+
		lamina->params_NNP.size()+sheath->params_NNP.size()+internode->params_NNP.size()+peduncle->params_NNP.size()+chaff->params_NNP.size();
	
	int temp=0;
	NNP_params_size.push_back(this->params_NNP.size());
	temp+=this->params_NNP.size();
	size_ladder.push_back(temp);
	NNP_params_size.push_back(environment->params_NNP.size());
	temp+=environment->params_NNP.size();
	size_ladder.push_back(temp);
	NNP_params_size.push_back(grain->params_NNP.size());
	temp+=grain->params_NNP.size();
	size_ladder.push_back(temp);
	NNP_params_size.push_back(root->params_NNP.size());
	temp+=root->params_NNP.size();
	size_ladder.push_back(temp);
	NNP_params_size.push_back(lamina->params_NNP.size());
	temp+=lamina->params_NNP.size();
	size_ladder.push_back(temp);
	NNP_params_size.push_back(sheath->params_NNP.size());
	temp+=sheath->params_NNP.size();
	size_ladder.push_back(temp);
	NNP_params_size.push_back(internode->params_NNP.size());
	temp+=internode->params_NNP.size();
	size_ladder.push_back(temp);
	NNP_params_size.push_back(peduncle->params_NNP.size());
	temp+=peduncle->params_NNP.size();
	size_ladder.push_back(temp);
	NNP_params_size.push_back(chaff->params_NNP.size());
	temp+=chaff->params_NNP.size();
	size_ladder.push_back(temp);

	return NbNNP;
}
//Always after use SetNNP function! put all the NNP parameters together in one vetor!
std::vector<double> Plant::SetAllParamsVector()
{
	std::vector<double> A;
	std::vector<std::string> Names;
	string buffer;
	int i;
	
	for(i=0;i<(this->NNP_params_size[0]);i++)
		{
			A.push_back(this->params[this->params_NNP[i]]);
			Names.push_back("Plant_"+this->paramNames[this->params_NNP[i]]);
		}

	for(i=0;i<(this->NNP_params_size[1]);i++)
		{
			A.push_back((this->environment)->params[(this->environment)->params_NNP[i]]);
			Names.push_back("Env_"+(this->environment)->paramNames[(this->environment)->params_NNP[i]]);
		}

	for(i=0;i<(this->NNP_params_size[2]);i++)
		{
			A.push_back((this->grain)->params[(this->grain)->params_NNP[i]]);
			Names.push_back("Grain_"+(this->grain)->paramNames[(this->grain)->params_NNP[i]]);
		}

	for(i=0;i<(this->NNP_params_size[3]);i++)
		{
			A.push_back((this->root)->params[(this->root)->params_NNP[i]]);
			Names.push_back("Root_"+(this->root)->paramNames[(this->root)->params_NNP[i]]);
		}

	for(i=0;i<(this->NNP_params_size[4]);i++)
		{
			A.push_back((this->lamina)->params[(this->lamina)->params_NNP[i]]);
			Names.push_back("Lamina_"+(this->lamina)->paramNames[(this->lamina)->params_NNP[i]]);
		}
	for(i=0;i<(this->NNP_params_size[5]);i++)
		{
			A.push_back((this->sheath)->params[(this->sheath)->params_NNP[i]]);
			Names.push_back("Sheath_"+(this->sheath)->paramNames[(this->sheath)->params_NNP[i]]);
		}
	for(i=0;i<(this->NNP_params_size[6]);i++)
		{
			A.push_back((this->internode)->params[(this->internode)->params_NNP[i]]);
			Names.push_back("Internode_"+(this->internode)->paramNames[(this->internode)->params_NNP[i]]);
		}
	for(i=0;i<(this->NNP_params_size[7]);i++)
		{
			A.push_back((this->peduncle)->params[(this->peduncle)->params_NNP[i]]);
			Names.push_back("Peduncle_"+(this->peduncle)->paramNames[(this->peduncle)->params_NNP[i]]);
		}
	for(i=0;i<(this->NNP_params_size[8]);i++)
		{
			A.push_back((this->chaff)->params[(this->chaff)->params_NNP[i]]);
			Names.push_back("Chaff_"+(this->chaff)->paramNames[(this->chaff)->params_NNP[i]]);
		}
	AllParamVectorNNP=A;
	AllParamNameVectorNNP=Names;
	//A.push_back(1);
	return A;
}

void Plant::Modify_Params(Plant* object,int i)
{
	int  entrance;

	if (0<(i+1)&&(i+1)<=this->size_ladder[0])
		entrance=1; //"plant"
	else
	{
		if (this->size_ladder[0]<(i+1)&&(i+1)<=this->size_ladder[1])
			entrance=2; //"environment";
		else 
		{
			if (this->size_ladder[1]<(i+1)&&(i+1)<=this->size_ladder[2])
				entrance=3; //"grain";
			else
			{
				if (this->size_ladder[2]<(i+1)&&(i+1)<=this->size_ladder[3])
					entrance=4; //"root";
				else
				{
					if (this->size_ladder[3]<(i+1)&&(i+1)<=this->size_ladder[4])
						entrance=5; //"lamina";
					else
					{
						if (this->size_ladder[4]<(i+1)&&(i+1)<=this->size_ladder[5])
							entrance=6; //"sheath";
						else
						{
							if (this->size_ladder[5]<(i+1)&&(i+1)<=this->size_ladder[6])
								entrance=7; //"internode";
							else
							{
								if (this->size_ladder[6]<(i+1)&&(i+1)<=this->size_ladder[7])
									entrance=8; //"peduncle";
								else
								{
									if (this->size_ladder[7]<(i+1)&&(i+1)<=this->size_ladder[8])
										entrance=9; //"chaff";
								}
							}
						}
					}
				}
			}
		}
	}

	switch(entrance)
	{
	    case 1: //"plant":
			{
				((ParamPlant*)this)->params[this->params_NNP[i]]=object->params[this->params_NNP[i]];
				((ParamPlant*)this)->SetParamPlant(this->params[0]);
				break;
			}
		case 2: //"environment":
			{
				i=i-size_ladder[0];
				(this->environment)->params[(this->environment)->params_NNP[i]]=(object->environment)->params[(object->environment)->params_NNP[i]];
				((ParamEnv*)(this->environment))->SetParamEnv((this->environment)->params[0],(this->environment)->params[1]);
				break;
			}
		case 3: //"grain":
			{
				i=i-size_ladder[1];
				(this->grain)->params[(this->grain)->params_NNP[i]]=(object->grain)->params[(object->grain)->params_NNP[i]];
				((ParamGrain*)(this->grain))->SetParamGrain((this->grain)->params[0],(this->grain)->params[1],(this->grain)->params[2],(this->grain)->params[3],(this->grain)->params[4],(this->grain)->params[5],(this->grain)->params[6],(this->grain)->params[7]);
				break;
			}
		case 4: //"root":
			{
				i=i-size_ladder[2];
				(this->root)->params[(this->root)->params_NNP[i]]=(object->root)->params[(object->root)->params_NNP[i]];
				((ParamRoot*)(this->root))->SetParamRoot((this->root)->params[0],(this->root)->params[1],(this->root)->params[2],(this->root)->params[3],(this->root)->params[4],(this->root)->params[5],(this->root)->params[6],
				(this->root)->params[7],(this->root)->params[8],(this->root)->params[9],(this->root)->params[10],(this->root)->params[11],(this->root)->params[12],(this->root)->params[13],(this->root)->params[14]);
				break;
			}
		case 5: //"lamina":
			{
				i=i-size_ladder[3];
				(this->lamina)->params[(this->lamina)->params_NNP[i]]=(object->lamina)->params[(object->lamina)->params_NNP[i]];
				(this->lamina)->SetParamEntity((this->lamina)->params[0],(this->lamina)->params[1],(this->lamina)->params[2],(this->lamina)->params[3],(this->lamina)->params[4],(this->lamina)->params[5],(this->lamina)->params[6],
				(this->lamina)->params[7],(this->lamina)->params[8],(this->lamina)->params[9],(this->lamina)->params[10],(this->lamina)->params[11],(this->lamina)->params[12],(this->lamina)->params[13],(this->lamina)->params[14]);
				break;
			}
		case 6: //"sheath":
			{
				i=i-size_ladder[4];
				(this->sheath)->params[(this->sheath)->params_NNP[i]]=(object->sheath)->params[(object->sheath)->params_NNP[i]];
				(this->sheath)->SetParamEntity((this->sheath)->params[0],(this->sheath)->params[1],(this->sheath)->params[2],(this->sheath)->params[3],(this->sheath)->params[4],(this->sheath)->params[5],(this->sheath)->params[6],
				(this->sheath)->params[7],(this->sheath)->params[8],(this->sheath)->params[9],(this->sheath)->params[10],(this->sheath)->params[11],(this->sheath)->params[12],(this->sheath)->params[13],(this->sheath)->params[14]);
				break;
			}
		case 7: //"internode":
			{
				i=i-size_ladder[5];
				(this->internode)->params[(this->internode)->params_NNP[i]]=(object->internode)->params[(object->internode)->params_NNP[i]];
				(this->internode)->SetParamEntity((this->internode)->params[0],(this->internode)->params[1],(this->internode)->params[2],(this->internode)->params[3],(this->internode)->params[4],(this->internode)->params[5],(this->internode)->params[6],
				(this->internode)->params[7],(this->internode)->params[8],(this->internode)->params[9],(this->internode)->params[10],(this->internode)->params[11],(this->internode)->params[12],(this->internode)->params[13],(this->internode)->params[14]);
				break;
			}
		case 8: //"peduncle":
			{
				i=i-size_ladder[6];
				(this->peduncle)->params[(this->peduncle)->params_NNP[i]]=(object->peduncle)->params[(object->peduncle)->params_NNP[i]];
				(this->peduncle)->SetParamEntity((this->peduncle)->params[0],(this->peduncle)->params[1],(this->peduncle)->params[2],(this->peduncle)->params[3],(this->peduncle)->params[4],(this->peduncle)->params[5],(this->peduncle)->params[6],
				(this->peduncle)->params[7],(this->peduncle)->params[8],(this->peduncle)->params[9],(this->peduncle)->params[10],(this->peduncle)->params[11],(this->peduncle)->params[12],(this->peduncle)->params[13],(this->peduncle)->params[14]);
				break;
			}
		case 9: //"chaff":
			{
				i=i-size_ladder[7];
				(this->chaff)->params[(this->chaff)->params_NNP[i]]=(object->chaff)->params[(object->chaff)->params_NNP[i]];
				(this->chaff)->SetParamEntity((this->chaff)->params[0],(this->chaff)->params[1],(this->chaff)->params[2],(this->chaff)->params[3],(this->chaff)->params[4],(this->chaff)->params[5],(this->chaff)->params[6],
				(this->chaff)->params[7],(this->chaff)->params[8],(this->chaff)->params[9],(this->chaff)->params[10],(this->chaff)->params[11],(this->chaff)->params[12],(this->chaff)->params[13],(this->chaff)->params[14]);
				break;
			}
	}

	/*if (0<(i+1)&&(i+1)<=this->size_ladder[0])
	{
		this->params[this->params_NNP[i]]=object->params[this->params_NNP[i]];
		this->SetParamPlant(this->params[0]);
		
	}

	if (this->size_ladder[0]<(i+1)&&(i+1)<=this->size_ladder[1])
	{
		i=i-size_ladder[0];
		(this->environment)->params[(this->environment)->params_NNP[i]]=(object->environment)->params[(object->environment)->params_NNP[i]];
		this->SetParamEnv((this->environment)->params[0],(this->environment)->params[1]);
		
	}

	if (this->size_ladder[1]<(i+1)&&(i+1)<=this->size_ladder[2])
	{
		i=i-size_ladder[1];
		(this->grain)->params[(this->grain)->params_NNP[i]]=(object->grain)->params[(object->grain)->params_NNP[i]];
		this->SetParamGrain((this->grain)->params[0],(this->grain)->params[1],(this->grain)->params[2],(this->grain)->params[3],(this->grain)->params[4],(this->grain)->params[5]);
		
	}

	if (this->size_ladder[2]<(i+1)&&(i+1)<=this->size_ladder[3])
	{
		i=i-size_ladder[2];
		(this->root)->params[(this->root)->params_NNP[i]]=(object->root)->params[(object->root)->params_NNP[i]];
		this->SetParamRoot((this->root)->params[0],(this->root)->params[1],(this->root)->params[2],(this->root)->params[3],(this->root)->params[4],(this->root)->params[5],(this->root)->params[6],
			 (this->root)->params[7],(this->root)->params[8],(this->root)->params[9],(this->root)->params[10]);
		
	}

	if (this->size_ladder[3]<(i+1)&&(i+1)<=this->size_ladder[4])
	{
		i=i-size_ladder[3];
		(this->lamina)->params[(this->lamina)->params_NNP[i]]=(object->lamina)->params[(object->lamina)->params_NNP[i]];
		(this->lamina)->SetParamEntity((this->lamina)->params[0],(this->lamina)->params[1],(this->lamina)->params[2],(this->lamina)->params[3],(this->lamina)->params[4],(this->lamina)->params[5],(this->lamina)->params[6],
			 (this->lamina)->params[7],(this->lamina)->params[8],(this->lamina)->params[9],(this->lamina)->params[10],(this->lamina)->params[11],(this->lamina)->params[12],(this->lamina)->params[13],(this->lamina)->params[14]);
		
	}

	if (this->size_ladder[4]<(i+1)&&(i+1)<=this->size_ladder[5])
	{
		i=i-size_ladder[4];
		(this->sheath)->params[(this->sheath)->params_NNP[i]]=(object->sheath)->params[(object->sheath)->params_NNP[i]];
		(this->sheath)->SetParamEntity((this->sheath)->params[0],(this->sheath)->params[1],(this->sheath)->params[2],(this->sheath)->params[3],(this->sheath)->params[4],(this->sheath)->params[5],(this->sheath)->params[6],
			 (this->sheath)->params[7],(this->sheath)->params[8],(this->sheath)->params[9],(this->sheath)->params[10],(this->sheath)->params[11],(this->sheath)->params[12],(this->sheath)->params[13],(this->sheath)->params[14]);
		
	}

	if (this->size_ladder[5]<(i+1)&&(i+1)<=this->size_ladder[6])
	{
		i=i-size_ladder[5];
		(this->internode)->params[(this->internode)->params_NNP[i]]=(object->internode)->params[(object->internode)->params_NNP[i]];
		(this->internode)->SetParamEntity((this->internode)->params[0],(this->internode)->params[1],(this->internode)->params[2],(this->internode)->params[3],(this->internode)->params[4],(this->internode)->params[5],(this->internode)->params[6],
			 (this->internode)->params[7],(this->internode)->params[8],(this->internode)->params[9],(this->internode)->params[10],(this->internode)->params[11],(this->internode)->params[12],(this->internode)->params[13],(this->internode)->params[14]);
		
	}
		
	if (this->size_ladder[6]<(i+1)&&(i+1)<=this->size_ladder[7])
	{
		i=i-size_ladder[6];
		(this->peduncle)->params[(this->peduncle)->params_NNP[i]]=(object->peduncle)->params[(object->peduncle)->params_NNP[i]];
		(this->peduncle)->SetParamEntity((this->peduncle)->params[0],(this->peduncle)->params[1],(this->peduncle)->params[2],(this->peduncle)->params[3],(this->peduncle)->params[4],(this->peduncle)->params[5],(this->peduncle)->params[6],
			 (this->peduncle)->params[7],(this->peduncle)->params[8],(this->peduncle)->params[9],(this->peduncle)->params[10],(this->peduncle)->params[11],(this->peduncle)->params[12],(this->peduncle)->params[13],(this->peduncle)->params[14]);
		
	}
		
	if (this->size_ladder[7]<(i+1)&&(i+1)<=this->size_ladder[8])
	{
		i=i-size_ladder[7];
		(this->chaff)->params[(this->chaff)->params_NNP[i]]=(object->chaff)->params[(object->chaff)->params_NNP[i]];
		(this->chaff)->SetParamEntity((this->chaff)->params[0],(this->chaff)->params[1],(this->chaff)->params[2],(this->chaff)->params[3],(this->chaff)->params[4],(this->chaff)->params[5],(this->chaff)->params[6],
			 (this->chaff)->params[7],(this->chaff)->params[8],(this->chaff)->params[9],(this->chaff)->params[10],(this->chaff)->params[11],(this->chaff)->params[12],(this->chaff)->params[13],(this->chaff)->params[14]);
		
	}*/

}

void Plant::Copy(Plant* source)
{	
	//for father class copy
	this->ParamPlant::Copy((ParamPlant*)source);    
	
	//for SA
	NNP_params_size=source->NNP_params_size;
	size_ladder=source->size_ladder;
	NbNNP=source->NbNNP;
	AllParamVectorNNP=source->AllParamVectorNNP;
	//AllParamNameVectorNNP=source->AllParamNameVectorNNP;

	//individual parameters
	Time=source->Time;
	TTime=source->TTime;
	TTimes=source->TTimes;
	dTT=source->dTT;
	dTTs=source->dTTs;
	plt_density=source->plt_density;

	((ParamGrain*)grain)->ParamGrain::Copy((ParamGrain*)(source->grain));
	((ParamRoot*)root)->ParamRoot::Copy((ParamRoot*)(source->root));
	((ParamEnv*)environment)->ParamEnv::Copy((ParamEnv*)(source->environment));

	((ParamEntity*)lamina)->ParamEntity::Copy((ParamEntity*)(source->lamina));
	((ParamEntity*)sheath)->ParamEntity::Copy((ParamEntity*)(source->sheath));
	((ParamEntity*)internode)->ParamEntity::Copy((ParamEntity*)(source->internode));
	((ParamEntity*)peduncle)->ParamEntity::Copy((ParamEntity*)(source->peduncle));
	((ParamEntity*)chaff)->ParamEntity::Copy((ParamEntity*)(source->chaff));
	
	
	types=source->types;
	Nmob=source->Nmob;
	Nmobs=source->Nmobs;
	DMGreenTot=source->DMGreenTot;
	DMGreenTots=source->DMGreenTots;
	Imports=source->Imports;
	
	gai_cumsss=source->gai_cumsss;
	PARss=source->PARss;
	PARsss=source->PARsss;
	AreaGreenTots=source->AreaGreenTots;
	HeightMaxss=source->HeightMaxss;
	visibleAreass=source->visibleAreass;
	rowss_exposed=source->rowss_exposed;

	Production=source->Production; 
	Productions=source->Productions;
	Demand=source->Demand;    
	Demands=source->Demands;

}
//WQL-SA


// END OF File : Plant.cpp

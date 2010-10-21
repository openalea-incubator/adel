/**
 * @file Entity.cpp
 * @brief Implementation of class Entity
 * Copyright: Jessica Bertheloot
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "Entity.h"

// System includes ////////////////////////////////////////// System include //
#include <iostream>
#include <fstream>
#include "math.h"
#include <algorithm>

// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
Entity::Entity()
{
  
  nbrow = 0;

  //vector of parameter values
  params.push_back(0);
  params.push_back(0);
  params.push_back(0);
  params.push_back(0);
  params.push_back(0);
  params.push_back(0);
  //params.push_back(pe2);
  //params.push_back(pe3);
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
  paramNames.push_back("p");
  paramNames.push_back("SynthRate");
  paramNames.push_back("DegRate");
  paramNames.push_back("k1");
  paramNames.push_back("k2");
  paramNames.push_back("peff");
  //paramNames.push_back("peff2");
  //paramNames.push_back("peff3");
  paramNames.push_back("pm1");
  paramNames.push_back("pm2");
  paramNames.push_back("pdeath");
  paramNames.push_back("sink");
  paramNames.push_back("alpha");	
  paramNames.push_back("beta");
  paramNames.push_back("TTexp");
  paramNames.push_back("InsertionAngle");
  paramNames.push_back("DegDMRate");
}

Entity::Entity(std::string fileNameParam, std::string fileNameVar0, int nrow, std::string outputdir):ParamEntity(fileNameParam, fileNameVar0, nrow)
{
	nbrow=nrow;
}

Entity::~Entity()
{
}

////////////////////////////// Private Methods ////////////////////////////////

void Entity::InitialEntity()
{
	//cout << "";
	for (int r=0;r<nbrow;r++)
	{
		// quand extension a la preflo, le push_bach sera une fonction 
		// individuelle (avec comme argument nrow) qui augmentera la taille du vecteur
		Nstructs.push_back(GetNstructs0()[r]); 
		Nphs.push_back(GetNphs0()[r]);
		Lengths.push_back(GetLengths0()[r]);
		Area0s.push_back(GetAreas0()[r]);
		TTinits.push_back(GetTTinits0()[r]);
		DMrems.push_back(GetDMrems0()[r]);
		DMstructs.push_back(GetDMstructs0()[r]);
		//cout << TTinits[r] << "\t"; 
	}
	//cout << GetTTExp() << endl;

	Nstructss.push_back(Nstructs); 
	Nphss.push_back(Nphs);
	TTinitss.push_back(TTinits);
	DMremss.push_back(DMrems);
	DMstructss.push_back(DMstructs);
	// pas s�r ici : il faut initialiser avec l'organogen�se // en supposant que la floraison a lieu ?1000 dd

	std::vector<double> NphSynthRates;
	std::vector<double> NphloemSynthRates;
	std::vector<double> DMremSynthRates;
	for (int r=0;r<nbrow;r++) // temporary for checking the outputs in Affich()
	{
		NphSynthRates.push_back(0.);
		NphloemSynthRates.push_back(0.);
		DMremSynthRates.push_back(0.);
	}
	NphSynthRatess.push_back(NphSynthRates);
	NphloemSynthRatess.push_back(NphloemSynthRates);
	DMremSynthRatess.push_back(DMremSynthRates);

	double a=GetAlpha();
	double b=GetBeta();
	if(a>1 && b>1) // �crire une fonction 
	{
		NormSink=pow((a-1)/(a+b-2),a-1) * pow((b-1)/(a+b-2),b-1);
	}
	else
	{
		NormSink=1.;
		cout << "Attention: coefficient fonction beta < 1"<<endl;
	}
}
void Entity::ResetEntity()
{
	Nstructs.clear();
	Nstructss.clear();
	Nphs.clear();
	Nphss.clear();
	Lengths.clear();
	Area0s.clear();// Area at entity maturity
	TTinits.clear();
	TTinitss.clear();
	NormSink=0;
	DMrems.clear();
	DMremss.clear();
	DMstructs.clear();
	DMstructss.clear();

	AreaGreens.clear();
	AreaGreenss.clear();
	NphM2Maxs.clear();
	NphM2Maxss.clear();

	NphSynthRatess.clear();
	NphloemSynthRatess.clear();
	DMremSynthRatess.clear();
	
	nbrow = 0;

}
//Add Elmer
int Entity::GetNbRow()
{
  return nbrow;
}

std::vector<std::vector<double> > Entity::GetNstructss() {
	return Nstructss;
}

std::vector<std::vector<double> > Entity::GetNphss() {
	return Nphss;
}

std::vector<std::vector<double> > Entity::GetTTinitss() {
	return TTinitss;
}
std::vector<std::vector<double> > Entity::GetDMremss() {
  return DMremss;
}
std::vector<std::vector<double> > Entity::GetDMstructss() {
	return DMstructss;
}
std::vector<std::vector<double> > Entity::GetAreaGreenss() {
	return AreaGreenss;
}
std::vector<std::vector<double> > Entity::GetNphM2Maxss() {
	return NphM2Maxss;
}
std::vector<std::vector<double> > Entity::GetNphSynthRatess() {
	return NphSynthRatess;
}
std::vector<std::vector<double> > Entity::GetNphloemSynthRatess() {
	return NphloemSynthRatess;
}
std::vector<std::vector<double> > Entity::GetDMremSynthRatess() {
	return DMremSynthRatess;
}
std::vector<std::vector<double> > Entity::GetPhotoss() {
	return Photoss;
}
//Add Elmer

double Entity::NstructSynthRate(double ConcNmob,double photo, double dDM)// ne pas l'appeler Nstruct
{
	if (dDM<=photo)
	{
		return 0;
	}else
	{
		return ConcNmob*GetP()*(dDM-photo); 
	}
}


double Entity::NphSynthRate(double dTT, double ConcNmob, double DMgreen, double PAR)
{
	double val=dTT * GetSYNTCOEF() * DMgreen * (ConcNmob / (ConcNmob + GetK1()))*(PAR / (PAR + GetK2()));
	if (val>0)
	{
		return val;
	}else
	{
		return 0;
	}
}

double Entity::NphDegRate(double dTT, double Nph)
{
  return dTT * GetDEGCOEF()*Nph;
}

double Entity::NphM2Min(double Nphm2max, int row)
{
	//return Nphm2max*GetPdeath();
	//cout << row << "\t" << GetNphs0()[row-1]/GetAreas0()[row-1] << endl;
	//if (row==1){
		//cout << "Fe4" << row << "\t" << GetAreas0()[row-1] << endl;
	//	return GetNphs0()[row-1]/GetAreas0()[row-1]*1.5; // c'est un peu du bidouillage pour que la feuille 4 commence à mourir dès la floraison (comme dans les observations)
	//}else
	//{
		//cout << row << "\t" << GetAreas0()[row-1] << endl;
		return GetNphs0()[row-1]/GetAreas0()[row-1]*GetPdeath();
	//}
}

double Entity::DeltaAreaGreen(double Areagreen, double Nph, double Nphinput, double Nphm2max, int row)
{
	if((Nph/Areagreen)<NphM2Min(Nphm2max,row) && Areagreen>0 && Nphinput<0)
	{
		return Nphinput*Areagreen/Nph;
	}else
	{
		return 0;
	}
}

/*double Entity::Peff(double Nph, double Areagreen)
{
	return GetPe1() * (Nph/Areagreen) * GetPe2() / (GetPe1() * (Nph/Areagreen) + GetPe2()) + GetPe3();
}*/

double Entity::Pmax(double Nph, double Areagreen)
{
	return GetPm1() + GetPm2() * (Nph/Areagreen);
}

double Entity::Photo(double Nph, double Areagreen, double PAR)
{
	if (Areagreen>0)
	{
		return GetPeff()*Pmax(Nph,Areagreen)*PAR / (GetPeff()*PAR+Pmax(Nph,Areagreen)) *Areagreen;
	}else
	{
		return 0;
	}
}

double Entity::ComputeSink(double TT, double TTin) // mettre TT ?
{
	double u=(TT-TTin)/GetTTExp();
	if(u<=0||u>=1)
	{
		return 0;
	}
	else
	{
		return GetSink()*pow(u,GetAlpha()-1)*pow(1-u,GetBeta()-1)/NormSink; // voir ?quoi correspond le param�tre Sink
		//return GetSink()*pow(u,GetAlpha()-1)*pow(1-u,GetBeta()-1);
	}
}

double Entity::DMDegRate(double dTT, double dm)
{
	return dTT * GetDEGDMCOEF()*dm;
}

/////////////////////////////// Public Methods ////////////////////////////////

double Entity::DMGreen(int row)
{
	return (DMrems[row-1]+DMstructs[row-1])*AreaGreens[row-1]/Area0s[row-1];//voir si pertinent de mettre DM total
}

double Entity::VisibleLength(int type, int row, double Height, double HeightPreviousEntity)
{
	if (type==LAMINA)
	{
		return Lengths[row-1]; // Here I suppose that the entire lamina length is visible
	}else
	{
		if (Height > HeightPreviousEntity)
		{
			if ( (Height-HeightPreviousEntity)>Lengths[row-1] )
			{
				return Lengths[row-1];
			}else
			{
				return Height-HeightPreviousEntity;
			}
		}else
		{
			return 0;
		}
	}
}


double Entity::VisibleArea(int row, double VisibleLength) // simple cross product
{
	if (VisibleLength>0)
	{
		//cout << Area0s[row-1] << "\t" << VisibleLength << "\t" << Lengths[row-1] << "\t";
		return Area0s[row-1]*VisibleLength/Lengths[row-1];
	}else
	{
		//cout << Area0s[row-1] << "\t" << VisibleLength << "\t" << Lengths[row-1] << "\t";
		return 0;
	}
}

void Entity::InitializeAreaGreensNm2Max(int nrow, std::vector<double> visibleAreas)
{
	for (int r=0;r<nrow;r++)
	{
		AreaGreens.push_back(visibleAreas[r]);
		if (AreaGreens[r]>0)
		{
			NphM2Maxs.push_back((GetNphs0()[r])/(AreaGreens[r]));
		}else
		{
			NphM2Maxs.push_back(0);
		}
		//cout << visibleAreas[r] << "-" << NphM2Maxs[r] << "\t";
	}
	//cout << endl;
	AreaGreenss.push_back(visibleAreas);
	NphM2Maxss.push_back(NphM2Maxs);
}

/*double Entity::NphSynthRate(int row, double dTT, double ConcNmob, double PAR)
{
	double val=dTT * GetSYNTCOEF() * DMGreen(row) * (ConcNmob / (ConcNmob + GetK1()))*(PAR / (PAR + GetK2()));
	if (val>0)
	{
		return val;
	}else
	{
		return 0;
	}
}

double Entity::NphDegRate(int row,double dTT)
{
  return dTT * GetDEGCOEF()*Nphs[row-1];
}*/

double Entity::NphInput(int row, double TT, double dTT, double ConcNmob, double ConcNmobMin, double PAR, double Q, double D)
{
	if ((ConcNmob-ConcNmobMin)>0)
	{
		//cout<<NphSynthRate(dTT,ConcNmob,DMGreen(row),PAR)<<endl;
		if (D<=0) return NphSynthRate(dTT,ConcNmob-ConcNmobMin,DMGreen(row),PAR)-NphDegRate(dTT,Nphs[row-1])+ NstructSynthRate(ConcNmob-ConcNmobMin,Photo(Nphs[row-1],AreaGreens[row-1],PAR),0);
		else return NphSynthRate(dTT,ConcNmob-ConcNmobMin,DMGreen(row),PAR)-NphDegRate(dTT,Nphs[row-1])+ NstructSynthRate(ConcNmob-ConcNmobMin,Photo(Nphs[row-1],AreaGreens[row-1],PAR),ComputeSink(TT, TTinits[row-1])*Q/D);
		//return NphSynthRate(row,dTT,ConcNmob-ConcNmobMin,PAR)-NphDegRate(row,dTT);
	}
	else
	{
		//cout<<NphDegRate(dTT,Nphs[row-1])<<endl;
		return (-1)*NphDegRate(dTT,Nphs[row-1]);
		//return (-1)*NphDegRate(row,dTT);
	}
}

void Entity::ActuNphM2Max(int nrow, double TT, double dTT, double ConcNmob, double ConcNmobMin, std::vector<double> PARs, double Q, double D)
{
	for (int r=0;r<nrow;r++)
	{
		if (AreaGreens[r]>0)
		{
			NphM2Maxs[r]=max(NphM2Maxs[r],(Nphs[r]+NphInput(r+1,TT,dTT,ConcNmob,ConcNmobMin, PARs[r],Q,D))/(AreaGreens[r]+DeltaAreaGreen(AreaGreens[r],Nphs[r],NphInput(r+1,TT,dTT,ConcNmob,ConcNmobMin,PARs[r],Q,D),NphM2Maxs[r],r+1)));
		}
	}
	//for (int r=0;r<(nrow-1);r++) // les rangs en dessous de la feuille drapeau sont supposés avoir une même valeur max que ce que rencontre la feuille drapeau.
	//{ // cela ne fonctionne pas !!
	//	NphM2Maxs[r]= (GetNphs0()[nrow-1])/(GetAreas0()[nrow-1]);
	//}
	NphM2Maxss.push_back(NphM2Maxs);
}

void Entity::ActuNstructs(int nrow, double ConcNmob, std::vector<double> PARs, double TT, double Q, double D) // ?fusionner avec ActuNph pour donner ActuN: attention: ?bien actualiser avant Nph (Cf c'est Nph au pas de temps pr�c�dent qui compte)
{
		for (int r=0;r<nrow;r++)
		{
			//Nstructs[r]+=NstructSynthRate(ConcNmob,Photo(Nphs[r],AreaGreens[r],PARs[r]),ComputeSink(TT, TTinits[r])*Q/D);
			Nstructs[r]+=0;
		}
	Nstructss.push_back(Nstructs);
	//cout<<Nstructs[0]<<"\t"<<Nstructs[1]<<"\t"<<Nstructs[2]<<"\t"<<Nstructs[3]<<endl;
}

void Entity::ActuNphs(int nrow, double TT, double dTT, double ConcNmob, double ConcNmobMin, std::vector<double> PARs, double Q, double D)
{
	std::vector<double> NphSynthRates;
	std::vector<double> NphloemSynthRates;
	for (int r=0;r<nrow;r++)
	{
		//cout << DMGreen(r+1) << endl;
		Nphs[r]+=NphInput(r+1,TT,dTT,ConcNmob,ConcNmobMin,PARs[r],Q,D);
		NphSynthRates.push_back(NphSynthRate(dTT,ConcNmob-ConcNmobMin,DMGreen(r+1),PARs[r]));
		if (D==0) NphloemSynthRates.push_back(NstructSynthRate(ConcNmob-ConcNmobMin,Photo(Nphs[r],AreaGreens[r],PARs[r]),0));
		else NphloemSynthRates.push_back(NstructSynthRate(ConcNmob-ConcNmobMin,Photo(Nphs[r],AreaGreens[r],PARs[r]),ComputeSink(TT, TTinits[r])*Q/D));
	}
	Nphss.push_back(Nphs);
	NphSynthRatess.push_back(NphSynthRates);
	NphloemSynthRatess.push_back(NphloemSynthRates);
	//cout<<Nphs[0]<<endl;
}



void Entity::ActuAreaGreens(int nrow, double TT, double dTT, double ConcNmob,  double ConcNmobMin, std::vector<double> PARs, double Q, double D)
{
	//cout << TT << "\t";
	std::vector<double> Photos;
	for (int r=0;r<nrow;r++)
	{
		AreaGreens[r]+= DeltaAreaGreen(AreaGreens[r],Nphs[r],NphInput(r+1,TT,dTT, ConcNmob,ConcNmobMin,PARs[r],Q,D),NphM2Maxs[r],r+1);// attention: cela prend ici la valeur de Nphm2max actualis�e (cela n'a a piori pas d'importance car la s�nescence se produit apr�s que le maximum soit atteint!)
	//temporaire:
		Photos.push_back(Photo(Nphs[r],AreaGreens[r],PARs[r])/AreaGreens[r]);
	}
	//temporaire:
	Photoss.push_back(Photos);
	//cout << "photo=" << Photo(Nphs[nrow-1],AreaGreens[nrow-1],PARs[nrow-1])/AreaGreens[nrow-1] << "\t" << "Nphm2=" << Nphs[nrow-1]/AreaGreens[nrow-1] << "\t" << "PAR=" << PARs[nrow-1] << "\t";
	//cout << endl;
	AreaGreenss.push_back(AreaGreens);
}

void Entity::ActuDryMass(int nrow, double TT, double dTT, double Q, double D)//il s'agit uniquement de DMrem
{
	// On doit bien avoir Masses.size()=TTinit
	std::vector<double> DMremSynthRates;
	//cout << TT << "\t" << Q << endl;
	for(int r=0; r<nrow; r++)
	{
		if (D==0){
			DMrems[r]+=(-1)*DMDegRate(dTT,DMrems[r]);
			DMremSynthRates.push_back(0);
		}else
		{
			DMrems[r]+=(ComputeSink(TT, TTinits[r])*Q/D)-DMDegRate(dTT,DMrems[r]);
			DMremSynthRates.push_back(ComputeSink(TT, TTinits[r])*Q/D);
		}
	}
	//cout<< GetPm1() << "\t" << GetPm2()<< "\t" << Q << "\t" << ComputeSink(TT, TTinits[0])*Q/D << "\t" << ComputeSink(TT, TTinits[1])*Q/D << "\t" << ComputeSink(TT, TTinits[2])*Q/D << "\t" << ComputeSink(TT, TTinits[3])*Q/D <<endl;
	DMremss.push_back(DMrems);
	DMremSynthRatess.push_back(DMremSynthRates);
	DMstructss.push_back(DMstructs);
}

double Entity::ComputeProduction(int nrow, std::vector<double> PARs, int T, int nbsubsteps, double TT) // voir si je dois mettre TT
{
	double production=0;
	for (int r=0;r<nrow;r++)
	{
		production+= Photo(Nphs[r],AreaGreens[r],PARs[r])/nbsubsteps;
	}
	if (TT>=1000) production=0;// faire quelque chose de plus robuste; ici c'est purement visuel.
	return production;
}

double Entity::ComputeRemobilizedDM(int nrow, double dTT, double TT)
{
	double remobilizedDM=0;
	for (int r=0;r<nrow;r++)
	{
		remobilizedDM+=DMDegRate(dTT,DMrems[r]);
	}
	if (TT>=800) remobilizedDM=0;// faire quelque chose de plus robuste; ici c'est purement visuel.
	return remobilizedDM;
}

double Entity::ComputeDemand(int nrow, double TT) 
{
	double demand=0;
	for(int r=0; r<nrow; r++)
	{
		demand+=ComputeSink(TT, TTinits[r]);
	}
	return demand;
}

// Pour le calcul du GAI dans la class Plant:

/*double Entity::GetAreaGreen0(int row) // area of the exposed part at anthesis
{
	return GetAgreens0()[row-1];
}*/


// For the output file of Plant:
double Entity::GetAreaGreenTot(int nrow)
{
	double AreaGreenTot=0;
	for (int r=0;r<nrow;r++)
	{
		AreaGreenTot+=AreaGreens[r];
	}
	return AreaGreenTot;
}



void Entity::Affich(std::vector<double> dTTs, std::vector<double> TTimes, std::string type,int nrow,std::string outputdir)
{	
	string buffer;
	ofstream file0,file1,file2,file3,file4,file5,file6,file7,file8,file9,file10;

	buffer=outputdir+"outputs/NstructPred";
	string name0=buffer+type+".txt";
	file0.open(name0.c_str());
	for (unsigned int j=0;j<Nstructss.size();j++) 
	{
		file0 << TTimes[j]; 
		for (int r=0;r<nrow;r++) file0 << "\t" << Nstructss[j][r];
		file0 << "\n";
	}
	file0.close();

	buffer=outputdir+"outputs/NphPred";
	string name1=buffer+type+".txt";
	file1.open(name1.c_str());
	for (unsigned int j=0;j<Nphss.size();j++)
	{
		file1 << TTimes[j];
		for (int r=0;r<nrow;r++) file1 << "\t" << Nphss[j][r];
		file1 << "\n";
	}
	file1.close();
	
	buffer=outputdir+"outputs/AgreenPred";
	string name2=buffer+type+".txt";
	file2.open(name2.c_str());
	for (unsigned int j=0;j<Nphss.size();j++)
	{
		file2 << TTimes[j];
		for (int r=0;r<nrow;r++) file2 << "\t" << AreaGreenss[j][r] ;
		file2 << "\n";
	}
	file2.close();

	buffer=outputdir+"outputs/Nphm2Pred";
	string name3=buffer+type+".txt";
	file3.open(name3.c_str());
	for (unsigned int j=0;j<Nphss.size();j++)
	{
		file3 << TTimes[j];
		for (int r=0;r<nrow;r++) file3 << "\t" << Nphss[j][r]/AreaGreenss[j][r] ;
		file3 << "\n";
	}
	file3.close();

	buffer=outputdir+"outputs/DMremPred";
	string name4=buffer+type+".txt";
	file4.open(name4.c_str());
	for (unsigned int j=0;j<DMremss.size();j++)
	{
		file4 << TTimes[j];
		for (int r=0;r<nrow;r++) file4 << "\t" << DMremss[j][r];
		file4 << "\n";
	}
	file4.close();

	buffer=outputdir+"outputs/DMstructPred";
	string name7=buffer+type+".txt";
	file7.open(name7.c_str());
	for (unsigned int j=0;j<DMstructss.size();j++)
	{
		file7 << TTimes[j];
		for (int r=0;r<nrow;r++) file7 << "\t" << DMstructss[j][r];
		file7 << "\n";
	}
	file7.close();


	buffer=outputdir+"outputs/NphDegRatePred";
	string name5=buffer+type+".txt";
	file5.open(name5.c_str());
	for (unsigned int j=0;j<DMremss.size();j++)
	{
		file5 << TTimes[j];
		for (int r=0;r<nrow;r++) file5 << "\t" << NphDegRate(dTTs[j],Nphss[j][r]);
		//for (int r=0;r<nrow;r++) file5 << "\t" << NphDegRate(r+1,dTT);
		file5 << "\n";
	}
	file5.close();


	buffer=outputdir+"outputs/NphSynthRatePred";
	string name6=buffer+type+".txt";
	file6.open(name6.c_str());
	for (unsigned int j=0;j<DMremss.size();j++)
	{
		file6 << TTimes[j];
		for (int r=0;r<nrow;r++) file6 << "\t" << NphSynthRatess[j][r] << "\t" << NphloemSynthRatess[j][r];
		file6 << "\n";
	}
	file6.close();

	buffer=outputdir+"outputs/DMremDegRatePred";
	string name8=buffer+type+".txt";
	file8.open(name8.c_str());
	for (unsigned int j=0;j<DMremss.size();j++)
	{
		file8 << TTimes[j];
		for (int r=0;r<nrow;r++) file8 << "\t" << DMDegRate(dTTs[j],DMremss[j][r]);
		file8 << "\n";
	}
	file5.close();

	buffer=outputdir+"outputs/DMremSynthRatePred";
	string name9=buffer+type+".txt";
	file9.open(name9.c_str());
	for (unsigned int j=0;j<DMremSynthRatess.size();j++)
	{
		file9 << TTimes[j];
		for (int r=0;r<nrow;r++) file9 << "\t" << DMremSynthRatess[j][r];
		file9 << "\n";
	}
	file9.close();

	buffer=outputdir+"outputs/PhotoPred";
	string name10=buffer+type+".txt";
	file10.open(name10.c_str());
	for (unsigned int j=0;j<Photoss.size();j++)
	{
		file10 << TTimes[j];
		for (int r=0;r<nrow;r++) file10 << "\t" << Photoss[j][r];
		file10 << "\n";
	}
	file10.close();

}


// END OF File : Entity.cpp

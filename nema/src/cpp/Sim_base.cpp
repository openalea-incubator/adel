/**
 * @file Sim_base.cpp
 * @brief Implementation of class Sim_base
 * Copyright: Qiongly Wu
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
 */

// Local includes /////////////////////////////////////////// Local Includes //
#include "Sim_base.h"

// System includes ////////////////////////////////////////// System include //
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
using namespace std;

// Namespace to use ///////////////////////////////////////////////////////////
using namespace std;

//////////////////////// Constructor Destructor ... ///////////////////////////
Sim_base::Sim_base()
{
}
Sim_base::Sim_base(std::string fileNameSim)
{
	fsim.open(fileNameSim.c_str());
	this->read_sim();
	fsim.close();
}

Sim_base::~Sim_base()
{
}

///////////////////////public method//////////////////////////////////
int Sim_base::Get_NbSim_i()
{
	return NbSim_i;
}
int Sim_base::Get_NbSim_r()
{
	return NbSim_r;
}

//
std::vector<int> Sim_base::Get_sim_i()
{
	return sim_i;
}
std::vector<double> Sim_base::Get_sim_r()
{
	return sim_r;
}
//	
int Sim_base::Get_randUG()
{
	return randUG;
}
int Sim_base::Get_SY_DM()
{
	return SY_DM;
}
int Sim_base::Get_write_SRC()
{
	return write_SRC;
}
int Sim_base::Get_write_R2()
{
	return write_R2;
}
int Sim_base::Get_write_Y()
{
	return write_Y;
}
int Sim_base::Get_sim_sobS_i()
{
	return sim_sobS_i;
}
int Sim_base::Get_sim_sobST_i()
{
	return sim_sobST_i;
}
int Sim_base::Get_sim_sobS_ij()
{
	return sim_sobS_ij;
}
int Sim_base::Get_sim_sobST_ij()
{
	return sim_sobST_ij;
}
int Sim_base::Get_sim_sig()
{
	return sim_sig;
}
int Sim_base::Get_sim_SRC()
{
	return sim_SRC;
}
int Sim_base::Get_sim_sob()
{
	return sim_sob;
}
//
double Sim_base::Get_sim_eps()
{
	return sim_eps;
}
double Sim_base::Get_sim_srceps()
{
	return sim_srceps;
}
///////////////////////private method////////////////////////////////
void Sim_base::read_sim()
{
	int linelength=200;
	fsim>>randUG;	
	sim_i.push_back(randUG);
	fsim.ignore(linelength, '\n');

	fsim>>SY_DM;	
	sim_i.push_back(SY_DM);
	fsim.ignore(linelength, '\n');

	fsim>>write_SRC;	
	sim_i.push_back(write_SRC);
	fsim.ignore(linelength, '\n');

	fsim>>write_R2;	
	sim_i.push_back(write_R2);
	fsim.ignore(linelength, '\n');

	fsim>>write_Y;	
	sim_i.push_back(write_Y);
	fsim.ignore(linelength, '\n');

	fsim>>sim_sobS_i;	
	sim_i.push_back(sim_sobS_i);
	fsim.ignore(linelength, '\n');

	fsim>>sim_sobST_i;	
	sim_i.push_back(sim_sobST_i);
	fsim.ignore(linelength, '\n');

	fsim>>sim_sobS_ij;	
	sim_i.push_back(sim_sobS_ij);
	fsim.ignore(linelength, '\n');

	fsim>>sim_sobST_ij;	
	sim_i.push_back(sim_sobST_ij);
	fsim.ignore(linelength, '\n');

	fsim>>sim_sig;	
	sim_i.push_back(sim_sig);
	fsim.ignore(linelength, '\n');

	fsim>>sim_SRC;	
	sim_i.push_back(sim_SRC);
	fsim.ignore(linelength, '\n');

	fsim>>sim_sob;	
	sim_i.push_back(sim_sob);
	fsim.ignore(linelength, '\n');

	NbSim_i=sizeof(sim_i);

	fsim>>sim_eps;	
	sim_r.push_back(sim_eps);
	fsim.ignore(linelength, '\n');

	fsim>>sim_srceps;	
	sim_r.push_back(sim_srceps);
	fsim.ignore(linelength, '\n');

	NbSim_r=sizeof(sim_r);

}

void Sim_base::read_sim(std::string fileNameSim)
{
	fsim.open(fileNameSim.c_str());
	this->read_sim();
	fsim.close();
} 
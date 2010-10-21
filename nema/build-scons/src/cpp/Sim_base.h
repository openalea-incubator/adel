/**
 * @file Sim_base.h
 * @brief Declaration of class Sim_base for Sensitivity analysis parameters
 * Copyright: Qiongli Wu
 * INRIA - EPI DigiPlant
 * File created for the digiplant Project.
*/
#ifndef Sim_base_H_
#define Sim_base_H_

// Local includes /////////////////////////////////////////// Local Includes //


// System includes ////////////////////////////////////////// System include //
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Sim_base
{
// Constructor Destructor ... /////////////////////////////////////////////////
public:
	Sim_base();
	//Default constructor

	Sim_base(std::string fileNameSim);
	//Constructor to read the class members from a file

	~Sim_base();
	//Destructor

// Public Methods /////////////////////////////////////////// Public Methods //
public:
	int Get_NbSim_i();
	int Get_NbSim_r();

	std::vector<int> Get_sim_i();
	std::vector<double> Get_sim_r();
	
	int Get_randUG();
	int Get_SY_DM();
	int Get_write_SRC();
	int Get_write_R2();
	int Get_write_Y();
	int Get_sim_sobS_i();
	int Get_sim_sobST_i();
	int Get_sim_sobS_ij();
	int Get_sim_sobST_ij();
	int Get_sim_sig();
	int Get_sim_SRC();
	int Get_sim_sob();

	double Get_sim_eps();
	double Get_sim_srceps();	

// Private Methods ///////////////////////////////////////// Private Methods //
private:
	void read_sim();
	//Reads the class members from fsim;

	void read_sim(std::string fileNameSim);
	//Reads the class members from file;
private:
	int NbSim_i;
	int NbSim_r;
	std::vector<int> sim_i;
	std::vector<double> sim_r;
	int randUG, SY_DM, write_SRC, write_R2, write_Y, sim_sobS_i, sim_sobST_i, sim_sobS_ij,sim_sobST_ij,sim_sig,
	    sim_SRC, sim_sob;
	double sim_eps, sim_srceps;

// Protected /////////////////////////////////////////////////// Protected //
protected:
	std::ifstream fsim;
	//The stream to read the sim_base parameters


};
#endif/*Sim_base_H_*/
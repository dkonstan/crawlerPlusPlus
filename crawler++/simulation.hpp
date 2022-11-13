#include <iostream>
#include <vector>
#include <random>
#include "forces.hpp"
#include "utils.hpp"


#ifndef SIMULATION_HPP
#define SIMULATION_HPP
class Simulation
{
public:
	Simulation(Topology& top, Coordinates& crd, Velocities& vel, Parameters& param, std::vector<Force>& forces, std::string trajFile, std::string logFile, std::string outfile): 
											top(top), crd(crd), vel(vel), param(param), forces(forces), trajectory(trajFile), log(logFile), outfile(outfile) {};
	~Simulation();
	void crawl();
	void writeCoordinates(std::ofstream& output);
private:
	Topology& top;
	Coordinates& crd;
	Parameters& param;
	Velocities& vel;
	std::vector<Force>& forces;
	std::string trajectory;
	std::string log;
	std::string outfile;
	Matrix currForces, newCurrForces;
	double currEnergy, currKineticEnergy, currTemperature;
	double nDegreesOfFreedom;
	std::vector<double> bondLengths;
	std::vector<double> bondAngles;
	std::vector<double> lambdas;
	std::vector<double> lambdasAngle;
	double getTotalEnergy(const Matrix& xyz); // actually potential energy...
	double getKineticEnergy(const Matrix& velxyz);
	Matrix getTotalForces(const Matrix& xyz);
	void integrateVV();
	void integrateLangevin();
	void integrateVVSHAKE();
	void integrateLangevinSHAKE();
	void minimizeEnergy(std::ofstream& logFile, std::ofstream& trajFile);
	void runVVDynamics(std::ofstream& logFile, std::ofstream& trajFile);
	void runLangevinDynamics(std::ofstream& logFile, std::ofstream& trajFile);
	void runVVDynamicsSHAKE(std::ofstream& logFile, std::ofstream& trajFile);
	void runLangevinDynamicsSHAKE(std::ofstream& logFile, std::ofstream& trajFile);
	void runVVDynamicsSHAKEAngles(std::ofstream& logFile, std::ofstream& trajFile);
	void integrateVVSHAKEAngles();
	void runLangevinDynamicsSHAKEAngles(std::ofstream& logFile, std::ofstream& trajFile);
	void integrateLangevinSHAKEAngles();
	void report(std::ofstream& logFile, std::ofstream& trajFile, int step);
	// void runLD(std::ofstream& logFile, std::ofstream &trajFile);
friend std::ostream& operator<<(std::ostream& os, const Simulation& sim);
};



#endif // SIMULATION_HPP

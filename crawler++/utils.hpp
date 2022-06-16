#include <vector>
#include <cmath>
#include <string>


#ifndef UTILS_HPP
#define UTILS_HPP

// #include "parameters.hpp"

typedef std::vector<std::vector<double> > Matrix; // 2D vector, tired of typing the whole thing

const double PI = 3.1415926535;
const double boltzmannK = 3.167e-6; // Hartree / K (https://physics.stackexchange.com/questions/635767/boltzmann-constant-in-atomic-units)

class Topology
{
public:
	Topology(std::string topFile);
	~Topology();
	void readTopology(std::string topFile);
	void checkTopology();
	std::vector<double> masses;
	std::vector<std::vector<int> > bondIdx;
	std::vector<double> bondLengths;
	std::vector<double> bondKs;
	std::vector<std::vector<int> > angleIdx;
	std::vector<double> angleValues;
	std::vector<double> angleThetas;
	std::vector<std::vector<int> > dihedralIdx;
	std::vector<double> dihedralValues;
	std::vector<double> dihedralThetas;
	std::vector<double> charges;
	std::vector<std::vector<int> > vdwIdx;
	std::vector<double> vdwSigmas1;
	std::vector<double> vdwSigmas2;
	std::vector<double> vdwEpsilons;
	std::vector<std::string> atomTypes;
	std::vector<double> box;
	int nAtoms;

friend std::ostream& operator<<(std::ostream& os, const Topology& topol);
};

class Coordinates
{
public:
	Matrix xyz;

	Coordinates(std::string crdFile);
	~Coordinates();
	void readCoordinates(std::string crdFile);
	void checkCoordinates(const Topology& top);
friend std::ostream& operator<<(std::ostream& os, const Coordinates& crd);
};


class Parameters
{
public:
	std::string integrator;
	int nSteps;
	double dt;
	bool constrainBonds;
	double temperature;
	int reportInterval;
	double collisionFreq = 0;

	Parameters(std::string paramFile);
	void readParameters(std::string paramFile);
	void checkParameters();
	~Parameters();
friend std::ostream& operator<<(std::ostream& os, const Parameters& prm);

};

class Velocities
{
public:
	Matrix xyz;

	Velocities(const Topology& top, const Parameters& prm);
	~Velocities();
	void setToTemperature(Topology& top, double temperature);
friend std::ostream& operator<<(std::ostream& os, const Velocities& vel);	
};


// I can't believe std::vector hasn't already overloaded <<, so annoying
std::ostream& operator<<(std::ostream& os, const Matrix& m);
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<int> >& m);
std::ostream& operator<<(std::ostream& os, const std::vector<double>& v);
std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& v);


namespace utils {
	std::vector<double> pbcAdjust(const Topology& top, std::vector<double>& vec);
	double norm(const std::vector<double>& vec);
	void fillForcesCoordinates(Matrix& forces, int nAtoms);
	void addForces(Matrix& forces, const Matrix&& forces2);
	Matrix subVecs(const Matrix& a, const Matrix& b);
	Matrix matMul(const Matrix& a, const Matrix& b);
	Matrix transpose(const Matrix& a);
	double matNorm(const Matrix& a);
	void growFillZeros(Matrix &m, int rows, int cols);
	void fillZeros(Matrix &m);
	Matrix invertMatrix(Matrix &m);
	void growFillEye(Matrix &m, int N);
	Matrix attachEyeToMatrix(Matrix& m, Matrix& id);
	double dot(const std::vector<double>& v1, const std::vector<double>& v2);
	std::vector<double> cross(const std::vector<double>& a, const std::vector<double>& b);
	double sqNorm(std::vector<double>& a, std::vector<double>& b);
	void growFillZeros3D(std::vector<std::vector<std::vector<double> > >& a, int ni, int nj, int nk);
	std::vector<double> flipSign(std::vector<double>& vec);
	std::vector<double> solveEquation(Matrix& mat, std::vector<double>&& vec);
	double max(std::vector<double> vec);
}
#endif // UTILS_HPP
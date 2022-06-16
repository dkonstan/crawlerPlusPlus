#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include "simulation.hpp"
#include "forces.hpp"
#include "utils.hpp"
#include "constraints.hpp"


void addForces(Matrix& to, Matrix&& from)
{
	for (int i = 0; i < to.size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			to[i][j] += from[i][j];	
		}
		
	}
}

Simulation::~Simulation()
{
	;
}


double Simulation::getTotalEnergy(const Matrix& xyz)
{

	double totalEnergy = 0.0;
	for (auto& force : forces)
	{
		totalEnergy += force.energyFunction(top, xyz);
		// if (force.type == "dihedral") break;
	}

	return totalEnergy;
}

double Simulation::getKineticEnergy(const Matrix& velxyz)
{
	double totalKE = 0.0;
	for (int i = 0; i < top.nAtoms; ++i)
	{
		double m = top.masses[i];
		double dotProduct = 0.0;
		for (int j = 0; j < 3; ++j) dotProduct += velxyz[i][j] * velxyz[i][j];
		totalKE += 0.5 * m * dotProduct;
	}

	return totalKE;
}

Matrix Simulation::getTotalForces(const Matrix& xyz)
{
	Matrix totalForces;
	utils::fillForcesCoordinates(totalForces, top.nAtoms);
	for (auto& force : forces)
	{
		utils::addForces(totalForces, force.forceFunction(top, xyz));
		// if (force.type == "dihedral") break;
	}
	return totalForces;
}

std::ostream& operator<<(std::ostream& os, const Simulation& sim)
{
	os << "Simulation:" << std::endl;
	os << "\t" << sim.top.nAtoms << " atoms" << std::endl;
	os << "\tintegrator: " << sim.param.integrator << std::endl;
	os << "\tnSteps: " << sim.param.nSteps << std::endl;
	os << "\tdt: " << sim.param.dt << std::endl;
	os << "\tconstrainBonds: " << sim.param.constrainBonds << std::endl;
	os << "\tforces included:" << std::endl;
	for (int i = 0; i < sim.forces.size(); ++i)
	{
		os << "\t\t" << sim.forces[i].type << std::endl;
	}
	return os;
}

void Simulation::writeCoordinates(std::ofstream& output)
{
	output << top.nAtoms << "\ncomment\n";
	for (int i = 0; i < top.nAtoms; ++i)
	{
		output << top.atomTypes[i] << " " << crd.xyz[i][0] << " " << crd.xyz[i][1] << " " << crd.xyz[i][2] << "\n";
	}
}

void Simulation::minimizeEnergy(std::ofstream& logFile, std::ofstream& trajFile)
{
	double newEnergy;
	int i, j, k, s;

	Matrix newCoords;
	Matrix newForces;
	utils::fillForcesCoordinates(newCoords, top.nAtoms);

	double stepSize = 1e-3; // just to start this off safely

	logFile << "starting minimization..." << std::endl;

	Matrix currForces = getTotalForces(crd.xyz);
	double energy = getTotalEnergy(crd.xyz);

	for (s = 0; s < param.nSteps; ++s)
	{
		logFile << "step: " << s + 1 << " of " << param.nSteps << ", step size: " << stepSize << ", U(x): " << energy << std::endl;
		for (i = 0; i < top.nAtoms; ++i)
		{
			for (j = 0; j < 3; ++j)
			{
				newCoords[i][j] = crd.xyz[i][j] + currForces[i][j] * stepSize;
			}
		}
		newForces = getTotalForces(newCoords);
		newEnergy = getTotalEnergy(newCoords);

		// update stepSize (check this more thoroughly later)
		// stepSize = (utils::matNorm(
		// 	utils::matMul(
		// 		utils::transpose(utils::subVecs(crd.xyz, newCoords)),
		// 		utils::subVecs(newForces, currForces)
		// 		)
		// 	)) /
		// 	(std::pow(utils::matNorm(
		// 		 	  utils::subVecs(newForces, currForces)
		// 		 	  ), 2)
		// 	);
		// if (isnan(stepSize))
		// {
		// 	stepSize = 1e-3;  // think about this
		// 
		crd.xyz = newCoords;
		currForces = newForces;
		energy = newEnergy;
		writeCoordinates(trajFile);
	}
	std::ofstream outFileObj;
	outFileObj.open(outfile);
	writeCoordinates(outFileObj);
}

void Simulation::report(std::ofstream& logFile, std::ofstream& trajFile, int step)
{
	logFile << "step " << step + 1 << " of " << param.nSteps << std::endl;
	logFile << "potential energy: " << currEnergy << std::endl;
	logFile << "kinetic energy: " << currKineticEnergy << std::endl;
	logFile << "total energy: " << currKineticEnergy + currEnergy << std::endl;
	logFile << "temperature: " << currTemperature << "\n" << std::endl;
	writeCoordinates(trajFile);
}

void Simulation::integrateVV()
{
	currForces = getTotalForces(crd.xyz);

	int i, j;
	double m;
	for (i = 0; i < top.nAtoms; ++i)
	{
		m = top.masses[i];
		crd.xyz[i][0] = crd.xyz[i][0] + (param.dt * vel.xyz[i][0]) + ((param.dt * param.dt) / (2 * m)) * currForces[i][0];
		crd.xyz[i][1] = crd.xyz[i][1] + (param.dt * vel.xyz[i][1]) + ((param.dt * param.dt) / (2 * m)) * currForces[i][1];
		crd.xyz[i][2] = crd.xyz[i][2] + (param.dt * vel.xyz[i][2]) + ((param.dt * param.dt) / (2 * m)) * currForces[i][2];
	}

	Matrix newForces = getTotalForces(crd.xyz);

	for (i = 0; i < top.nAtoms; ++i)
	{
		m = top.masses[i];
		vel.xyz[i][0] = vel.xyz[i][0] + (1 / m) * (param.dt / 2) * (currForces[i][0] + newForces[i][0]);
		vel.xyz[i][1] = vel.xyz[i][1] + (1 / m) * (param.dt / 2) * (currForces[i][1] + newForces[i][1]);
		vel.xyz[i][2] = vel.xyz[i][2] + (1 / m) * (param.dt / 2) * (currForces[i][2] + newForces[i][2]) ;
	}

}


double randNormal(double mean, double stddev)
{	
	// Box muller method
	// https://stackoverflow.com/questions/19944111/creating-a-gaussian-random-generator-with-a-mean-and-standard-deviation
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do {
            x = 2.0 * rand() / RAND_MAX - 1;
            y = 2.0 * rand() / RAND_MAX - 1;

            r = x * x + y * y;
        } while (r == 0.0 || r > 1.0);

        {
            double d = sqrt(-2.0 * log(r) / r);
            double n1 = x * d;
            n2 = y * d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2 * stddev + mean;
    }
}

void Simulation::integrateLangevin()
{
	currForces = getTotalForces(crd.xyz);

	int i, j;
	double mi;
	double temp = param.temperature;
	double gamma = param.collisionFreq;
	double dt = param.dt;
	Matrix At, xit, thetat;
	utils::growFillZeros(At, top.nAtoms, 3);
	utils::growFillZeros(xit, top.nAtoms, 3);
	utils::growFillZeros(thetat, top.nAtoms, 3);
	
	std::vector<double> sigmais(top.nAtoms);

	for (i = 0; i < top.nAtoms; ++i)
	{
		mi = top.masses[i];
		sigmais[i] = sqrt(2 * boltzmannK * temp * gamma / mi);
		for (j = 0; j < 3; ++j)
		{
			xit[i][j] = randNormal(0.0, 1.0);
			thetat[i][j] = randNormal(0.0, 1.0);
		}
		for (j = 0; j < 3; ++j)
		{
			At[i][j] = 0.5 * (dt * dt) * (currForces[i][j] - gamma * vel.xyz[i][j]) + sigmais[i] * pow(dt, 1.5) * ((0.5) * xit[i][j] + (1.0 / (2 * sqrt(3))) * thetat[i][j]);
			crd.xyz[i][j] += dt * vel.xyz[i][j] + At[i][j];
		}
	}

 	Matrix newForces = getTotalForces(crd.xyz);

	for (i = 0; i < top.nAtoms; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			vel.xyz[i][j] += (0.5 * dt * (newForces[i][j] + currForces[i][j])) - (dt * gamma * vel.xyz[i][j]) + (sigmais[i] * sqrt(dt) * xit[i][j]) - (gamma * At[i][j]);
		}
	}
}

void Simulation::integrateVVSHAKE()
{
	currForces = getTotalForces(crd.xyz);
	/*
	// make sure bond lengths stay the same (they do)
	std::vector<double> bond(3);
	for (int i = 0; i < top.bondIdx.size(); ++i)
	{

		for (int j = 0; j < 3; ++j)
		{
			bond[j] = crd.xyz[top.bondIdx[i][0]][j] - crd.xyz[top.bondIdx[i][1]][j];
		}
		double norm = utils::norm(bond);
		bondLengths[i] = norm;

	}
	std::cout << bondLengths << std::endl;
	*/

	int i, j;
	double m;
	double dt = param.dt;
	for (i = 0; i < top.nAtoms; ++i)
	{
		m = top.masses[i];
		crd.xyz[i][0] = crd.xyz[i][0] + (dt * vel.xyz[i][0]) + ((dt * dt) / (2 * m)) * currForces[i][0];
		crd.xyz[i][1] = crd.xyz[i][1] + (dt * vel.xyz[i][1]) + ((dt * dt) / (2 * m)) * currForces[i][1];
		crd.xyz[i][2] = crd.xyz[i][2] + (dt * vel.xyz[i][2]) + ((dt * dt) / (2 * m)) * currForces[i][2];
		vel.xyz[i][0] = vel.xyz[i][0] + (1 / m) * (0.5 * dt) * currForces[i][0]; // initial velocity update (dt / 2)
		vel.xyz[i][1] = vel.xyz[i][1] + (1 / m) * (0.5 * dt) * currForces[i][1];
		vel.xyz[i][2] = vel.xyz[i][2] + (1 / m) * (0.5 * dt) * currForces[i][2];
	}

	shakePositions(top, crd.xyz, bondLengths, lambdas); // updates positions and lambdas

	Matrix newForces = getTotalForces(crd.xyz);
	shakeVelocities(top, crd.xyz, vel.xyz, bondLengths, lambdas, dt, newForces);  // updates velocities
}

void Simulation::runVVDynamics(std::ofstream& logFile, std::ofstream& trajFile)
{

	nDegreesOfFreedom = 3 * top.nAtoms;  // 3N (no constraints)

	for (int s = 0; s < param.nSteps; ++s)
	{
		currEnergy = getTotalEnergy(crd.xyz);
		currKineticEnergy = getKineticEnergy(vel.xyz);
		currTemperature = (2.0 * currKineticEnergy) / (nDegreesOfFreedom * boltzmannK);
		if (s % param.reportInterval == 0) report(logFile, trajFile, s);
		integrateVV();
	}
}

void Simulation::runLangevinDynamics(std::ofstream& logFile, std::ofstream &trajFile)
{
	nDegreesOfFreedom = 3 * top.nAtoms; // 3N (no constraints)

	for (int s = 0; s < param.nSteps; ++s)
	{
		currEnergy = getTotalEnergy(crd.xyz);
		currKineticEnergy = getKineticEnergy(vel.xyz);
		currTemperature = (2.0 * currKineticEnergy) / (nDegreesOfFreedom * boltzmannK);
		if (s % param.reportInterval == 0) report(logFile, trajFile, s);
		integrateLangevin();
	}
}

void Simulation::runVVDynamicsSHAKE(std::ofstream& logFile, std::ofstream &trajFile)
{
	nDegreesOfFreedom = 3 * top.nAtoms - top.bondIdx.size();  // 3N - # constraints

	double norm;
	std::vector<double> bond(3);
	for (int i = 0; i < top.bondIdx.size(); ++i)
	{
		lambdas.push_back(0.0);
		
		for (int j = 0; j < 3; ++j)
		{
			bond[j] = crd.xyz[top.bondIdx[i][0]][j] - crd.xyz[top.bondIdx[i][1]][j];
		}
		norm = utils::norm(bond);
		bondLengths.push_back(norm);
	}
	for (int s = 0; s < param.nSteps; ++s)
	{
		currEnergy = getTotalEnergy(crd.xyz);
		currKineticEnergy = getKineticEnergy(vel.xyz);
		currTemperature = (2.0 * currKineticEnergy) / (nDegreesOfFreedom * boltzmannK);
		if (s % param.reportInterval == 0) report(logFile, trajFile, s);
		integrateVVSHAKE();
	}
}

void Simulation::crawl()
{
	std::ofstream logFile;
	logFile.open(log);
	std::ofstream trajFile;
	trajFile.open(trajectory);
	const std::string beginMessage = "--- CRAWLER++ will now crawl. ---";
	const std::string endMessage = "---- CRAWLER++ has crawled.---\n";

	std::cout << beginMessage << std::endl;

	if (param.integrator == "min")
	{
		std::cout << "minimizing energy..." << std::endl;
		minimizeEnergy(logFile, trajFile);
	}
	else if (param.integrator == "vv")
	{
		std::cout << "running dynamics (NVE)..." << std::endl;
		if (param.constrainBonds) runVVDynamicsSHAKE(logFile, trajFile);
		else runVVDynamics(logFile, trajFile);
	}
	else if (param.integrator == "langevin")
	{
		std::cout << "running dynamics (NVT)..." << std::endl;
		if (param.constrainBonds);
		else runLangevinDynamics(logFile, trajFile);
	}

	logFile << endMessage << std::endl;
	std::cout << endMessage << std::endl;

	trajFile.close();
	logFile.close();
}



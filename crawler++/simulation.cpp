#include <vector>
#include <cmath>
#include <fstream>
#include "simulation.hpp"
#include "forces.hpp"
#include "utils.hpp"

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

	double stepSize = 1e-5; // just to start this off safely

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
		stepSize = (utils::matNorm(
			utils::matMul(
				utils::transpose(utils::subVecs(crd.xyz, newCoords)),
				utils::subVecs(newForces, currForces)
				)
			)) /
			(std::pow(utils::matNorm(
				 	  utils::subVecs(newForces, currForces)
				 	  ), 2)
			);
		if (isnan(stepSize))
		{
			stepSize = 1e-3;  // think about this
		}
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

	newCurrForces = getTotalForces(crd.xyz);

	for (i = 0; i < top.nAtoms; ++i)
	{
		m = top.masses[i];
		vel.xyz[i][0] = vel.xyz[i][0] + (1 / m) * (param.dt / 2) * (currForces[i][0] + newCurrForces[i][0]);
		vel.xyz[i][1] = vel.xyz[i][1] + (1 / m) * (param.dt / 2) * (currForces[i][1] + newCurrForces[i][1]);
		vel.xyz[i][2] = vel.xyz[i][2] + (1 / m) * (param.dt / 2) * (currForces[i][2] + newCurrForces[i][2]) ;
	}

}
void Simulation::runVVDynamics(std::ofstream& logFile, std::ofstream& trajFile)
{

	nDegreesOfFreedom = 3 * top.nAtoms;  // 3N (no constraints)
	for (int s = 0; s < param.nSteps; ++s)
	{
		currEnergy = getTotalEnergy(crd.xyz);
		currKineticEnergy = getKineticEnergy(vel.xyz);
		currTemperature = (2.0 * currKineticEnergy) / (nDegreesOfFreedom * boltzmannK);

		if (s % param.reportInterval == 0)
		{
			report(logFile, trajFile, s);
		}
		integrateVV();
	}
}

void Simulation::crawl()
{
	std::ofstream logFile;
	logFile.open(log);
	std::ofstream trajFile;
	trajFile.open(trajectory);

	if (param.integrator == "min")
	{	
		minimizeEnergy(logFile, trajFile);
	}
	else if (param.integrator == "vv")
	{
		if (param.constrainBonds)
		{
			std::cout << "constaining bonds placeholder" << std::endl;
		}
		else
		{
			runVVDynamics(logFile, trajFile);
		}
	}
	// else if (param.integrator == "langevin")
	// {
	// 	if (param.constrainBonds)
	// 	{
	// 		;
	// 	}
	// 	else
	// 	{
	// 		;
	// 	}
	// }

	const std::string endMessage = "\n\n---- CRAWLER++ has crawled.---";
	logFile << endMessage << std::endl;
	std::cout << endMessage << std::endl;

	trajFile.close();
	logFile.close();
}



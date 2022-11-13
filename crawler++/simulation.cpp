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
		for (int j = 0; j < 3; ++j)
		{
			dotProduct += velxyz[i][j] * velxyz[i][j];
		}
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
		// for output put it back into Å from meters
		output << top.atomTypes[i] << " " << crd.xyz[i][0] * 1e10 << " " << crd.xyz[i][1] * 1e10 << " " << crd.xyz[i][2] * 1e10 << "\n";
	}
}

void Simulation::minimizeEnergy(std::ofstream& logFile, std::ofstream& trajFile)
{
	double newEnergy;
	int i, j, k, s;

	Matrix newCoords;
	Matrix newForces;
	utils::fillForcesCoordinates(newCoords, top.nAtoms);

	double stepSize = 0.001; // just to start this off safely

	logFile << "starting minimization..." << std::endl;

	// std::cout << top << std::endl;
	// exit(0);
	currForces = getTotalForces(crd.xyz);
	// exit(0);
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
		// std::cout << currForces[0] << std::endl;

		energy = newEnergy;
		if (s % 100 == 0)
		{
			std::cout << energy << std::endl;
			writeCoordinates(trajFile);	
		}
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


void Simulation::integrateVVSHAKE()
{

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

	shakePositions(top, crd.xyz, bondLengths, lambdas);
	shakeVelocities(top, crd.xyz, vel.xyz, bondLengths, lambdas, param.dt, newForces);

	currForces = newForces;
}


void Simulation::integrateVVSHAKEAngles()
{

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
	shakeAngles(top, crd.xyz, bondAngles, lambdasAngle);
	shakePositions(top, crd.xyz, bondLengths, lambdas);
	shakeAnglesVelocities(top, crd.xyz, vel.xyz, bondAngles, lambdasAngle, param.dt, newForces);
	shakeVelocities(top, crd.xyz, vel.xyz, bondLengths, lambdas, param.dt, newForces);

	currForces = newForces;

	std::vector<double> cog(3, 0.0);
	for (i = 0; i < top.nAtoms; ++i)
	{
		cog[0] += crd.xyz[i][0];
		cog[1] += crd.xyz[i][1];
		cog[2] += crd.xyz[i][2];
	}

	cog[0] /= top.nAtoms;
	cog[1] /= top.nAtoms;
	cog[2] /= top.nAtoms;
	for (i = 0; i < top.nAtoms; ++i)
	{
		crd.xyz[i][0] -= cog[0];
		crd.xyz[i][1] -= cog[1];
		crd.xyz[i][2] -= cog[2];
	}
}

void Simulation::integrateLangevinSHAKEAngles()
{
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
			xit[i][j] = utils::randNormal(0.0, 1.0);
			thetat[i][j] = utils::randNormal(0.0, 1.0);
		}
		for (j = 0; j < 3; ++j)
		{
			At[i][j] = 0.5 * pow(dt, 2) * (currForces[i][j] / mi - gamma * vel.xyz[i][j]) + sigmais[i] * pow(dt, 1.5) * ((0.5) * xit[i][j] + (1.0 / (2 * sqrt(3))) * thetat[i][j]);
			crd.xyz[i][j] += dt * vel.xyz[i][j] + At[i][j];
		}
	}

	shakeAngles(top, crd.xyz, bondAngles, lambdasAngle);
	shakePositions(top, crd.xyz, bondLengths, lambdas);

 	Matrix newForces = getTotalForces(crd.xyz);

	for (i = 0; i < top.nAtoms; ++i)
	{
		mi = top.masses[i];
		for (j = 0; j < 3; ++j)
		{
			vel.xyz[i][j] += (0.5 * dt * (1 / mi) * (newForces[i][j] + currForces[i][j])) - (dt * gamma * vel.xyz[i][j]) + (sigmais[i] * sqrt(dt) * xit[i][j]) - (gamma * At[i][j]);
		}
	}

	shakeAnglesVelocities(top, crd.xyz, vel.xyz, bondAngles, lambdasAngle, param.dt, newForces);
	shakeVelocities(top, crd.xyz, vel.xyz, bondLengths, lambdas, param.dt, newForces);

	currForces = newForces;

	std::vector<double> cog(3, 0.0);
	std::vector<double> cov(3, 0.0);
	for (i = 0; i < top.nAtoms; ++i)
	{
		cog[0] += crd.xyz[i][0];
		cog[1] += crd.xyz[i][1];
		cog[2] += crd.xyz[i][2];
		cov[0] += vel.xyz[i][0];
		cov[1] += vel.xyz[i][1];
		cov[2] += vel.xyz[i][2];
	}

	cog[0] /= top.nAtoms;
	cog[1] /= top.nAtoms;
	cog[2] /= top.nAtoms;
	cov[0] /= top.nAtoms;
	cov[1] /= top.nAtoms;
	cov[2] /= top.nAtoms;
	for (i = 0; i < top.nAtoms; ++i)
	{
		crd.xyz[i][0] -= cog[0];
		crd.xyz[i][1] -= cog[1];
		crd.xyz[i][2] -= cog[2];
		vel.xyz[i][0] -= cov[0];
		vel.xyz[i][1] -= cov[1];
		vel.xyz[i][2] -= cov[2];
	}
}


void Simulation::integrateLangevin()
{
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

	// J = kg-m^2/s^2
	// N = kg-m/s^2
	for (i = 0; i < top.nAtoms; ++i)
	{
		mi = top.masses[i];
		sigmais[i] = sqrt(2 * boltzmannK * temp * gamma / mi);
		for (j = 0; j < 3; ++j)
		{
			xit[i][j] = utils::randNormal(0.0, 1.0);
			thetat[i][j] = utils::randNormal(0.0, 1.0);
		}
		for (j = 0; j < 3; ++j)
		{
			At[i][j] = 0.5 * pow(dt, 2) * (currForces[i][j] / mi - gamma * vel.xyz[i][j]) + sigmais[i] * pow(dt, 1.5) * ((0.5) * xit[i][j] + (1.0 / (2 * sqrt(3))) * thetat[i][j]);
			crd.xyz[i][j] += dt * vel.xyz[i][j] + At[i][j];
		}
	}

 	Matrix newForces = getTotalForces(crd.xyz);

	for (i = 0; i < top.nAtoms; ++i)
	{
		gamma = param.collisionFreq;
		mi = top.masses[i];
		for (j = 0; j < 3; ++j)
		{
			vel.xyz[i][j] += (0.5 * dt * (1 / mi) * (newForces[i][j] + currForces[i][j])) - (dt * gamma * vel.xyz[i][j]) + (sigmais[i] * sqrt(dt) * xit[i][j]) - (gamma * At[i][j]);
		}
	}

	currForces = newForces;

	std::vector<double> cog(3, 0.0);
	for (i = 0; i < top.nAtoms; ++i)
	{
		cog[0] += crd.xyz[i][0];
		cog[1] += crd.xyz[i][1];
		cog[2] += crd.xyz[i][2];
	}

	cog[0] /= top.nAtoms;
	cog[1] /= top.nAtoms;
	cog[2] /= top.nAtoms;
	for (i = 0; i < top.nAtoms; ++i)
	{
		crd.xyz[i][0] -= cog[0];
		crd.xyz[i][1] -= cog[1];
		crd.xyz[i][2] -= cog[2];
	}
}


void Simulation::runVVDynamics(std::ofstream& logFile, std::ofstream& trajFile)
{

	nDegreesOfFreedom = 3 * top.nAtoms;  // 3N (no constraints)
	currForces = getTotalForces(crd.xyz);

	for (int s = 0; s < param.nSteps; ++s)
	{
		currEnergy = getTotalEnergy(crd.xyz);
		currKineticEnergy = getKineticEnergy(vel.xyz);
		currTemperature = (2.0 * currKineticEnergy) / (nDegreesOfFreedom * boltzmannK);
		if (s % param.reportInterval == 0) report(logFile, trajFile, s);
		integrateVV();
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

	currForces = getTotalForces(crd.xyz);
	for (int s = 0; s < param.nSteps; ++s)
	{
		currEnergy = getTotalEnergy(crd.xyz);
		currKineticEnergy = getKineticEnergy(vel.xyz);
		currTemperature = (2.0 * currKineticEnergy) / (nDegreesOfFreedom * boltzmannK);
		if (s % param.reportInterval == 0) report(logFile, trajFile, s);
		integrateVVSHAKE();
	}
}


void Simulation::runVVDynamicsSHAKEAngles(std::ofstream& logFile, std::ofstream& trajFile)
{
	nDegreesOfFreedom = 3 * top.nAtoms - top.bondIdx.size() - top.angleIdx.size();  // 3N - # constraints

	double norm;
	std::vector<double> bond(3), bond1(3), bond2(3);
	for (int i = 0; i < top.bondIdx.size(); ++i)
	{
		lambdas.push_back(0.0);
		
		for (int j = 0; j < 3; ++j)
		{
			bond[j] = crd.xyz[top.bondIdx[i][0]][j] - crd.xyz[top.bondIdx[i][1]][j];
		}
		bond = utils::pbcAdjust(top, bond);
		bondLengths.push_back(utils::norm(bond));
	}

	for (int i = 0; i < top.angleIdx.size(); ++i)
	{
		lambdasAngle.push_back(0.0);
		bondAngles.push_back(utils::getAngle(top, crd.xyz[top.angleIdx[i][0]], crd.xyz[top.angleIdx[i][1]], crd.xyz[top.angleIdx[i][2]]));
	}
	currForces = getTotalForces(crd.xyz);
	for (int s = 0; s < param.nSteps; ++s)
	{
		currEnergy = getTotalEnergy(crd.xyz);
		currKineticEnergy = getKineticEnergy(vel.xyz);
		currTemperature = (2.0 * currKineticEnergy) / (nDegreesOfFreedom * boltzmannK);
		if (s % param.reportInterval == 0) report(logFile, trajFile, s);
		integrateVVSHAKEAngles();
	}
}


void Simulation::runLangevinDynamics(std::ofstream& logFile, std::ofstream &trajFile)
{
	nDegreesOfFreedom = 3 * top.nAtoms; // 3N (no constraints)

	currForces = getTotalForces(crd.xyz);
	for (int s = 0; s < param.nSteps; ++s)
	{
		currEnergy = getTotalEnergy(crd.xyz);
		currKineticEnergy = getKineticEnergy(vel.xyz);
		currTemperature = (2.0 * currKineticEnergy) / (nDegreesOfFreedom * boltzmannK);
		if (s % param.reportInterval == 0) report(logFile, trajFile, s);
		integrateLangevin();
	}
}


void Simulation::runLangevinDynamicsSHAKEAngles(std::ofstream& logFile, std::ofstream& trajFile)
{
	currForces = getTotalForces(crd.xyz);

	nDegreesOfFreedom = 3 * top.nAtoms - top.bondIdx.size() - top.angleIdx.size() - 3;  // 3N - # constraints - 3 (for com stationary)

	double norm;
	std::vector<double> bond(3), bond1(3), bond2(3);
	for (int i = 0; i < top.bondIdx.size(); ++i)
	{
		lambdas.push_back(0.0);
		
		for (int j = 0; j < 3; ++j)
		{
			bond[j] = crd.xyz[top.bondIdx[i][0]][j] - crd.xyz[top.bondIdx[i][1]][j];
		}
		bond = utils::pbcAdjust(top, bond);
		bondLengths.push_back(utils::norm(bond));
	}

	for (int i = 0; i < top.angleIdx.size(); ++i)
	{
		lambdasAngle.push_back(0.0);
		bondAngles.push_back(utils::getAngle(top, crd.xyz[top.angleIdx[i][0]], crd.xyz[top.angleIdx[i][1]], crd.xyz[top.angleIdx[i][2]]));
	}

	currForces = getTotalForces(crd.xyz);
	for (int s = 0; s < param.nSteps; ++s)
	{
		currEnergy = getTotalEnergy(crd.xyz);
		currKineticEnergy = getKineticEnergy(vel.xyz);
		currTemperature = (2.0 * currKineticEnergy) / (nDegreesOfFreedom * boltzmannK);
		if (s % param.reportInterval == 0) report(logFile, trajFile, s);
		integrateLangevinSHAKEAngles();
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
		if (param.constrainBonds && !param.constrainAngles)
		{
			runVVDynamicsSHAKE(logFile, trajFile);
		}
		else if (param.constrainBonds && param.constrainAngles)
		{
			runVVDynamicsSHAKEAngles(logFile, trajFile);
		}
		else if (!param.constrainBonds && param.constrainAngles)
		{
			std::cout << "can't constrain angles without constraining bonds!" << std::endl;
		}
		else
		{
			runVVDynamics(logFile, trajFile);
		}
	}
	else if (param.integrator == "langevin")
	{
		std::cout << "running dynamics (NVT)..." << std::endl;
		if (param.constrainBonds && !param.constrainAngles)
		{
			// runLangevinDynamicsSHAKE(logFile, trajFile);
		}
		else if (param.constrainBonds && param.constrainAngles)
		{
			runLangevinDynamicsSHAKEAngles(logFile, trajFile);
		}
		else if (!param.constrainBonds && param.constrainAngles)
		{
			std::cout << "can't constrain angles without constraining bonds!" << std::endl;
		}
		else
		{
			runLangevinDynamics(logFile, trajFile);
		}
	}

	logFile << endMessage << std::endl;
	std::cout << endMessage << std::endl;

	trajFile.close();
	logFile.close();
}


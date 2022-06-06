#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include "utils.hpp"



inline int nthSubstr(int n, const std::string& s, const std::string& p) {
	/* https://www.oreilly.com/library/view/c-cookbook/0596007612/ch04s11.html */

   std::string::size_type i = s.find(p);     // Find the first occurrence

   int j;
   for (j = 1; j < n && i != std::string::npos; ++j)
      i = s.find(p, i + 1); // Find the next occurrence

   if (j == n) {
     return(i);
   }
   else {
     return(-1);
   }
}

Topology::Topology(std::string topFile)
{
	readTopology(topFile);
}

Topology::~Topology()
{
	;
}

void Topology::checkTopology()
{
	if ((masses.size() != charges.size()) ||
		(masses.size() != nAtoms) ||
			(charges.size() != nAtoms) ||
			(atomTypes.size() != masses.size()) ||
			(atomTypes.size() != charges.size()))
	{
		std::cout << "ERROR: inconsistent number of atoms in topology" << std::endl;
		exit(1);
	}
}

void Topology::readTopology(std::string topFile)
{
	std::ifstream top(topFile);
	if (!top.good()) {
		std::cout << "can't open topology file" << std::endl;
		exit(1);
	}
	std::string line, atomType, number;
	// int bondNumber1, bondNumber2, angleNumber1, angleNumber2, angleNumber3;
	// double bondNumber3, bondNumber4, angleNumber4, angleNumber5;
	std::vector<int> bondIdxRow, angleIdxRow, dihedralIdxRow, vdwIdxRow;
	std::string delim = " ";

	try {
		while(std::getline(top, line))
		{
			if (line == "<atomtypes>")
			{
				std::getline(top, line);
				while (line != "<end>")
				{
					number = line.substr(0, line.find(delim));
					atomType = line.substr(number.length() + delim.length(), line.length() - line.find(delim));
					atomTypes.push_back(atomType);
					std::getline(top, line);
				}
			}

			if (line == "<masses>")
			{
				std::getline(top, line);
				while (line != "<end>")
				{
					std::string number = line.substr(0, line.find(delim));
					double mass = std::stod(line.substr(number.length() + delim.length(), line.length() - line.find(delim)));
					masses.push_back(mass);
					std::getline(top, line);
				}
			}

			if (line == "<charges>")
			{
				std::getline(top, line);
				while (line != "<end>")
				{
					std::string number = line.substr(0, line.find(delim));
					double charge = std::stod(line.substr(number.length() + delim.length(), line.length() - line.find(delim)));
					charges.push_back(charge);
					std::getline(top, line);
				}
			}

			if (line == "<bonds>")
			{
				std::getline(top, line);
				while (line != "<end>")
				{
					int bondNumber1 = std::stoi(line.substr(0, line.find(delim)));
					int bondNumber2 = std::stoi(line.substr(nthSubstr(1, line, delim)));
					double bondNumber3 = std::stod(line.substr(nthSubstr(2, line, delim)));
					double bondNumber4 = std::stod(line.substr(nthSubstr(3, line, delim)));
					bondIdxRow.push_back(bondNumber1);
					bondIdxRow.push_back(bondNumber2);
					bondIdx.push_back(bondIdxRow);
					bondLengths.push_back(bondNumber3);
					bondKs.push_back(bondNumber4);
					bondIdxRow.clear();
					std::getline(top, line);
				}
			}

			if (line == "<angles>")
			{
				std::getline(top, line);
				while (line != "<end>")
				{
					int angleNumber1 = std::stoi(line.substr(0, line.find(delim)));
					int angleNumber2 = std::stoi(line.substr(nthSubstr(1, line, delim)));
					int angleNumber3 = std::stoi(line.substr(nthSubstr(2, line, delim)));
					double angleNumber4 = std::stod(line.substr(nthSubstr(3, line, delim)));
					double angleNumber5 = std::stod(line.substr(nthSubstr(4, line, delim)));
					angleIdxRow.push_back(angleNumber1);
					angleIdxRow.push_back(angleNumber2);
					angleIdxRow.push_back(angleNumber3);
					angleIdx.push_back(angleIdxRow);
					angleValues.push_back(angleNumber4 * (PI / 180.0)); // convert to radians for internal use
					angleThetas.push_back(angleNumber5);
					angleIdxRow.clear();
					std::getline(top, line);
				}
			}

			if (line == "<dihedrals>")
			{
				std::getline(top, line);
				while (line != "<end>")
				{
					int dihedralNumber1 = std::stoi(line.substr(0, line.find(delim)));
					int dihedralNumber2 = std::stoi(line.substr(nthSubstr(1, line, delim)));
					int dihedralNumber3 = std::stoi(line.substr(nthSubstr(2, line, delim)));
					int dihedralNumber4 = std::stoi(line.substr(nthSubstr(3, line, delim)));
					double dihedralNumber5 = std::stod(line.substr(nthSubstr(4, line, delim)));
					double dihedralNumber6 = std::stod(line.substr(nthSubstr(5, line, delim)));
					dihedralIdxRow.push_back(dihedralNumber1);
					dihedralIdxRow.push_back(dihedralNumber2);
					dihedralIdxRow.push_back(dihedralNumber3);
					dihedralIdxRow.push_back(dihedralNumber4);
					dihedralIdx.push_back(dihedralIdxRow);
					dihedralValues.push_back(dihedralNumber5 * (PI / 180.0));  // convert to radians for internal use
					dihedralThetas.push_back(dihedralNumber6);
					dihedralIdxRow.clear();
					std::getline(top, line);
				}
			}
			if (line == "<vdw>")
			{
				std::getline(top, line);
				while (line != "<end>")
				{
					int vdwNumber1 = std::stoi(line.substr(0, line.find(delim)));
					int vdwNumber2 = std::stoi(line.substr(nthSubstr(1, line, delim)));
					double vdwNumber3 = std::stod(line.substr(nthSubstr(2, line, delim)));
					double vdwNumber4 = std::stod(line.substr(nthSubstr(3, line, delim)));
					double vdwNumber5 = std::stod(line.substr(nthSubstr(4, line, delim)));
					
					vdwIdxRow.push_back(vdwNumber1);
					vdwIdxRow.push_back(vdwNumber2);
					vdwSigmas1.push_back(vdwNumber3);
					vdwSigmas2.push_back(vdwNumber4);
					vdwEpsilons.push_back(vdwNumber5);
					vdwIdx.push_back(vdwIdxRow);
					vdwIdxRow.clear();
					std::getline(top, line);
				}
			}
			if (line == "<box>")
			{
				std::getline(top, line);
				while (line != "<end>")
				{
					double boxNumber1 = std::stod(line.substr(0, line.find(delim)));
					double boxNumber2 = std::stod(line.substr(nthSubstr(1, line, delim)));
					double boxNumber3 = std::stod(line.substr(nthSubstr(2, line, delim)));
					box.push_back(boxNumber1);
					box.push_back(boxNumber2);
					box.push_back(boxNumber3);
					std::getline(top, line);
				}
			}

		}
		nAtoms = atomTypes.size();
	}
	catch (...)
	{
		std::cout << "ERROR: fix topology!" << std::endl;
		exit(1);
	}
	top.close();
}

std::ostream& operator<<(std::ostream& os, const Topology& topol)
{
	os << "Topology:";

	os << "\n\tatomtypes: ";
	for (auto& elem: topol.atomTypes)
	{
		os << elem << " ";
	}

	os << "\n\tmasses: ";
	for (auto& elem : topol.masses)
	{
		os << elem << " ";
	}
	os << "\n\tcharges: ";
	for (auto& elem : topol.charges)
	{
		os << elem << " ";
	}

	os << "\n\tbonds:";
	for (int i = 0; i < topol.bondIdx.size(); ++i)
	{
		os << "\n\t";
		for (int j = 0; j < topol.bondIdx[i].size(); ++j)
		{
			os << topol.bondIdx[i][j] << " ";
		}
		os << "L: " << topol.bondLengths[i] << " k: " << topol.bondKs[i];
	}
	os << "\n\tangles:";
	for (int i = 0; i < topol.angleIdx.size(); ++i)
	{
		os << "\n\t";
		for (int j = 0; j < topol.angleIdx[i].size(); ++j)
		{
			os << topol.angleIdx[i][j] << " ";
		}
		// convert back to degrees for printing
		os << "v: " << topol.angleValues[i] * (180.0 / PI) << " θ: " << topol.angleThetas[i];
	}
	os << "\n\tdihedrals:";
	for (int i = 0; i < topol.dihedralIdx.size(); ++i)
	{
		os << "\n\t";
		for (int j = 0; j < topol.dihedralIdx[i].size(); ++j)
		{
			os << topol.dihedralIdx[i][j] << " ";
		}
		// convert back to degrees for printing
		os << "v: " << topol.dihedralValues[i] * (180.0 / PI) << " θ: " << topol.dihedralThetas[i];
	}

	os << "\n\tVDW:";
	for (int i = 0; i < topol.vdwIdx.size(); ++i)
	{
		os << "\n\t";
		for (int j = 0; j < topol.vdwIdx[i].size(); ++j)
		{
			os << topol.vdwIdx[i][j] << " ";
		}
		// convert back to degrees for printing
		os << "σ1: " << topol.vdwSigmas1[i] << " σ2: " << topol.vdwSigmas2[i] << " ε: " << topol.vdwEpsilons[i];
	}
	os << "\n\tbox dimensions:";
	os << "\n\t";
	os << topol.box[0] <<  " " << topol.box[1] << " " << topol.box[2];

	os << "\n";

	return os;
}


Coordinates::Coordinates(std::string crdFile)
{
	readCoordinates(crdFile);
}

Coordinates::~Coordinates()
{
	;
}

void Coordinates::readCoordinates(std::string crdFile)
{
	xyz.clear();  // prepare to read in new coordinates
	std::ifstream crd(crdFile);
	if (!crd.good()) {
		std::cout << "can't open coordinates file" << std::endl;
		exit(1);
	}
	std::string line;
	std::string delim = " ";
	std::vector<double> xyzRow;

	std::getline(crd, line); // skip the number of atoms line
	std::getline(crd, line); // skip the obligatory comment line

	try {
		while(std::getline(crd, line))
		{
			double x = std::stod(line.substr(nthSubstr(1, line, delim)));
			double y = std::stod(line.substr(nthSubstr(2, line, delim)));
			double z = std::stod(line.substr(nthSubstr(3, line, delim)));
			xyzRow.push_back(x);
			xyzRow.push_back(y);
			xyzRow.push_back(z);
			xyz.push_back(xyzRow);
			xyzRow.clear();
		}
	}
	catch(...)
	{
		std::cout << "ERROR: fix coordinate xyz file" << std::endl;
		exit(1);
	}
	crd.close();
}

std::ostream& operator<<(std::ostream& os, const Coordinates& crd)
{
	os << "Coordinates (xyz):";
	os << "\n";
	for (int i = 0; i < crd.xyz.size(); ++i)
	{
		os << i << "> ";
		for (int j = 0; j < 3; ++j)
		{
			os << crd.xyz[i][j] << " ";
		}
		os << "\n";
		
	}
	return os;
}


void Coordinates::checkCoordinates(const Topology& top)
{
	if (top.nAtoms != xyz.size())
	{
		std::cout << "ERROR: number of atoms in coordinate file doesn't match topology" << std::endl;
		exit(1);
	}
}

Velocities::Velocities(const Topology& top, const Parameters& prm)
{
	utils::growFillZeros(xyz, top.nAtoms, 3);
	// initialize xyz to temperature in parameters
}

Velocities::~Velocities()
{
	;
}

std::ostream& operator<<(std::ostream& os, const Velocities& vel)
{
	os << "Velocities (xyz):";
	os << "\n";
	for (int i = 0; i < vel.xyz.size(); ++i)
	{
		os << i << "> ";
		for (int j = 0; j < 3; ++j)
		{
			os << vel.xyz[i][j] << " ";
		}
		os << "\n";
		
	}
	return os;
}

void Velocities::setToTemperature(Topology& top, double temperature)
{
	
	double randNormal;
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, 1.0);

	for (int i = 0; i < top.nAtoms; ++i)
	{
		double m = top.masses[i];
		double prefactor = sqrt((boltzmannK * temperature) / (2 * m));
		randNormal = distribution(generator);
		xyz[i][0] = prefactor * randNormal;
		randNormal = distribution(generator);
		xyz[i][1] = prefactor * randNormal;
		randNormal = distribution(generator);
		xyz[i][2] = prefactor * randNormal;
	}
}

Parameters::Parameters(std::string paramFile)
{
	Parameters::readParameters(paramFile);
}

Parameters::~Parameters()
{
	;
}

void Parameters::readParameters(std::string paramFile)
{
	std::ifstream prm(paramFile);
	if (!prm.good())
	{
		std::cout << "can't open input file" << std::endl;
		exit(1);
	}
	std::string line, delim = " ";

	while (getline(prm, line))
	{
		if (line.rfind("nSteps:") == 0)  // if line starts with nSteps:
		{
			nSteps = std::stoi(line.substr(nthSubstr(1, line, delim)));
		}
		else if (line.rfind("integrator:") == 0)
		{
			integrator = line.substr(nthSubstr(1, line, delim));
			integrator = integrator.substr(1);  // remove extra space in the beginning
		}
		else if (line.rfind("dt:") == 0)
		{
			dt = std::stod(line.substr(nthSubstr(1, line, delim)));
		}
		else if (line.rfind("constrainBonds:") == 0)
		{
			constrainBonds = bool(std::stoi(line.substr(nthSubstr(1, line, delim))));
		}
		else if (line.rfind("temperature:") == 0)
		{
			temperature = std::stod(line.substr(nthSubstr(1, line, delim)));
		}
		else if (line.rfind("reportInterval:") == 0)
		{
			reportInterval = std::stod(line.substr(nthSubstr(1, line, delim)));
		}
	}
	prm.close();

}

void Parameters::checkParameters()
{
	if (nSteps < 1)
	{
		std::cout << "ERROR: nSteps < 1" << std::endl;
		exit(1);
	}

	if (dt <= 0)
	{
		std::cout << "ERROR: dt <= 0" << std::endl;
		exit(1);
	}

	if (integrator != "min" && integrator != "vv" && integrator != "langevin")
	{
		std::cout << "ERROR: integrator must be 'min', 'vv', or 'langevin'" << std::endl;
	}
}

std::ostream& operator<<(std::ostream& os, const Parameters& prm)
{
	os << "Parameters:" << std::endl;
	os << "\tnSteps: " << prm.nSteps << std::endl;
	os << "\tdt: " << prm.dt << std::endl;
	os << "\tintegrator: " << prm.integrator << std::endl;
	os << "\tconstrainBonds: " << prm.constrainBonds << std::endl;
	os << "\ttemperature: " << prm.temperature << std::endl;
	os << "\treportInterval" << prm.reportInterval << std::endl;
	return os;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m)
{
	for (int i = 0; i < m.size(); ++i)
	{
		for (int j = 0; j < m[i].size(); ++j)
		{
			std::cout << m[i][j] << "\t";
		}
		std::cout << std::endl;
		
	}

	return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<double>& v)
{
	for (int i = 0; i < v.size(); ++i)
	{	
		std::cout << v[i] << "\t";
	}
	std::cout << std::endl;

	return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& v)
{
	for (int i = 0; i < v.size(); ++i)
	{	
		std::cout << v[i] << "\n";
	}
	std::cout << std::endl;

	return os;
}


std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<int> >& m)
{
	for (int i = 0; i < m.size(); ++i)
	{
		for (int j = 0; j < m[i].size(); ++j)
		{
			std::cout << m[i][j] << "\t";
		}
		std::cout << std::endl;
		
	}

	return os;
}

namespace utils
{

std::vector<double> pbcAdjust(const Topology& top, std::vector<double>& vec) // only for 3-vectors!
{
	std::vector<double> newVec(3);
	for (int k = 0; k < 3; ++k)
	{
		if (vec[k] > top.box[k] / 2)
		{
			newVec[k] = vec[k] - top.box[k];
		}
		else if (vec[k] < -top.box[k] / 2)
		{
			newVec[k] = vec[k] + top.box[k];
		}
		else
		{
			newVec[k] = vec[k];
		}
	}

	return newVec;
}

double norm(const std::vector<double>& vec) // only for 3-vectors!
{
	double length = 0.0;
	for (int k = 0; k < 3; ++k)
	{
		length += vec[k] * vec[k];
	}

	return sqrt(length);
}


void fillForcesCoordinates(Matrix& forces, int nAtoms)
{
	std::vector<double> force = {0.0, 0.0, 0.0};
	for (int i = 0; i < nAtoms; ++i)
	{
		forces.push_back(force);
	}
}

void addForces(Matrix& forces, const Matrix&& forces2)
{

	for (int i = 0; i < forces.size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			forces[i][j] += forces2[i][j];
		}
	}
}

Matrix subVecs(const Matrix& a, const Matrix& b)
{
	Matrix c;
	fillForcesCoordinates(c, a.size());
	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			c[i][j] = a[i][j] - b[i][j];
	
		}
	}

	return c;
}

Matrix matMul(const Matrix& a, const Matrix& b)
{
	Matrix c;
	fillForcesCoordinates(c, 3);  // 3 x 3 matrix of zeros
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < b.size(); ++k)
			{
				c[i][j] += a[i][k] * b[k][j];
			}
	
		}
	}

	return c;
}

Matrix transpose(const Matrix& a)
{
	Matrix b(3, std::vector<double>(a.size()));  // 3 x 3 matrix
	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			b[j][i] = a[i][j];
		}
	}

	return b;
}

double matNorm(const Matrix& a)
{
	double norm; // www.mathworld.wolfram.com/FrobeniusNorm.html
	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			norm += a[i][j] * a[i][j];
		}
	}
	return std::sqrt(norm);
}

void growFillZeros(Matrix &m, int rows, int cols)
{
	for (int i = 0; i < rows; ++i)
	{
		std::vector<double> row(cols);
		for (int j = 0; j < cols; ++j)
		{
			row[j] = 0.0;
		}
		m.push_back(row);
	}
}

void growFillEye(Matrix &m, int N)
{
	growFillZeros(m, N, N);
	for (int i = 0; i < N; ++i)
	{
		m[i][i] = 1.0;
	}
}

void fillZeros(Matrix &m)
{
	for (int i = 0; i < m.size(); ++i)
	{
		for (int j = 0; j < m[i].size(); ++j)
		{
			m[i][j] = 0.0;
		}
	}
}

Matrix attachEyeToMatrix(Matrix& m, Matrix& id)
{
	/* attaches an identity matrix to the right of an existing matrix */
	Matrix n;
	for (int i = 0; i < m.size(); ++i)
	{	
		std::vector<double> row;
		for (int j = 0; j < m[0].size(); ++j)
		{
			row.push_back(m[i][j]);
		}
		for (int j = 0; j < id.size(); ++j)
		{
			row.push_back(id[i][j]);
		}
		n.push_back(row);
	}

	return n;
}

Matrix detachEyeFromMatrix(Matrix& n, int nEye)
{
	// detaches the left half of an augmented matrix (for the end of a matrix inversion) 
	Matrix clean;

	for (int i = 0; i < n.size(); ++i)
	{	
		std::vector<double> row;
		for (int j = nEye; j < n[0].size(); ++j)  // the inversion result is on the right half
		{
			row.push_back(n[i][j]);
		}
		clean.push_back(row);
	}

	return clean;
}

void switchRows(Matrix& m, int row1, int row2)
{
	int nCols = m[0].size();

	std::vector<double> tmp(nCols);
	for (int i = 0; i < nCols; ++i)
	{
		tmp[i] = m[row1][i]; // read row1 into a tmp
	}
	for (int i = 0; i < nCols; ++i)
	{
		m[row1][i] = m[row2][i]; // replace row1 with row2
	}
	for (int i = 0; i < nCols; ++i)
	{
		m[row2][i] = tmp[i]; 
	}
}

Matrix invertMatrix(Matrix &m)
{
	/*
	 * matrix inversion by row reduction of an augmented matrix (Gaussian elimination)
	 */

	Matrix id;  // the right half of the augemented matrix is initially an identity matrix
	growFillEye(id, m.size());  // identity matrix
	Matrix augm = attachEyeToMatrix(m, id);  // now the original matrix is augmented with an identity matrix
	int nRows = augm.size(); 
	int nCols = augm[0].size(); // nCols

	int i, j, k, l;
	std::vector<double> tempRow1(nCols), tempRow2(nCols);

	/* reduce to reduced row echelon form */

	// sweep forwards
	for (i = 0; i < nRows; ++i)
	{
		j = i;
		// find the first nonzero coefficient after j = i
		while (j < nRows - 1)
		{
			if (augm[i][j] != 0)
			{
				break;
			}
			else
			{
				j += 1;
			}
		}

		if (i != j) switchRows(augm, i, j);

		for (k = i + 1; k < nRows; ++k)
		{
			if (augm[k][i] != 0.0)
			{
				for (l = 0; l < nCols; ++l)
				{
					tempRow1[l] = augm[i][l] * augm[k][i];
					tempRow2[l] = augm[k][l] * augm[i][i];
				}

				if (tempRow1[i] == tempRow2[i])
				{
					for (l = 0; l < nCols; ++l)
					{
						augm[k][l] = tempRow2[l] - tempRow1[l];
					}
				}
				else if (tempRow1[i] == -tempRow2[i])
				{
					for (l = 0; l < nCols; ++l)
					{
						augm[k][l] = tempRow2[l] + tempRow1[l];
					}	
				}
			}
		}

	}

	// sweep backwards
	for (i = nRows - 1; i >= 0; --i)
	{
		j = i;
		while (j > 0)
		{
			if (augm[i][j] != 0.0)
			{
				break;
			}
			else
			{
				j -= 1;
			}
		}

		if (i != j) switchRows(augm, i, j);

		for (k = i - 1; k >= 0; --k)
		{
			if (augm[k][i] != 0)
			{

				for (l = 0; l < nCols; ++l)
				{
					tempRow1[l] = augm[i][l] * augm[k][i];
					tempRow2[l] = augm[k][l] * augm[i][i];
				}

				if (tempRow1[i] == tempRow2[i])
				{
					for (l = 0; l < nCols; ++l)
					{
						augm[k][l] = tempRow2[l] - tempRow1[l];
					}
				}
				else if (tempRow1[i] == -tempRow2[i])
				{
					for (l = 0; l < nCols; ++l)
					{
						augm[k][l] = tempRow2[l] + tempRow1[l];
					}	
				}
			}
		}
	}

	for (i = 0; i < nRows; ++i)
	{
		// reduce to smallest terms
		for (j = nRows; j < nCols; ++j)
		{
			augm[i][j] /= augm[i][i];
		}
	}

	// clean up the left hand side into an identity matrix (for OCD reasons only)
	for (i = 0; i < nRows; ++i)
	{
		augm[i][i] /= augm[i][i];
	}

	Matrix result = detachEyeFromMatrix(augm, id.size()); // now the augmented matrix is broken to reveal the inverted matrix in the right half

	for (int i = 0; i < nRows; ++i)
	{
		for (int j = 0; j < nRows; ++j)
		{
			if (isnan(result[i][j]))
			{
				 throw std::runtime_error("utils.cpp invertMatrix() - matrix not invertible");
			}
		}
	}
	return result;
}

double dot(const std::vector<double>& v1, const std::vector<double>& v2)
{
	double dotProduct = 0.0;
	for (int k = 0; k < 3; ++k)
	{
		dotProduct += v1[k] * v2[k];
	}
	return dotProduct;
}

std::vector<double> cross(const std::vector<double>& a, const std::vector<double>& b)
{
	std::vector<double> prod(3);
	prod[0] = a[1] * b[2] - a[2] * b[1];
	prod[1] = a[2] * b[0] - a[0] * b[2];
	prod[2] = a[0] * b[1] - a[1] * b[0];
	return prod;
}
} // namespace utils

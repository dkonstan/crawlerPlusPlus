#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "forces.hpp"
#include "utils.hpp"



Force::~Force()
{
	;
}


std::ostream& operator<<(std::ostream& os, const Force& force)
{
	os << "Force of type " << force.type;

	return os;
}


void Force::setEnergyFunction()
{
	if (type == "bond")
	{
		energyFunction = getBondEnergy;
	}
	else if (type == "angle")
	{
		energyFunction = getAngleEnergy;
	}
	else if (type == "dihedral")
	{
		energyFunction = getDihedralEnergy;
	}
	else if (type == "vdw")
	{
		energyFunction = getVDWEnergy;
	}
	else if (type == "coulomb")
	{
		energyFunction = getCoulombEnergy;
	}
	else
	{
		std::cout << "energy type not recognized" << std::endl;
		exit(1);
	}
}


void Force::setForceFunction()
{
	if (type == "bond")
	{
		forceFunction = getBondForces;
	}
	else if (type == "angle")
	{
		forceFunction = getAngleForces;
	}
	else if (type == "dihedral")
	{
		forceFunction = getDihedralForces;
	}
	else if (type == "vdw")
	{
		forceFunction = getVDWForces;
	}
	else if (type == "coulomb")
	{
		forceFunction = getCoulombForces;
	}
	else
	{
		std::cout << "force type not recognized" << std::endl;
		exit(1);
	}
}


double getBondEnergy(const Topology& top, const Matrix& xyz)
{
	double bondEnergy = 0.0;
	for (int i = 0; i < top.bondIdx.size(); ++i)
	{
		int atom1Idx = top.bondIdx[i][0];
		int atom2Idx = top.bondIdx[i][1];
		double k = top.bondKs[i];
		double bondLength = top.bondLengths[i];
		std::vector<double> r1r2(3);
		for (int j = 0; j < 3; ++j)
		{
			r1r2[j] = xyz[atom1Idx][j] - xyz[atom2Idx][j];
		}
		utils::pbcAdjust(top, r1r2);
		double distance = utils::norm(r1r2);
		double stretch = distance - bondLength;
		bondEnergy += 0.5 * k * stretch * stretch;
	}
	return bondEnergy;
}


Matrix getBondForces(const Topology& top, const Matrix& xyz)
{
	Matrix bondForces;
	utils::growFillZeros(bondForces, top.nAtoms, 3);

	int i, j;
	for (i = 0; i < top.bondIdx.size(); ++i)
	{
		int atom1Idx = top.bondIdx[i][0];
		int atom2Idx = top.bondIdx[i][1];
		double k = top.bondKs[i];
		double bondLength = top.bondLengths[i];

		std::vector<double> r1r2(3);
		for (j = 0; j < 3; ++j)
		{
			r1r2[j] = xyz[atom1Idx][j] - xyz[atom2Idx][j];
		}

		r1r2 = utils::pbcAdjust(top, r1r2);

		double distance = utils::norm(r1r2);

		double stretch = distance - bondLength;
		std::vector<double> unitVector(3);
		for (j = 0; j < 3; ++j)
		{
			unitVector[j] = r1r2[j] / distance;
		}

		for (j = 0; j < 3; ++j)
		{
			bondForces[atom1Idx][j] += -k * stretch * unitVector[j];
			bondForces[atom2Idx][j] += -k * stretch * (-unitVector[j]);
		}
	}

	return bondForces;
}


double getAngleEnergy(const Topology& top, const Matrix& xyz)
{
	double angleEnergy = 0.0;
	int i, j; 
	for (i = 0; i < top.angleIdx.size(); ++i)
	{

		int atom1Idx = top.angleIdx[i][0];
		int atom2Idx = top.angleIdx[i][1];
		int atom3Idx = top.angleIdx[i][2];
		const double kTheta = top.angleThetas[i];
		std::vector<double> atom1(3), atom2(3), atom3(3);
		for (j = 0; j < 3; ++j)
		{
			atom1[j] = xyz[atom1Idx][j];
			atom2[j] = xyz[atom2Idx][j];
			atom3[j] = xyz[atom3Idx][j];
		}

		std::vector<double> vec12(3);
		std::vector<double> vec32(3);
		for (j = 0; j < 3; ++j)
		{
			vec12[j] = atom1[j] - atom2[j];
			vec32[j] = atom3[j] - atom2[j];
		}
		vec12 = utils::pbcAdjust(top, vec12);
		vec32 = utils::pbcAdjust(top, vec32);
		double norm12 = utils::norm(vec12);
		double norm32 = utils::norm(vec32);
		for (j = 0; j < 3; ++j)
		{
			vec12[j] /= norm12;
			vec32[j] /= norm32;
		}

		double dot1232 = utils::dot(vec12, vec32);
		// correct minor numerical errors that break acos()
		if (dot1232 < -1.0)
		{
			dot1232 = -1.0;
		}
		else if (dot1232 > 1.0)
		{
			dot1232 = 1.0;
		}

		const double angle = acos(dot1232);
		const double angle0 = top.angleValues[i];

		angleEnergy += 0.5 * kTheta * (angle - angle0) * (angle - angle0);

	}

	return angleEnergy;
}


Matrix getAngleForces(const Topology& top, const Matrix& xyz)
{
	Matrix angleForces;
	utils::growFillZeros(angleForces, top.nAtoms, 3);
	const double limitAngle = 0.001;  // how close to 0 or pi before doing something special

	int i, j;
	for (i = 0; i < top.angleIdx.size(); ++i)
	{
		int atom1Idx = top.angleIdx[i][0];
		int atom2Idx = top.angleIdx[i][1];
		int atom3Idx = top.angleIdx[i][2];
		const double kTheta = top.angleThetas[i];
		std::vector<double> atom1(3), atom2(3), atom3(3);
		for (j = 0; j < 3; ++j)
		{
			atom1[j] = xyz[atom1Idx][j];
			atom2[j] = xyz[atom2Idx][j];
			atom3[j] = xyz[atom3Idx][j];
		}

		std::vector<double> vec12(3);
		std::vector<double> vec32(3);
		for (j = 0; j < 3; ++j)
		{
			vec12[j] = atom1[j] - atom2[j];
			vec32[j] = atom3[j] - atom2[j];
		}
		vec12 = utils::pbcAdjust(top, vec12);
		vec32 = utils::pbcAdjust(top, vec32);

		double norm12 = utils::norm(vec12);
		double norm32 = utils::norm(vec32);
		for (j = 0; j < 3; ++j)
		{
			vec12[j] /= norm12;
			vec32[j] /= norm32;
		}

		double dot1232 = utils::dot(vec12, vec32);
		// correct minor numerical errors that break acos()
		if (dot1232 < -1.0)
		{
			dot1232 = -1.0;
		}
		else if (dot1232 > 1.0)
		{
			dot1232 = 1.0;
		}

		const double angle = acos(dot1232);
		const double angle0 = top.angleValues[i];
		std::vector<double> forcei(3);
		std::vector<double> forcek(3);
		std::vector<double> rijVec(3);
		std::vector<double> rkjVec(3);
		std::vector<double> atom1minus2(3);
		std::vector<double> atom3minus2(3);
		std::vector<double> force1(3);
		std::vector<double> force2(3);
		std::vector<double> force3(3);
		double rij;
		double rkj;
		double prefactor = (1.0 / sqrt(1 - (cos(angle) * cos(angle)) ) );

		if (fabs(angle) < limitAngle || fabs(PI - angle) < limitAngle)
		{
			// if very close to 0 or pi -> turn off angle forces (singularity!)
			std::fill(forcei.begin(), forcei.end(), 0.0);
			std::fill(forcek.begin(), forcek.end(), 0.0);
		}
		else
		{
			for (j = 0; j < 3; ++j)
			{
				atom1minus2[j] = atom1[j] - atom2[j];
				atom3minus2[j] = atom3[j] - atom2[j];
			}
			rijVec = utils::pbcAdjust(top, atom1minus2);
			rkjVec = utils::pbcAdjust(top, atom3minus2);
			rij = utils::norm(rijVec);
			rkj = utils::norm(rkjVec);
			for (j = 0; j < 3; ++j)
			{
				forcei[j] = prefactor * (1 / rij) * ((rijVec[j] / rij) * cos(angle) - (rkjVec[j] / rkj));
				forcek[j] = prefactor * (1 / rkj) * ((rkjVec[j] / rkj) * cos(angle) - (rijVec[j] / rij));
			}

		}

		for (j = 0; j < 3; ++j)
		{
			force1[j] = -kTheta * (angle - angle0) * forcei[j];
			force3[j] = -kTheta * (angle - angle0) * forcek[j];
			force2[j] = -(force1[j] + force3[j]);  // sum of forces on the 3 atoms must be zero
		}

		for (j = 0; j < 3; ++j)
		{
			angleForces[atom1Idx][j] += force1[j];
			angleForces[atom2Idx][j] += force2[j];
			angleForces[atom3Idx][j] += force3[j];  // sum of forces on the 3 atoms must be zero
		}

	}
	return angleForces;
}

double getDihedralEnergy(const Topology& top, const Matrix& xyz)
{
	double dihedralEnergy = 0.0;
	int i, j;

	std::vector<double> atom1(3), atom2(3), atom3(3), atom4(3);
	double atom1Idx, atom2Idx, atom3Idx, atom4Idx, kTheta, dott;
	std::vector<double> normal123(3), normal234(3), vec12(3), vec32(3), vec23(3), vec43(3);
	double norm123, norm234;
	for (i = 0; i < top.dihedralIdx.size(); ++i)
	{
		atom1Idx = top.dihedralIdx[i][0];
		atom2Idx = top.dihedralIdx[i][1];
		atom3Idx = top.dihedralIdx[i][2];
		atom4Idx = top.dihedralIdx[i][3];
		for (j = 0; j < 3; ++j)
		{
			atom1[j] = xyz[atom1Idx][j];
			atom2[j] = xyz[atom2Idx][j];
			atom3[j] = xyz[atom3Idx][j];
			atom4[j] = xyz[atom4Idx][j];
		}
		kTheta = top.dihedralThetas[i];

		for (j = 0; j < 3; ++j)
		{
			vec12[j] = atom1[j] - atom2[j];
			vec32[j] = atom3[j] - atom2[j];
			vec23[j] = -vec32[j];
			vec43[j] = atom4[j] - atom3[j];
		}
		vec12 = utils::pbcAdjust(top, vec12);
		vec32 = utils::pbcAdjust(top, vec32);
		vec23 = utils::pbcAdjust(top, vec23);
		vec43 = utils::pbcAdjust(top, vec43);
		normal123 = utils::cross(vec12, vec32);
		normal234 = utils::cross(vec43, vec23);

		norm123 = utils::norm(normal123);
		norm234 = utils::norm(normal234);
		for (j = 0; j < 3; ++j)
		{
			normal123[j] /= norm123;
			normal234[j] /= norm234;
		}

		dott = utils::dot(normal123, normal234);
		if (dott < -1.0) dott = -1.0;
		else if (dott > 1.0) dott = 1.0;

		double dihedralAngle = acos(dott);
		double dihedralAngle0 = top.dihedralValues[i];
		dihedralEnergy += 0.5 * kTheta * (dihedralAngle - dihedralAngle0) * (dihedralAngle - dihedralAngle0);
	}
	return 0.0;
}

Matrix getDihedralForces(const Topology& top, const Matrix& xyz)
{
	Matrix dihedralForces;
	utils::growFillZeros(dihedralForces, top.nAtoms, 3);

	int i, j;

	std::vector<double> atom1(3), atom2(3), atom3(3), atom4(3);
	double atom1Idx, atom2Idx, atom3Idx, atom4Idx, kTheta, dott;
	std::vector<double> normal123(3), normal234(3), vec12(3), vec32(3), vec23(3), vec43(3);
	std::vector<double> force1(3), force2(3), force3(3), force4(3);
	double norm123, norm234, dotprod;
	double dihedralAngle0 = top.dihedralValues[i];
	double dihedralAngle;
	std::vector<double> rijVec(3), rkjVec(3), rklVec(3), rmjVec(3), rnkVec(3);
	std::vector<double> forcei(3), forcej(3), forcek(3), forcel(3);
	double rij, rkj, rmj, rnk;
	for (i = 0; i < top.dihedralIdx.size(); ++i)
	{
		atom1Idx = top.dihedralIdx[i][0];
		atom2Idx = top.dihedralIdx[i][1];
		atom3Idx = top.dihedralIdx[i][2];
		atom4Idx = top.dihedralIdx[i][3];
		for (j = 0; j < 3; ++j)
		{
			atom1[j] = xyz[atom1Idx][j];
			atom2[j] = xyz[atom2Idx][j];
			atom3[j] = xyz[atom3Idx][j];
			atom4[j] = xyz[atom4Idx][j];
		}
		kTheta = top.dihedralThetas[i];

		for (j = 0; j < 3; ++j)
		{
			vec12[j] = atom1[j] - atom2[j];
			vec32[j] = atom3[j] - atom2[j];
			vec23[j] = -vec32[j];
			vec43[j] = atom4[j] - atom3[j];
		}
		vec12 = utils::pbcAdjust(top, vec12);
		vec32 = utils::pbcAdjust(top, vec32);
		vec23 = utils::pbcAdjust(top, vec23);
		vec43 = utils::pbcAdjust(top, vec43);
		normal123 = utils::cross(vec12, vec32);
		normal234 = utils::cross(vec43, vec23);
		norm123 = utils::norm(normal123);
		norm234 = utils::norm(normal234);

		for (j = 0; j < 3; ++j)
		{
			normal123[j] /= norm123;
			normal234[j] /= norm234;
		}

		dott = utils::dot(normal123, normal234);
		if (dott < -1.0) dott = -1.0;
		else if (dott > 1.0) dott = 1.0;

		dihedralAngle = acos(dott);

		if (norm123 < 1e-100 || norm234 < 1e-100)  // TODO is this small enough? if all atoms in dihedral are aligned (rare but possible)
		{
			// std::cout << "frog" << std::endl;
			std::fill(force1.begin(), force1.end(), 0.0);
			std::fill(force2.begin(), force2.end(), 0.0);
			std::fill(force3.begin(), force3.end(), 0.0);
			std::fill(force4.begin(), force4.end(), 0.0);
		}
		else
		{
			rijVec = vec12;  // actually vec21, so need to multiply by -1

			// for (j = 0; j < 3; ++j)
			// {
			// 	rijVec[j] = -rijVec[j];
			// }
			rkjVec = vec32;
			rklVec = vec43;
			rmjVec = utils::cross(rijVec, rkjVec);
			rnkVec = utils::cross(rkjVec, rklVec);

			rij = utils::norm(rijVec);
			rkj = utils::norm(rkjVec);
			rmj = utils::norm(rmjVec);
			rnk = utils::norm(rnkVec);

			double sign = utils::dot(rkjVec, utils::cross(utils::cross(rijVec, rkjVec), utils::cross(rkjVec, rklVec)));
			if (sign > 0.0)
			{
				sign = 1.0;
			}
			else
			{
				sign = -1.0;
			}
			for (j = 0; j < 3; ++j)
			{
				forcei[j] = (rkj / (rmj * rmj)) * rmjVec[j];
				forcel[j] = -(rkj / (rnk * rnk)) * rnkVec[j];
				forcej[j] = ( (utils::dot(rijVec, rkjVec) / (rkj * rkj)) - 1.0 ) * forcei[j] - (utils::dot(rklVec, rkjVec) / (rkj * rkj)) * forcel[j];
				forcek[j] = ( (utils::dot(rklVec, rkjVec) / (rkj * rkj)) - 1.0 ) * forcel[j] - (utils::dot(rijVec, rkjVec) / (rkj * rkj)) * forcei[j];
			}

			for (j = 0; j < 3; ++j)
			{
				force1[j] = sign * -kTheta * (dihedralAngle - dihedralAngle0) * forcei[j];
				force2[j] = sign * -kTheta * (dihedralAngle - dihedralAngle0) * forcej[j];
				force3[j] = sign * -kTheta * (dihedralAngle - dihedralAngle0) * forcek[j];
				force4[j] = sign * -kTheta * (dihedralAngle - dihedralAngle0) * forcel[j];
			}
		}
		for (j = 0; j < 3; ++j)
		{
			dihedralForces[atom1Idx][j] += force1[j];
			dihedralForces[atom2Idx][j] += force2[j];
			dihedralForces[atom3Idx][j] += force3[j];
			dihedralForces[atom4Idx][j] += force4[j];
		}

  	}
	return dihedralForces;
}



bool bonded(const Topology& top, int i, int j)
{
	bool areBonded = false;
	for (int n = 0; n < top.bondIdx.size(); ++n)
	{
		if (top.bondIdx[n][0] == i && top.bondIdx[n][1] == j)
		{
			areBonded = true;
			break;
		}
		else if (top.bondIdx[n][1] == i && top.bondIdx[n][0] == j)
		{
			areBonded = true;
			break;
		}
	}

	return areBonded;
}

bool nonBondedExclude(const Topology& top, int i, int j)
{
	// are the two atoms a pair that is excluded from nonbonded interactions?
	// e.g. the two hydrogens on a water molecule
	bool sameMolecule = false;
	for (int n = 0; n < top.nonBondedExclusionIdx.size(); ++n)
	{
		if (top.nonBondedExclusionIdx[n][0] == i && top.nonBondedExclusionIdx[n][1] == j)
		{
			sameMolecule = true;
			break;
		}
		else if (top.nonBondedExclusionIdx[n][1] == i && top.nonBondedExclusionIdx[n][0] == j)
		{
			sameMolecule = true;
			break;
		}
	}

	return sameMolecule;
}


double getVDWEnergy(const Topology& top, const Matrix& xyz)
{
	double vdwEnergy = 0.0, r;
	double sigma1, sigma2, eps1, eps2;
	int i, j, k;
	std::vector<double> r1r2(3);
	for (i = 0; i < top.nAtoms; ++i)
	{
		sigma1 = top.vdwSigmas[i];
		eps1 = top.vdwEpsilons[i];
		for (j = i + 1; j < top.nAtoms; ++j)
		{
			if (nonBondedExclude(top, i, j))
			{
				continue;
			}
			else
			{
				sigma2 = top.vdwSigmas[j];
				eps2 = top.vdwEpsilons[j];

				double sig = (sigma1 + sigma2) / 2;
				double eps = sqrt(eps1 * eps2);
				for (k = 0; k < 3; ++k) r1r2[k] = xyz[i][k] - xyz[j][k];
				r1r2 = utils::pbcAdjust(top, r1r2);
				r = utils::norm(r1r2);
				vdwEnergy += 4 * eps * (pow(sig / r, 12) - pow(sig / r, 6));	
			}
		}
	}
	return vdwEnergy;
}


Matrix getVDWForces(const Topology& top, const Matrix& xyz)
{
	Matrix vdwForces;
	utils::growFillZeros(vdwForces, top.nAtoms, 3);
 	double r;
	double sigma1, sigma2, eps1, eps2;
	int i, j, k;
	std::vector<double> r1r2(3);
	std::vector<double> r1r2Unit(3);

	for (i = 0; i < top.nAtoms; ++i)
	{
		sigma1 = top.vdwSigmas[i];
		eps1 = top.vdwEpsilons[i];
		for (j = i + 1; j < top.nAtoms; ++j)
		{
			if (nonBondedExclude(top, i, j))
			{
				continue;
			}
			else
			{
				sigma2 = top.vdwSigmas[j];
				eps2 = top.vdwEpsilons[j];

				double sig = (sigma1 + sigma2) / 2;
				double eps = sqrt(eps1 * eps2);
				for (k = 0; k < 3; ++k) r1r2[k] = xyz[i][k] - xyz[j][k];
				r1r2 = utils::pbcAdjust(top, r1r2);
				r = utils::norm(r1r2);

				double r1r2SqNorm = 0.0;
				for (k = 0; k < 3; ++k) r1r2SqNorm += r1r2[k] * r1r2[k];

				vdwForces[i][0] += (24 * r1r2[0] * eps * pow(sig, 6) * (-pow(r1r2SqNorm, 3) + 2 * pow(sig, 6))) / pow(r1r2SqNorm, 7);
				vdwForces[i][1] += (24 * r1r2[1] * eps * pow(sig, 6) * (-pow(r1r2SqNorm, 3) + 2 * pow(sig, 6))) / pow(r1r2SqNorm, 7);
				vdwForces[i][2] += (24 * r1r2[2] * eps * pow(sig, 6) * (-pow(r1r2SqNorm, 3) + 2 * pow(sig, 6))) / pow(r1r2SqNorm, 7);

				vdwForces[j][0] -= (24 * r1r2[0] * eps * pow(sig, 6) * (-pow(r1r2SqNorm, 3) + 2 * pow(sig, 6))) / pow(r1r2SqNorm, 7);
				vdwForces[j][1] -= (24 * r1r2[1] * eps * pow(sig, 6) * (-pow(r1r2SqNorm, 3) + 2 * pow(sig, 6))) / pow(r1r2SqNorm, 7);
				vdwForces[j][2] -= (24 * r1r2[2] * eps * pow(sig, 6) * (-pow(r1r2SqNorm, 3) + 2 * pow(sig, 6))) / pow(r1r2SqNorm, 7);
			}
		}
	}

	return vdwForces;
}


double getCoulombEnergy(const Topology& top, const Matrix& xyz)
{
	double coulombEnergy = 0.0;
	double r;
	int i, j, k;
	std::vector<double> rij(3);
	for (i = 0; i < top.charges.size(); ++i)
	{
		for (j = i + 1; j < top.charges.size(); ++j)
		{
			if (nonBondedExclude(top, i, j))
			{
				continue;
			}
			for (k = 0; k < 3; ++k)
			{
				rij[k] = xyz[i][k] - xyz[j][k];
			}
			rij = utils::pbcAdjust(top, rij);
			r = utils::norm(rij);

			coulombEnergy += COULOMB_K * (top.charges[i] * top.charges[j]) / r;
		}
	}
	return coulombEnergy;
}


Matrix getCoulombForces(const Topology& top, const Matrix& xyz)
{
	// TODO check this
	Matrix coulombForces;
	utils::growFillZeros(coulombForces, top.nAtoms, 3);
	double r;
	int i, j, k;
	std::vector<double> rij(3);
	for (i = 0; i < top.charges.size(); ++i)
	{
		for (j = i + 1; j < top.charges.size(); ++j)
		{
			if (nonBondedExclude(top, i, j))  // TODO this is right, right?
			{
				continue;
			}
			else
			{
				for (k = 0; k < 3; ++k)
				{
					rij[k] = xyz[i][k] - xyz[j][k];
				}
				rij = utils::pbcAdjust(top, rij);
				r = utils::norm(rij);

				// for (k = 0; k < 3; ++k)
				// {
				// 	rij[k] /= r;
				// }

				for (k = 0; k < 3; ++k)
				{
					coulombForces[i][k] += ((COULOMB_K * top.charges[i] * top.charges[j]) * rij[k]) / pow(r, 3);
					coulombForces[j][k] -= ((COULOMB_K * top.charges[i] * top.charges[j]) * rij[k]) / pow(r, 3);
				}			
			}

		}
	}
	return coulombForces;
}



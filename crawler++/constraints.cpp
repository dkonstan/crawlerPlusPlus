#include <iostream>
#include "utils.hpp"
#include "constraints.hpp"


void shakePositions(Topology& top, Matrix &xyz, std::vector<double>& bondLengths, std::vector<double>& lambdas){
	std::vector<std::vector<int> > idx = top.bondIdx;
	std::vector<double> m = top.masses;
	std::vector<double> d = bondLengths;
	const int N = top.nAtoms;
	const int nC = bondLengths.size();
	std::vector<double> sigma(nC);
	int i, j, k, l;
	std::vector<double> ri(3), rj(3);
	std::vector<std::vector<std::vector<double> > > deltaSigma;
	Matrix deltaSigmaMatrix;
	std::vector<double> adjustment(3);
	const double tolerance = 1e-10;
	std::vector<double> adjustmentMagnitudes(N);
	double convergenceIndicator = 1000.0;  // just a random large value to get started

	utils::growFillZeros3D(deltaSigma, nC, N, 3);
	utils::growFillZeros(deltaSigmaMatrix, nC, nC);

	// evaluate constraints at original position
	for (k = 0; k < nC; ++k)
	{
		i = idx[k][0];
		j = idx[k][1];
		sigma[k] = utils::sqNorm(xyz[i], xyz[j]) - (d[k] * d[k]);
	}

	// calculate each constraint's gradient wrt each atom's position
	for (k = 0; k < nC; ++k)
	{
		i = idx[k][0];
		j = idx[k][1];
		ri = xyz[i];
		rj = xyz[j];
		for (l = 0; l < 3; ++l)
		{	
			deltaSigma[k][i][l] = 2 * (ri[l] - rj[l]);  // gradient of constraint wrt i
			deltaSigma[k][j][l] = 2 * (rj[l] - ri[l]);  // gradient of constraint wrt j
		}
	}
	// calculate initial adjustment based on guesses for Lagrange multipliers
	for (i = 0; i < N; ++i)
	{
		std::fill(adjustment.begin(), adjustment.end(), 0.0);
		for (k = 0; k < nC; ++k)
		{
			for (l = 0; l < 3; ++l)
			{
				adjustment[l] += (1 / m[i]) * (lambdas[k] * deltaSigma[k][i][l]);
			}
		}
		// make adjustment to positions based on guess for lambda
		for (l = 0; l < 3; ++l)
		{
			xyz[i][l] += adjustment[l];
		}
	}

	// now time to calculate the remaining correction to the position to within the tolerance
	while (convergenceIndicator > tolerance)
	{
		// calculate constraints
		for (k = 0; k < nC; ++k)
		{
			i = idx[k][0];
			j = idx[k][1];
			sigma[k] = utils::sqNorm(xyz[i], xyz[j]) - (d[k] * d[k]);
		}

		for (i = 0; i < deltaSigma.size(); ++i) utils::fillZeros(deltaSigma[i]);
		utils::fillZeros(deltaSigmaMatrix);

		// calculate constraint gradients
		for (k = 0; k < nC; ++k)
		{
			i = idx[k][0];
			j = idx[k][1];
			ri = xyz[i];
			rj = xyz[j];
			for (l = 0; l < 3; ++l)
			{
				deltaSigma[k][i][l] = 2 * (ri[l] - rj[l]);  // gradient of constraint wrt i
				deltaSigma[k][j][l] = 2 * (rj[l] - ri[l]);  // gradient of constraint wrt j
			}
		}

		// build matrix for solving linear equations
		for (l = 0; l < nC; ++l)
		{
			for (k = 0; k < nC; ++k)
			{
				for (i = 0; i < N; ++i)
				{
					deltaSigmaMatrix[l][k] += (1 / m[i]) * utils::dot(deltaSigma[l][i], deltaSigma[k][i]);
				}
		
			}
		}

		// solve for remaining adjustment
		std::vector<double> dlambdas = utils::solveEquation(deltaSigmaMatrix, utils::flipSign(sigma));
		// std::cout << deltaSigmaM << std::endl;
		// std::cout << deltaSigmaMatrix << std::endl;
		// update lambdas to get a better guess next time
		for (k = 0; k < nC; ++k) lambdas[k] += dlambdas[k];

		for (i = 0; i < N; ++i)
		{
			std::fill(adjustment.begin(), adjustment.end(), 0.0);
			for (k = 0; k < nC; ++k)
			{
				for (l = 0; l < 3; ++l) adjustment[l] += (1 / m[i]) * (dlambdas[k] * deltaSigma[k][i][l]);
			}
			for (l = 0; l < 3; ++l)
			{
				xyz[i][l] += adjustment[l];	
			}
			adjustmentMagnitudes[i] = utils::norm(adjustment);
		}
		convergenceIndicator = utils::max(adjustmentMagnitudes);
	}
}


void shakeVelocities(Topology& top, Matrix& xyz, Matrix& vel, std::vector<double>& bondLengths, std::vector<double>& lambdas, double dt, Matrix& newForces)
{
	std::vector<std::vector<int> > idx = top.bondIdx;
	std::vector<double> m = top.masses;
	std::vector<double> d = bondLengths;
	const int N = top.nAtoms;
	const int nC = bondLengths.size();
	std::vector<double> sigma(nC);
	int i, j, k, l;
	std::vector<double> ri(3), rj(3);
	Matrix deltaSigmaMatrix;
	std::vector<double> adjustment(3), adjustment2(3);
	std::vector<std::vector<std::vector<double> > > deltaSigma;
	std::vector<double> deltaSigmaVector(nC);
	std::vector<double> mus;

	utils::growFillZeros3D(deltaSigma, nC, N, 3);
	utils::growFillZeros(deltaSigmaMatrix, nC, nC);

	// calculate constraint gradients
	for (k = 0; k < nC; ++k)
	{
		i = idx[k][0];
		j = idx[k][1];
		ri = xyz[i];
		rj = xyz[j];
		for (l = 0; l < 3; ++l)
		{
			deltaSigma[k][i][l] = 2 * (ri[l] - rj[l]);  // gradient of constraint wrt i
			deltaSigma[k][j][l] = 2 * (rj[l] - ri[l]);  // gradient of constraint wrt j
		}
	}

	for (i = 0; i < N; ++i)
	{
		std::fill(adjustment.begin(), adjustment.end(), 0.0);
		for (k = 0; k < nC; ++k)
		{
			for (l = 0; l < 3; ++l)
			{
				adjustment[l] += (1 / m[i]) * (1 / dt) * lambdas[k] * deltaSigma[k][i][l];
			}
		}

		for (l = 0; l < 3; ++l)
		{
			adjustment[l] += ((dt / 2) * (1 / m[i]) * newForces[i][l]);
		}

		for (l = 0; l < 3; ++l)
		{
			vel[i][l] += adjustment[l];
		}
	}

	std::fill(deltaSigmaVector.begin(), deltaSigmaVector.end(), 0.0);
	for (k = 0; k < nC; ++k)
	{
		for (i = 0; i < N; ++i)
		{
			deltaSigmaVector[k] += utils::dot(deltaSigma[k][i], vel[i]);
		}
	}

	utils::fillZeros(deltaSigmaMatrix);
	for (k = 0; k < nC; ++k)
	{
		for (l = 0; l < nC; ++l)
		{
			for (i = 0; i < N; ++i)
			{
				std::vector<double> val(N);
				for (int n = 0; n < N; ++n)
				{
					val[n] = (1 / m[i]) * deltaSigma[l][i][n];
				}
				deltaSigmaMatrix[k][l] += utils::dot(deltaSigma[k][i], val);
			}
		}
	}
	
	mus = utils::solveEquation(deltaSigmaMatrix, utils::flipSign(deltaSigmaVector));

	for (i = 0; i < N; ++i)
	{
		std::fill(adjustment2.begin(), adjustment2.end(), 0.0);
		for (k = 0; k < nC; ++k)
		{
			for (l = 0; l < 3; ++l)
			{
				adjustment2[l] += (1 / m[i]) * (mus[k] * deltaSigma[k][i][l]);  // Lagrange adjustment
			}
		}
		for (l = 0; l < 3; ++l)
		{
			vel[i][l] += adjustment2[l];
		}
	}

}

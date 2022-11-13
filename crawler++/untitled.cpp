// void Simulation::runLD(std::ofstream& logFile, std::ofstream &trajFile)
// {
// 	nDegreesOfFreedom = 3 * top.nAtoms - top.bondIdx.size() - top.angleIdx.size();  // 3N - # constraints
// 	int i, j;
// 	double mi;
// 	double temp = param.temperature;
// 	double gamma = param.collisionFreq;
// 	double dt = param.dt;
// 	// Matrix At, xit, thetat;
// 	// utils::growFillZeros(At, top.nAtoms, 3);
// 	// utils::growFillZeros(xit, top.nAtoms, 3);
// 	// utils::growFillZeros(thetat, top.nAtoms, 3);
// 	Matrix vPr;;
// 	Matrix deltaV;
// 	utils::growFillZeros(vPr, top.nAtoms, 3);
// 	utils::growFillZeros(deltaV, top.nAtoms, 3);


// 	double alpha = 1 - exp(-gamma * dt);
// 	double sqrtConstant;

// 	for (int s = 0; s < param.nSteps; ++s)
// 	{
// 		// std::cout << randNormal(0.0, 1.0) << std::endl;
// 		if (s % param.reportInterval == 0) report(logFile, trajFile, s);

// 		currEnergy = getTotalEnergy(crd.xyz);
// 		currKineticEnergy = getKineticEnergy(vel.xyz);
// 		currTemperature = (2.0 * currKineticEnergy) / (nDegreesOfFreedom * boltzmannK);
// 		currForces = getTotalForces(crd.xyz);

// 		for (i = 0; i < top.nAtoms; ++i)
// 		{
// 			mi = top.masses[i];
// 			for (j = 0; j < 3; ++j)
// 			{
// 				vPr[i][j] = vel.xyz[i][j] + (1 / mi) * currForces[i][j] * dt;
// 			}
// 		}

// 		for (i = 0; i < top.nAtoms; ++i)
// 		{
// 			mi = top.masses[i];
// 			sqrtConstant = sqrt(((boltzmannK * temp) / mi) * alpha * (2 - alpha));
// 			for (j = 0; j < 3; ++j)
// 			{
// 				deltaV[i][j] = -alpha * vPr[i][j] + sqrtConstant * randNormal(0.0, 1.0);
// 			}
// 		}

// 		for (i = 0; i < top.nAtoms; ++i)
// 		{
// 			for (j = 0; j < 3; ++j)
// 			{
// 				crd.xyz[i][j] += (vPr[i][j] + 0.5 * deltaV[i][j]) * dt;
// 			}
// 		}

// 		for (i = 0; i < top.nAtoms; ++i)
// 		{
// 			for (j = 0; j < 3; ++j)
// 			{
// 				vel.xyz[i][j] = vPr[i][j] + deltaV[i][j];
// 			}
// 		}
// 		// std::cout << deltaV << std::endl;
// 		// std::cout << vPr << std::endl;
// 		// std::cout << crd << std::endl;
// 		// std::cout << vel << std::endl;
// 		// exit(0);
// 		std::vector<double> cog(3, 0.0);
// 		for (i = 0; i < top.nAtoms; ++i)
// 		{
// 			cog[0] += crd.xyz[i][0];
// 			cog[1] += crd.xyz[i][1];
// 			cog[2] += crd.xyz[i][2];
// 		}

// 		cog[0] /= top.nAtoms;
// 		cog[1] /= top.nAtoms;
// 		cog[2] /= top.nAtoms;
// 		for (i = 0; i < top.nAtoms; ++i)
// 		{
// 			crd.xyz[i][0] -= cog[0];
// 			crd.xyz[i][1] -= cog[1];
// 			crd.xyz[i][2] -= cog[2];
// 		}
// 	}
// }
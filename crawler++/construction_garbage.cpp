	// for (i = 0; i < top.vdwIdx.size(); ++i)
	// {
	// 	double atom1Idx, atom2Idx, sigma, epsilon, r;
	// 	std::vector<double> r1r2(3);

	// 	atom1Idx = top.vdwIdx[i][0];
	// 	atom2Idx = top.vdwIdx[i][1];
	// 	sigma = top.vdwSigmas1[i]; // fix sigmas1 sigmas2
	// 	epsilon = top.vdwEpsilons[i];

	// 	for (j = 0; j < 3; ++j)
	// 	{
	// 		r1r2[j] = xyz[atom1Idx][j] - xyz[atom2Idx][j];
	// 	}
	// 	r1r2 = utils::pbcAdjust(top, r1r2);
	// 	r = utils::norm(r1r2);

	// 	vdwEnergy += 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
	// }


	// for (i = 0; i < top.vdwIdx.size(); ++i)
	// {
	// 	double atom1Idx, atom2Idx, sigma, epsilon, r;
	// 	std::vector<double> r1r2(3), unitVector(3);

	// 	atom1Idx = top.vdwIdx[i][0];
	// 	atom2Idx = top.vdwIdx[i][1];
	// 	sigma = top.vdwSigmas1[i]; // fix sigmas1 sigmas2
	// 	epsilon = top.vdwEpsilons[i];
	// 	for (j = 0; j < 3; ++j)
	// 	{
	// 		r1r2[j] = xyz[atom1Idx][j] - xyz[atom2Idx][j];
	// 	}
	// 	r1r2 = utils::pbcAdjust(top, r1r2);
	// 	r = utils::norm(r1r2);

	// 	for (j = 0; j < 3; ++j)
	// 	{
	// 		unitVector[j] = r1r2[j] / r;	
	// 	}
	// 	for (j = 0; j < 3; ++j)
	// 	{
	// 		vdwForces[atom1Idx][j] += -4 * epsilon * (
	// 							(6 * pow(sigma, 6)) / pow(r, 7) - (12 * pow(sigma, 12)) / pow(r, 13)
	// 						) * unitVector[j];
	// 		vdwForces[atom2Idx][j] += -4 * epsilon * (
	// 							(6 * pow(sigma, 6)) / pow(r, 7) - (12 * pow(sigma, 12)) / pow(r, 13)
	// 						) * (-unitVector[j]);	
	// 	}
	// }
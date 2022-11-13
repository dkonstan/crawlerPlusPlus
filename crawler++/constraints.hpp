#include <iostream>
#include "utils.hpp"

#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

// SHAKE constraints for bond lengths
void shakePositions(Topology& top, Matrix &xyz, std::vector<double>& bondLengths, std::vector<double>& lambdas);
void shakeVelocities(Topology& top, Matrix& xyz, Matrix& vel, std::vector<double>& bondLengths, std::vector<double>& lambdas, double dt, Matrix& newForces);

// SHAKE constraints for bond angles
void shakeAngles(Topology& top, Matrix &xyz, std::vector<double>& bondAngles, std::vector<double>& lambdas);
void shakeAnglesVelocities(Topology& top, Matrix& xyz, Matrix& vel, std::vector<double>& bondAngles, std::vector<double>& lambdas, double dt, Matrix& newForces);

void getAngleGradients(Topology& top, std::vector<double>& ri, std::vector<double>& rj, std::vector<double>& rk,
					   double angle, std::vector<double>& angleGradient1, std::vector<double>& angleGradient2,
					   std::vector<double>& angleGradient3);
#endif // CONSTRAINTS.HPP
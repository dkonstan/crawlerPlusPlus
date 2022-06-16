#include <iostream>
#include "utils.hpp"

#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

void shakePositions(Topology& top, Matrix &xyz, std::vector<double>& bondLengths, std::vector<double>& lambdas);
void shakeVelocities(Topology& top, Matrix& xyz, Matrix& vel, std::vector<double>& bondLengths, std::vector<double>& lambdas, double dt, Matrix& newForces);


#endif // CONSTRAINTS.HPP
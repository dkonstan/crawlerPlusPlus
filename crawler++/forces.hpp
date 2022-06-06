#include <string>
#include <functional>
#include "utils.hpp"


#ifndef CRAWLER_FORCES
#define CRAWLER_FORCES


double getBondEnergy(const Topology& top, const Matrix& xyz);
double getAngleEnergy(const Topology& top, const Matrix& xyz);
double getDihedralEnergy(const Topology& top, const Matrix& xyz);
double getVDWEnergy(const Topology& top, const Matrix& xyz);
double getCoulombEnergy(const Topology& top, const Matrix& xyz);

Matrix getBondForces(const Topology& top, const Matrix& xyz);
Matrix getAngleForces(const Topology& top, const Matrix& xyz);
Matrix getDihedralForces(const Topology& top, const Matrix& xyz);
Matrix getVDWForces(const Topology& top, const Matrix& xyz);
Matrix getCoulombForces(const Topology& top, const Matrix& xyz);

class Force
{
public:
	Force(std::string type): type(type)
	{
		setEnergyFunction();
		setForceFunction();
	};
	~Force();
	std::function<double(const Topology&, const Matrix&)> energyFunction;
	std::function<Matrix(const Topology&, const Matrix&)> forceFunction;
	std::string type;
private:
	void setEnergyFunction();
	void setForceFunction();
friend std::ostream& operator<<(std::ostream& os, const Force& force);
};

#endif // CRAWLER_FORCES
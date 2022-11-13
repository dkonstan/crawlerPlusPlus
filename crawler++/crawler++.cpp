#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include "utils.hpp"
#include "forces.hpp"
#include "simulation.hpp"


int main(int argc, char* argv[])
{
	const std::string helpMessage = \
"\n\nMISSING ARGUMENTS: input must be\n\
crawler++ -t [topology file] -c [coordinates file] -i [parameter file] -l [log file] -o [trajectory file] -x [outfile]\n\
or crawler won't crawl\n\n";

	if (argc < 13)
	{
		std::cout << helpMessage << std::endl;
		return 1;
	}
	else 
	{
		std::vector<std::string> args;
		std::string topologyFilename = "";
		std::string coordFilename = "";
		std::string inputFilename = "";
		std::string trajectoryFilename = "";
		std::string logFilename = "";
		std::string outFilename = "";
		std::vector<std::string> inputErr;
		int i;
		for (i = 1; i < argc; ++i) args.push_back(argv[i]);

		i = 0;
		while (i < argc)
		{
			if (args[i] == "-t")
			{
				if (args[i + 1][0] != '-') topologyFilename = args[i + 1]; else inputErr.push_back("argument parser problem, check input around topology");
			}
			else if (args[i] == "-i")
			{
				if (args[i + 1][0] != '-') inputFilename = args[i + 1]; else inputErr.push_back("argument parser problem, check input around input file");
			}
			else if (args[i] == "-c")
			{
				if (args[i + 1][0] != '-') coordFilename = args[i + 1]; else inputErr.push_back("argument parser problem, check input around coordinate file");
			}
			else if (args[i] == "-o")
			{
				if (args[i + 1][0] != '-') trajectoryFilename = args[i + 1]; else inputErr.push_back("argument parser problem, check input around trajectory file");
			}
			else if (args[i] == "-l")
			{
				if (args[i + 1][0] != '-') logFilename = args[i + 1]; else inputErr.push_back("argument parser problem, check input around log file");
			}
			else if (args[i] == "-x")
			{
				if (args[i + 1][0] != '-') outFilename = args[i + 1]; else inputErr.push_back("argument parser problem, check input around outfile");
			}
			++i;
		}
		// if (topologyFilename.length() == 0) inputErr.push_back("missing topology (-t)");
		// if (coordFilename.length() == 0) inputErr.push_back("missing coordinates (-c)");
		// if (inputFilename.length() == 0) inputErr.push_back("missing input file (-i)");
		// if (trajectoryFilename.length() == 0) inputErr.push_back("missing trajectory (-o)");
		// if (logFilename.length() == 0) inputErr.push_back("missing log file (-l)");
		// if (outFilename.length() == 0) inputErr.push_back("missing outfile file (-x)");
		// std::cout << inputErr.size() << std::endl;
		// if (inputErr.size() > 0)
		// {
		// 	std::cout << helpMessage << std::endl;
		// 	return 1;
		// }

		std::ofstream logfile;
		logfile.open(logFilename);

		Topology top(topologyFilename);
		top.checkTopology();

		logfile << top << std::endl;
		Coordinates crd(coordFilename);
		crd.checkCoordinates(top);
		logfile << crd << std::endl;
		Parameters param(inputFilename);
		param.checkParameters();
		logfile << param << std::endl;
		Velocities vel(top, param);
		srand((unsigned) time(NULL));  // unique seed for the random number generator
		vel.setToTemperature(top, param.temperature);

		std::vector<Force> forces;
		if (top.bondIdx.size() > 0) forces.push_back(Force("bond"));
		if (top.angleIdx.size() > 0) forces.push_back(Force("angle"));
		// if (top.dihedralIdx.size() > 0) forces.push_back(Force("dihedral"));
		if (top.vdwSigmas.size() > 0) forces.push_back(Force("vdw"));
		if (top.charges.size() > 0) forces.push_back(Force("coulomb"));

		Simulation sim(top, crd, vel, param, forces, trajectoryFilename, logFilename, outFilename);
		logfile << sim << std::endl;
		logfile.close();

		sim.crawl();
	}

	return 0;
}

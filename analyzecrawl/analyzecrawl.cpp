#include <iostream>
#include <fstream>
#include <string>
#include <vector>


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

double mean(std::vector<double> data)
{
	double sum = 0.0;
	for (int i = 0; i < data.size(); ++i)
	{
		sum += data[i];
	}
	return sum / data.size();
}

int main(int argc, char* argv[])
{
	std::string helpMessage = "\n./analyzecrawl [whatToExtract] [logfile] > [.csv]\n";

	if (argc < 3)
	{
		std::cout << helpMessage << std::endl;
	}

	std::string arg = argv[1];
	std::string log = argv[2];
	std::ifstream logfile(log);
	std::string line, prefix;
	const std::string delim = ": ";
	int position;
	if (arg == "potentialEnergy") position = 0;
	else if (arg == "kineticEnergy") position = 1;
	else if (arg == "totalEnergy") position = 2;
	else if (arg == "temperature") position = 3;
	else
	{
		std::cout << helpMessage << std::endl;
		exit(1);
	}
	std::vector<double> data;
	double value;
	if (logfile.good())
	{
		try {
			while (getline(logfile, line)) 
			{
				if (line.rfind("step") == 0)
				{

					for (int j = 0; j < position; ++j)
					{
						getline(logfile, line);
					}
					while (getline(logfile, line))
					{
						prefix = line.substr(0, line.find(delim));
						value = std::stod(line.substr(prefix.length() + delim.length()));
						if (line == "---- CRAWLER++ has crawled.---") exit(0);
						std::cout << value << std::endl;
						data.push_back(value);
						for (int j = 0; j < 5; ++j) getline(logfile, line);
						if (!logfile.good()){
							std::cout << "mean: " << mean(data) << std::endl;
							exit(0);
						}
					}
				}
			}
		}
		catch(std::out_of_range& e) { std::cout << "mean: " << mean(data) << std::endl; exit(0); } ;
	}
	return 0;
}
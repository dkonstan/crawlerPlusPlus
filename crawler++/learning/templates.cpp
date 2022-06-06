#include <iostream>
#include <vector>
#include "templates.hpp"

int main()
{
	double a = 0, b = 1;
	std::cout << starwars::max(a, b) << std::endl;

	std::vector<double> v(10);
	for (auto i = v.begin(); i != v.end(); ++i)
	{
		*i = 1.0;
	}

	// v.push_back(2.0);
	// for (auto& elem : v)
	// {
	// 	std::cout << elem << std::endl;
	// }

	// auto it = std::find(v.begin(), v.end(), 4);
	// std::cout << "After " << *it << " comes " << *(it + 1) << std::endl;
	// auto it2 = v.insert(it + 1, 5);
	// std::vector<matrix> v2;
	// std::cout << v2 << std::endl;
	// std::vector<std::vector<double>> V;
	std::fill(v.begin(), v.end(), 3.0);
	// V.push_back(v);
	// V.push_back(v);========

	for (auto& elem : v)
	{
		std::cout << elem << std::endl;
	}
	// for (auto& elem : V)========
	// {
	// 	std::cout << elem << std::endl;
	// }
	return 0;
}


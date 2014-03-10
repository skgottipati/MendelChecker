#ifndef CL_H
#define CL_H
#include "cart.h"
#include "optparse/OptionParser.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <unordered_map>
#include <ctime>

using namespace std;

string itoa(int x);



class founders_mean_GLs{

	public:
	string snpid;
	std::map<std::string,long double> GLs;
	void set_snpid(string );
	void set_GLs(std::map<std::string,long double> GL);
};


void new_compute_likelihood(vector< vector <vector <LINE>>> snps, vector< vector< LINE>> founders, unordered_map<std::string, double> Penetrance, string filename, double alpha, string unfFLAG);

#endif
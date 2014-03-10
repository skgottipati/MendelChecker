#ifndef CART_H
#define CART_H
#include "GenotypeInfo.h"

#include "genotypeLikelihood.h"

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

typedef std::vector<GLPROB> Vi;
typedef std::vector<Vi> Vvi;

typedef std::map<std::string,long double> Mgl;

struct Digits {
	Mgl::const_iterator begin;
	Mgl::const_iterator end;
	Mgl::const_iterator me;
	bool founder;
};


struct offd {
	std::map<std::string,long double> GLs;
	string xlg;
};

typedef std::vector<Digits> Vd;


long double cart_product(Vi& in, std::unordered_map<std::string, double> Penetrance, string chrom);

#endif
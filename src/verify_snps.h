#ifndef VP_H
#define VP_H

#include "split.h"
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

void verify_nuclear_family(std::vector<std::string> fam);

void verify_pedigrees(std::map<std::string, std::string> pedigree);
//void verify_snps(vector< vector <vector <LINE>>> snps, vector< vector< LINE>> founders);


#endif
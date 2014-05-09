#ifndef FR_H
#define FR_H

#include "split.h"
#include "genoped.h"
#include "GenotypeInfo.h"
#include "genotypeLikelihood.h"
#include "cart.h"

#include "computeLikelihood.h"

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

std::string getFileName(const std::string strPath);

void fileread(string fname, int bufsize, string phredFLAG, double alpha, string unfFLAG);

#endif
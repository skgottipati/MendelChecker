#ifndef GL_H
#define GL_H

#include "GenotypeInfo.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <sstream>


using namespace std;


std::map<std::string,long double> phred2prob(LINE snp);

class GLPROB {

	public:
	std::map<std::string,long double> GLs;
	void setelem(LINE&, map<string, long double>&, string);
	string chromosome;
	int position;
	string snpid;
	string familyid;
	string individualid;
	string mom;
	string dad;
	string xlg;	
};

//void GLPROB::setelem(LINE& line, map<string, long double>& fm, string likelihoodFLAG);

#endif
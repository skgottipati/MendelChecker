
#ifndef GI_H
#define GI_H
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <sstream>

using namespace std;

class LINE {

	public:
	std::map<std::string,int> GLs;	
	void setelem(string);	
	string chromosome;
	int position;
	string snpid;
	string familyid;
	string individualid;
	string mom;
	string dad;
	string sex;	
};

//void LINE::setelem(string line);

#endif
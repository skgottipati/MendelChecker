#ifndef __GENOPED_H
#define __GENOPED_H

// functions to split a string by a specific delimiter
#include "split.h"
#include "GenotypeInfo.h"
#include "computeLikelihood.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <iterator>
#include <iomanip>
#include <unordered_map>
#include <sstream>
#include <algorithm>
using namespace std;

vector<string> splitdelim(string line);

std::map<std::string, std::string> pedigree_reader(std::string pedfilename);

vector<string> get_individuals_from_pedigreemap(const std::map<std::string, std::string>* pedigree);

std::unordered_map<std::string, int> generateGTmap(vector<string> alleles);

std::unordered_map<std::string, double> computePenetrance();

std::string get_FileName(const std::string strPath);

void read_geno(string genofield, int bufsize, std::string phredFLAG, std::string vcfname, const std::map<std::string, std::string>* pedigree, std::unordered_map<std::string, double> Penetrance, double alpha, string unfFLAG);

//void parse(string& line, bool parseSamples);

#endif
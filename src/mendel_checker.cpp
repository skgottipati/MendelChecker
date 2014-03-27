#include "fileread.h"
#include "genoped.h"


#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include <cmath>
#include <iterator>
#include <iomanip>
#include <ctime>
#include <unordered_map>


#include "GenotypeInfo.h"
#include "genotypeLikelihood.h"
#include "cart.h"
#include "computeLikelihood.h"



#include "../optparse/OptionParser.h"

//#include <boost/math/distributions/chi_squared.hpp>
//#include <boost/cstdint.hpp>






using namespace std;

using namespace optparse;


string getelem(std::istringstream linestream){
	string w;
	if (getline (linestream, w, '\t'))
		return w;
	else
		return "";
}


int main(int argc, char* argv[]) {
	if (argc < 2) 
	{
		cerr << "Usage: input file name is required\n";
		return EXIT_FAILURE;
	}

#ifndef DISABLE_USAGE
  const string usage = "usage: %prog [OPTION]... DIR [FILE]...";
#else
  const string usage = SUPPRESS_USAGE;
#endif
  const string version = "%prog 1.0\nCopyright (C) 2014 Srikanth Gottipati\n"
    "License GPLv3+: GNU GPL version 3 or later "
    "<http://gnu.org/licenses/gpl.html>.\n"
    "This is free software: you are free to change and redistribute it.\n"
    "There is NO WARRANTY, to the extent permitted by law.";
  const string desc = "Mendel checker";
  const string epilog = "This program is great....please use it!!!";

  OptionParser parser = OptionParser()
    .usage(usage)
    .version(version)
    .description(desc)
    .epilog(epilog)	
#ifdef DISABLE_INTERSPERSED_ARGS
    .disable_interspersed_args()
#endif
  ;
	parser.add_option("-f", "--genoped") .dest("filename") .type("string") .help("input geno-ped file name with path") .metavar("FILE");
	parser.add_option("-g", "--vcf") .dest("vcf") .type("string") .help("input vcf file name with path") .metavar("FILE");
	parser.add_option("-e", "--ped") .dest("ped") .type("string") .help("input ped file name with path") .metavar("FILE");
	parser.add_option("-n", "--snpsperloop") .dest("snpsperloop") .type("int") .set_default(10000) .help("number of snps compute per loop");
	parser.add_option("-d", "--genofield") .dest("genofield") .type("string") .help("VCF genotype field, options: PL, GL, GP") .set_default("PL");
	parser.add_option("-p", "--sexPrior") .action("store") .dest("sexPrior") .type("double") .set_default(0.05) .help("default: %default sexPrior");
	parser.add_option("-u", "--uniform") .dest("uniformFLAG") .type("string").help("default: %default (population), true (uniform)") .set_default("false") .metavar("STRING");
	
	Values& options = parser.parse_args(argc, argv);
	vector<string> args = parser.args();
	clock_t sec;
	sec = clock();
	cout << "argc:" << argc << " : " << options["alpha"] << ":" << options["unif"] << endl;
	string vcf = (string) options.get("vcf");
	if (vcf != ""){
		std::unordered_map<std::string, double> Penetrance = computePenetrance();
		std::map<std::string, std::string> pedigree = pedigree_reader((string) options.get("ped"));
		read_geno((string) options.get("genofield"), (int) options.get("snpsperloop"),(string) options.get("vcf"), &pedigree, Penetrance, (double) options.get("sexPrior"), (string) options.get("uniformFLAG"));
	}
	else
		fileread(options["filename"], (double) options.get("sexPrior"), (string) options.get("uniformFLAG"));
	cout << options["alpha"] << ":" << sizeof(double) << endl;
	sec = clock() - sec;
	cout << "Compute time : "  << "\t" << ((long double)sec/CLOCKS_PER_SEC) << " seconds" <<  endl;

}


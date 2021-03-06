#include "genoped.h"
#include "verify_snps.h"
#include "fileread.h"

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
	//if (argc < 2) 
	//{
	//	cerr << "\nA vcf & pedigree file or a genoped file are required.\nEnter './MendelChecker --help' for available options.\nFor more information, please visit http://code.google.com/p/mendelchecker/wiki/Documentation\n\n";
	//	exit (EXIT_FAILURE);
	//}

#ifndef DISABLE_USAGE
  const string usage = "\nMendelChecker -f genoPedFile [-o outputFile] [-g genotypeProbabilityFormat ] [-a sexPrior] [-u useUniformPrior] \nMendelChecker -v vcfFile -p pedFile [-o outputFile] [-g genotypeProbabilityFormat ] -a sexPrior [-u useUniformPrior]";
#else
  const string usage = SUPPRESS_USAGE;
#endif
  const string version = "%prog 1.0\nCopyright (C) 2014 Srikanth Gottipati\n"
    "License GPLv3+: GNU GPL version 3 or later "
    "<http://gnu.org/licenses/gpl.html>.\n"
    "This is free software: you are free to change and redistribute it.\n"
    "There is NO WARRANTY, to the extent permitted by law.";
  const string desc = "MendelChecker: A C++ software for quality control in SNP discovery from next generation sequencing data using Mendelian inheritance in pedigrees";
  const string epilog = "For more information, please visit https://sourceforge.net/p/mendelchecker/wiki/Documentation/\n";

//Usage:
//MendelChecker -f genoPedFile [-o outputFile] [-g genotypeProbabilityFormat ] -a sexPrior [-u useUniformPrior]
//MendelChecker -v vcfFile -p pedFile [-o outputFile] [-g genotypeProbabilityFormat ] -a sexPrior [-u useUniformPrior]
// 
//MendelChecker: A C++ software for quality control in SNP discovery from next-generation sequencing data
//using Mendelian inheritance in pedigrees
// 
//Options:
//  -h, --help        	        	        	show this help message and exit
//  --version                     	        	show program's version number and exit
//  -f FILE, --genoped=FILE        	        	input geno-ped file name with path
//  -v FILE, --vcf=FILE           	        	input vcf file name with path
//  -p FILE, --ped=FILE           	        	input ped file name with path
//  -o FILE, --out=FILE           	        	output file name with path
//  -g STRING, --genoProb=STRING        	genotype probability format, options: PL, GL, GP (default: PL)
//  -a DOUBLE, --sexPrior=DOUBLE        	prior probability of sex-linkage (default: 0.05)
//  -u STRING, --uniform=STRING        	use uniform prior instead of expected population genotype frequencies (default: false)
// 
//For more information, please visit https://sourceforge.net/p/mendelchecker/wiki/Documentation/
  
  
  
  
  OptionParser parser = OptionParser()
    .usage(usage)
    .version(version)
    .description(desc)
    .epilog(epilog)
#ifdef DISABLE_INTERSPERSED_ARGS
    .disable_interspersed_args()
#endif
  ;
	parser.add_option("-f", "--genoped") .dest("filename") .type("string") .set_default("") .set_default("") .help("input geno-ped file name with path") .metavar("FILE");
	parser.add_option("-v", "--vcf") .dest("vcf") .type("string") .set_default("") .help("input vcf file name with path") .metavar("FILE");
	parser.add_option("-p", "--ped") .dest("ped") .type("string") .set_default("") .help("input ped file name with path") .metavar("FILE");
	parser.add_option("-o", "--out") .dest("outdir") .type("string") .set_default("") .help("output file name with path") .metavar("FILE");
	parser.add_option("-m", "--memoryAlloc") .dest("memoryAlloc") .type("string") .set_default("1GB") .help("size of usable memory (default: %default) by the program");
	parser.add_option("-g", "--genoProb") .dest("genofield") .type("string") .set_default("PL") .help("genotype probability format, options: PL, GL, GP (default: %default)");
	parser.add_option("-s", "--phredScore") .dest("phredScore") .type("string") .set_default("true") .help("genotype quality as phred score (default: %default)");
	parser.add_option("-a", "--sexPrior") .action("store") .dest("sexPrior") .type("double") .set_default(0.05) .help("prior probability of sex-linkage (default: %default)");
	parser.add_option("-u", "--uniform") .dest("uniformFLAG") .type("string").help("use uniform prior instead of expected population genotype frequencies (default: %default)") .set_default("false") .metavar("STRING");
	Values& options = parser.parse_args(argc, argv);
	if (argc <2)
	{
		const char* argvs[] = {"MendelChecker" , "--help"};
		Values& options2 = parser.parse_args(2, argvs);
		cout <<  (string) options2.get("phredScore") << endl;
		
	}
	//else{
	//	Values& options = parser.parse_args(argc, argv);
	//}
	vector<string> args = parser.args();
	clock_t sec;
	sec = clock();
	//cout << "argc:" << argc << " : " << options["alpha"] << ":" << options["unif"] << endl;
	string vcf = (string) options.get("vcf");
	string pedfilename = (string) options.get("ped");
	string genofilename = (string) options.get("filename");
	string outfname = (string) options.get("outdir");
	double alp = (double) options.get("sexPrior");
	string GF = (string) options.get("genofield");
	string uniFLAG = (string) options.get("uniformFLAG");
	string memAlloc = (string) options.get("memoryAlloc");
	string phreds = (string) options.get("phredScore");
	
	std::size_t foundat = memAlloc.find("GB");
	//cout << foundat << endl;
	if (foundat == 0)
	{
		cout << "No memory has been allocated." << endl;
		exit (EXIT_FAILURE);
	}
	if (foundat == std::string::npos)
	{
		cout << "Error in memory allocation input. Input values can only be in GB units, ex: 8GB" << endl;
		exit (EXIT_FAILURE);
	}
	string memAllocSubstr = memAlloc.substr(0,foundat);
	int bufsize = atoi(memAllocSubstr.c_str());
	//cout << "buf: " << bufsize << endl;
	cout << "\n" << "Inputs:" << endl;	
	if (vcf != "")
	{
		if (pedfilename == "")
		{
			cout << "Pedigree file(.ped) is required." << endl;
			exit (EXIT_FAILURE);
		}
		cout << "VCF file : " << vcf << endl;
		cout << "Pedigree file : " << pedfilename << endl;
	}
	if (vcf == "" && genofilename == "")
	{
		cout << "A vcf file or a genoped file is required." << endl;
		exit (EXIT_FAILURE);
	}
	if (vcf != "" && genofilename != "")
	{
		cout << "Both vcf and genoped files are provided, the default mode is to run MendelChecker on vcf+ped." << endl;
	}
	if (vcf == "" && genofilename !="")
	{
		cout << "Genoped file : " << genofilename << endl;
	}
	//cout << "Output file name prefix : " << outfname << endl;
	cout << "Genotype probability : " << GF << endl;
	cout << "Sex-prior alpha : " << alp << endl;
	cout << "Uniform prior : " << uniFLAG << endl;
	cout << "Memory allocated for buffer : " << memAlloc << endl;
	cout << "Phred score based genotype quality : " << phreds << "\n" << endl;
	cout << "Outputs:" << endl;
	string fname = "";
	if (outfname != "")
		fname = outfname;
	else if (vcf !="")
		fname = vcf;
	else if (vcf=="" && genofilename !="" )
		fname = genofilename;
		
	cout << "Pedigree likelihoods : " << fname << ".pedigreelikelihoods" << endl;
	cout << "SNP scores : " << fname << ".snpScores\n" << endl;
	if (alp*(1-alp)<=0)
	{
		cout << "Sex prior should be set to a value between 0 and 1." << endl;
		exit (EXIT_FAILURE);
	}
	std::unordered_map<std::string, double> Penetrance = computePenetrance();

	
	if (vcf != ""){
		//cout << "vcf" << endl;
		std::map<std::string, std::string> pedigree = pedigree_reader((string) options.get("ped"));
		verify_pedigrees(pedigree);		
		read_geno((string) options.get("genofield"), bufsize, phreds, vcf, &pedigree, Penetrance, (double) options.get("sexPrior"), (string) options.get("uniformFLAG"), fname);
	}
	else if (vcf == "" && genofilename != "")
	{
		//cout << "geno" << endl;
		fileread(genofilename, bufsize, phreds, (double) options.get("sexPrior"), (string) options.get("uniformFLAG"), fname);
		//cout << "end" << endl;
	}

	sec = clock() - sec;
	cout << "Compute time : "  << "\t" << ((long double)sec/CLOCKS_PER_SEC) << " seconds" <<  endl;
	
	return 0;

}


#include "genotypeLikelihood.h"


std::map<std::string,long double> phred2prob(LINE snp){
	long double total = 0.0;
	long double base = 10.0;
	long double div = 10.0;
	std::map<std::string, long double> GL;

	for (auto pos=snp.GLs.begin(); pos != snp.GLs.end(); pos++)
	{
		total += std::pow( base , (long double) -(pos->second)/10 );
//		cout << total << ":" << std::pow( base , (long double) -(pos->second)/10 ) << "\t";
	}
//	cout << endl;
	for (auto pos=snp.GLs.begin(); pos != snp.GLs.end(); pos++)
	{
		GL.insert(std::make_pair(pos->first, std::pow( base , -((long double) (pos->second))/div )/total));
//		cout << pos->first << ":" << std::pow( base , -((long double) (pos->second))/div )/total << "\t";
	}
//	cout << "total:" << total << endl;
	return GL;
}

void GLPROB::setelem(LINE& line, map<string, long double>& fm, string likelihoodFLAG){
	this->chromosome = line.chromosome;
	this->position = line.position;
	this->snpid = line.snpid;
	this->familyid = line.snpid;
	this->individualid = line.individualid;
	this->mom = line.mom;
	this->dad = line.dad;
	this->xlg = line.sex;
	// if (line.sex == 1)
		// this->xlg = "HOM";
	// else if (line.sex == 2)
		// this->xlg = "HET";
	// else if (line.sex == 0)
		// this->xlg = "UNK";
		
		
	if (likelihoodFLAG == "INF")
		if (((line.GLs.begin()))->second == -1)
			this->GLs = fm;
		else		
			this->GLs = phred2prob(line);
	else if (likelihoodFLAG == "UNINF")
		this->GLs = fm;
//	for (auto pos=this->GLs.begin(); pos!=this->GLs.end(); pos++){
//		cout << pos->first << ":" << pos->second << "\t";
//	}
//	cout << endl;	

}
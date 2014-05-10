#include "genotypeLikelihood.h"


std::map<std::string,long double> phred2prob(LINE snp, string phredFLAG){
	long double total = 0.0;
	long double base = 10.0;
	long double div = 10.0;
	long double eps = 0.0001;
	std::map<std::string, long double> GL;

	for (auto pos=snp.GLs.begin(); pos != snp.GLs.end(); pos++)
	{
		if (phredFLAG == "true")
		{
			total += std::pow( base , (long double) -(pos->second)/10 );
		}
		else
		{
			total += pos->second;
		}
//		cout << total << ":" << std::pow( base , (long double) -(pos->second)/10 ) << "\t";
	}
//	cout << endl;


	if (phredFLAG == "false")
	{
		try
		{
			if (total > 1+eps || total < 1-eps)
			{
				cout << "Sum of genotype probabilities for individual " << snp.individualid << " on chromosome " << snp.chromosome << " at loci " << snp.position << " is " << total << endl;
				throw "Sum of genotype probabilities does not add to 1.";
			}		
		}
		catch(const char* Message)
		{
			cerr << "Error: " << Message << endl;
			exit (EXIT_FAILURE);
		}
	}
	
	
	for (auto pos=snp.GLs.begin(); pos != snp.GLs.end(); pos++)
	{
		if (phredFLAG == "true")
		{
			GL.insert(std::make_pair(pos->first, std::pow( base , ((long double) -(pos->second))/div )/total));
		}
		else
		{
			GL.insert(std::make_pair(pos->first, (long double) pos->second/total));
		}
//		cout << pos->first << ":" << std::pow( base , -((long double) (pos->second))/div )/total << "\t";
	}
//	cout << "total:" << total << endl;

	return GL;
}

void GLPROB::setelem(LINE& line, map<string, long double>& fm, string likelihoodFLAG, string phredFLAG){
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
		
		
	try
	{
		for (auto pos=this->GLs.begin(); pos!=this->GLs.end(); pos++){
			if (pos->second != -1 && pos->second <0)
				throw "Genotype probability is less than 1.";
		}		
	}
	catch(const char* Message)
	{
		cerr << "Error: " << Message << endl;
		exit (EXIT_FAILURE);
	}		
		

	if (likelihoodFLAG == "INF")
		if (((line.GLs.begin()))->second == -1)
			this->GLs = fm;
		else		
			this->GLs = phred2prob(line, phredFLAG);
	else if (likelihoodFLAG == "UNINF")
		this->GLs = fm;

//	for (auto pos=this->GLs.begin(); pos!=this->GLs.end(); pos++){
//		cout << pos->first << ":" << pos->second << "\t";
//	}
//	cout << endl;

}
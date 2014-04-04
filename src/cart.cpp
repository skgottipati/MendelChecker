#include "cart.h"

long double cart_product(Vi& in, std::unordered_map<std::string, double> Penetrance, string chrom){

	Vd founders_GLiter;
	Digits het_founder;
	Digits hom_founder;
	vector<offd> offsprings_GLs;
	vector<offd> nonzeroTrans_offsprings_GL;
	Vd nonzeroTrans_offsprings_GLiter;

// 	Start founder iterators and store offspring GLs in a vector.

	for (Vi::const_iterator it = in.begin(); it != in.end(); ++it)
	{
        Digits d = {(*it).GLs.begin(), (*it).GLs.end(), (*it).GLs.begin(), ((*it).mom=="0" && (*it).dad=="0")};
		offd od = {(*it).GLs, (*it).xlg};
		if ((*it).mom=="0" && (*it).dad=="0" && (*it).xlg == "1")
			hom_founder = d;
		else if ((*it).mom=="0" && (*it).dad=="0" && (*it).xlg == "2")
			het_founder = d;
		else if ((*it).mom=="0" && (*it).dad=="0" && (*it).xlg == "0")
			cout << "Found parent with unknown gametic type" << endl;			
		else
			offsprings_GLs.push_back(od);
	}		
	founders_GLiter.push_back(hom_founder);
	founders_GLiter.push_back(het_founder);

	long double Likelihood = 0.0;
	int fnder = 0;
	while(1)
	{

		long double founder_prod = 1;
		string founder_genotype = "";
		//int run =0;


		for (Vd::const_iterator it = founders_GLiter.begin(); it != founders_GLiter.end(); it++)
		{
			founder_prod *= it->me->second;
			founder_genotype += it->me->first;
		}
		long double prod0 = founder_prod;
		for (auto noit = offsprings_GLs.begin(); noit != offsprings_GLs.end(); ++noit)
		{
			string chrom_offspring_gametic;		
			if (chrom == "A")
				chrom_offspring_gametic = chrom;
			else if (chrom == "X")
				chrom_offspring_gametic = chrom+noit->xlg;	
			long double sum = 0;
			for (auto it = (noit->GLs).begin(); it != (noit->GLs).end(); it++)
			{
				sum += it->second * Penetrance[founder_genotype+it->first+chrom_offspring_gametic];
			}
			prod0 *= sum;
		}
		Likelihood += prod0;

		nonzeroTrans_offsprings_GL.clear();
		nonzeroTrans_offsprings_GLiter.clear();
		for(Vd::iterator it = founders_GLiter.end()-1; ; )
		{
//			cout << "was here" << endl;
		    	// okay, I started at the left instead. sue me
			++(it->me);
//			cout << "was here" << endl;
			if(it->me == it->end)
			{
				if(it == founders_GLiter.begin())
				{
			    		// I'm the last digit, and I'm about to roll
//					cout << founder_genotype << "\t" << Likelihood << endl;
			    		return Likelihood;
				} 
				else 
				{
					// cascade
					it->me = it->begin;
					--it;
					fnder++;
//					cout <<  founder_genotype << " : " << Likelihood <<  endl;
				}
			} 
			else 
			{
				// normal
				fnder++;
//				cout <<  founder_genotype << " : " << Likelihood <<  endl;
				break;
			}

		}
//		cout << "was here too" << endl;
//		cout << founder_genotype << "\t" << Likelihood << endl;
	}



}

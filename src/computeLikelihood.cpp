#include "computeLikelihood.h"
#include "cart.h"

string itoa(int x){
	ostringstream ss;
	ss << x;
//	string str = ss.str();
	return ss.str();
}


void founders_mean_GLs::set_snpid(string snpid){
	this->snpid = snpid;
}

void founders_mean_GLs::set_GLs(std::map<std::string,long double> GL){
	this->GLs = GL;
}


void new_compute_likelihood(vector< vector <vector <LINE>>> snps, vector< vector< LINE>> founders, unordered_map<std::string, double> Penetrance, string filename, double alpha, string unfFLAG, string phredFLAG){

	static int func_call = 0;
	func_call++;

	std::vector<string> genotypes;
	genotypes.reserve(10);
	genotypes.insert(genotypes.end(),  {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"});




//	mean of founder alleles

//	cout << "Computing mean of founder genotype likelihoods" << endl;
//
//	vector<founders_mean_GLs> FMGLs;
//	FMGLs.reserve(founders.size());
//	
//	LINE temp_founder;
//	std::map<std::string,long double> GL;
//
//	std::map<std::string,long double> meanGL;
//	for (int i=0; i<10; i++)
//	{
//		meanGL.insert(std::make_pair(genotypes[i], 0.0));
//	}
//
//
////	ofstream mfglfile ("/home/chaitanya/Desktop/nancy/meanFounderGL.txt", ios::out );
////	if (func_call == 1)
////	{
//////		ofstream mfglfile ("/home/chaitanya/Desktop/nancy/meanFounderGL.txt", ios::out );
////		mfglfile << "SNPID" << "\t" << "#FOUNDERSwData" << "\t" ;
////		for (auto gt = genotypes.begin(); gt != genotypes.end(); gt++)
////		{
////			mfglfile << *gt << "\t";
////		}
////		mfglfile << "SUM_GL" << endl;
////	}
////	else
////	{
////		mfglfile.close();
//		ofstream mfglfile (filename+"_meanFounderGL.txt", ios::out|ios::app  );
////	}
//	for (vector< vector<LINE> >::size_type i = 0; i < founders.size(); i++)
//	{
//		founders_mean_GLs fmgl;
//		int num_founders = 0;
//		long double check_norm = 0;
//		for (auto mpos = meanGL.begin(); mpos != meanGL.end(); mpos++)
//			mpos->second = 0;
//		for (vector< LINE> ::size_type j = 0; j < founders[i].size(); j++)
//		{
//
//			temp_founder = founders[i][j];
//
////			cout << temp_founder.snpid << "\t";
//			auto cpos = temp_founder.GLs.begin();
//			if (cpos->second == -1) continue;
//			else
//			{
//				GL = phred2prob(temp_founder);
//				for (auto pos = GL.begin(); pos != GL.end(); pos++)
//				{
//					meanGL[pos->first]+=pos->second;
//		
//				}
//				num_founders++;
//			}					
//			
//			
//		}
//		if (num_founders !=0)
//		{
//			for (auto mpos = meanGL.begin(); mpos != meanGL.end(); mpos++)
//			{
//				mpos->second = mpos->second/num_founders;
//				check_norm += mpos->second;
//			}
//		}
//		else
//		{
//			for (auto mpos = meanGL.begin(); mpos != meanGL.end(); mpos++)
//			{
//				mpos->second = 0.1;
//				check_norm += mpos->second;
//			}
//
//		}
//		
//		
//		fmgl.set_snpid(temp_founder.snpid);
//		fmgl.set_GLs(meanGL);
//		FMGLs.emplace_back(fmgl);
////		cout << founders[i].size() << "\t" << fmgl.snpid << endl;
//		mfglfile <<  fmgl.snpid << "\t" << num_founders <<  "\t" ;
//		for (auto gt = genotypes.begin(); gt != genotypes.end(); gt++)
//		{
//			if (meanGL.find(*gt) == meanGL.end())
//				mfglfile << "0" << "\t";
//			else
//				mfglfile << meanGL[*gt] << "\t";
//		}
//		mfglfile << check_norm << endl;
//		
//	}
//	mfglfile.close();
//	cout << "Done!!!" << endl;


//	population genotype vector from founder genotypes

//	cout << "Computing population genotype vector from founders" << endl;

	std::vector<string> alleles;
	alleles.reserve(4);
	alleles.insert(alleles.end(),  {"A", "C", "G", "T"});

	vector<founders_mean_GLs> popGLs;
	popGLs.reserve(founders.size());
	
	LINE temp_popfounder;
	std::map<std::string,long double> pGL;

	std::map<std::string,long double> allelePr;
	for (int i=0; i<4; i++)
	{
		allelePr.insert(std::make_pair(alleles[i], 0.0));
	}


	std::map<std::string, long double> popGT;

////	ofstream popGTfile ("/home/chaitanya/Desktop/nancy/populationGL.txt", ios::out );
////	if (func_call == 1)
////	{
//////		ofstream popGTfile ("/home/chaitanya/Desktop/nancy/populationGL.txt", ios::out );
////		popGTfile << "SNPID" << "\t" << "#FOUNDERSwData" << "\t";
////		for (auto gt = genotypes.begin(); gt != genotypes.end(); gt++)
////		{
////			popGTfile << *gt << "\t";
////		}
////		popGTfile <<  "SUM_GL" << endl;
////	}
////	else
////	{	
////		popGTfile.close();
//		ofstream popGTfile (filename+"_populationGL.txt", ios::out|ios::app  );
////	}

	//ofstream popGTfile (filename+"_populationGL.txt", ios::out|ios::app  );
	for (vector< vector<LINE> >::size_type i = 0; i < founders.size(); i++)
	{
		founders_mean_GLs popgl;
		int num_founders = 0;
		long double check_norm = 0;
		for (auto mpos = allelePr.begin(); mpos != allelePr.end(); mpos++)
			mpos->second = 0;
		for (vector< LINE> ::size_type j = 0; j < founders[i].size(); j++)
		{

			temp_popfounder = founders[i][j];

//			cout << temp_founder.snpid << "\t";
			auto cpos = temp_popfounder.GLs.begin();
			if (cpos->second == -1) continue;
			else
			{
				pGL = phred2prob(temp_popfounder, phredFLAG);
				allelePr["A"] += 2*pGL["AA"] + pGL["AC"] + pGL["AG"] + pGL["AT"];
				allelePr["C"] += 2*pGL["CC"] + pGL["AC"] + pGL["CG"] + pGL["CT"];				
				allelePr["G"] += 2*pGL["GG"] + pGL["AG"] + pGL["CG"] + pGL["GT"];
				allelePr["T"] += 2*pGL["TT"] + pGL["AT"] + pGL["CT"] + pGL["GT"];
				num_founders++;
			}		
			
			
		}
		if (num_founders !=0)
		{
			for (auto mpos = allelePr.begin(); mpos != allelePr.end(); mpos++)
			{
				mpos->second = mpos->second/(2*num_founders);
				check_norm += mpos->second;
			}
		}
		else
		{
			for (auto mpos = allelePr.begin(); mpos != allelePr.end(); mpos++)
			{
				mpos->second = 0.1;
				check_norm += mpos->second;
			}
		}
		long double popGTtot = 0;
		for (auto mpos = allelePr.begin(); mpos != allelePr.end(); mpos++)
		{
			mpos->second = mpos->second/check_norm;
			popGTtot += mpos->second;
		}		

		long double prodallelePr;		
		for (auto gtpos = genotypes.begin(); gtpos != genotypes.end(); gtpos++)
		{
			if (gtpos->substr(0,1) == gtpos->substr(1,1))
				prodallelePr = allelePr[gtpos->substr(0,1)] * allelePr[gtpos->substr(1,1)];
			else
				prodallelePr = 2 * allelePr[gtpos->substr(0,1)] * allelePr[gtpos->substr(1,1)];
			if (prodallelePr !=0 && unfFLAG == "false")
			{
				popGT.insert(std::make_pair(*gtpos, prodallelePr));				
			}
			else if (unfFLAG == "true" )
			{
				popGT.insert(std::make_pair(*gtpos, 0.1));
			}
		}

		popgl.set_snpid(temp_popfounder.snpid);
		popgl.set_GLs(popGT);
		popGLs.emplace_back(popgl);
//		cout << founders[i].size() << "\t" << fmgl.snpid << endl;
		//popGTfile <<  popgl.snpid << "\t" << num_founders <<  "\t" ;
//		for (auto mpos = popGT.begin(); mpos != popGT.end(); mpos++)
//		{
//			popGTfile << mpos->first << ":" << mpos->second << ",";
//		}
		//for (auto gt = genotypes.begin(); gt != genotypes.end(); gt++)
		//{
		//	if (popGT.find(*gt) == popGT.end())
		//		popGTfile << "0" << "\t";
		//	else
		//		popGTfile << popGT[*gt] << "\t";
		//}
		//popGTfile <<  popGTtot << endl;
		popGT.clear();
		
	}
	//popGTfile.close();
//	cout << "Done!!!" << endl;



//	computing un-informative likelihoods


//	cout << "Computing Un-informative likelihoods" << endl;

//	ofstream ulfile ("/home/chaitanya/Desktop/nancy/uninformativelikelihoods.txt", ios::out );	
//	if (func_call == 1)
//	{
////		ofstream ulfile ("/home/chaitanya/Desktop/nancy/uninformativelikelihoods.txt", ios::out );	
//		ulfile << "SNPID" << "\t" << "#OFFSPRINGS" << "\t"  << "uninf-LIKELIHOOD" << endl;
//	}
//	else
//	{
//		ulfile.close();
//		ofstream ulfile (filename+"_uninformativelikelihoods.txt", ios::out|ios::app  );
//	}

	//ofstream ulfile (filename+"_uninformativelikelihoods.txt", ios::out|ios::app  );
	vector < std::map<string, long double>> UninfLikelihoods;
//	auto ufm=FMGLs.begin();
	auto ufm=popGLs.begin();
	long double ulikelihood;
	string ulsnpid;
	for (auto snp = snps.begin() ; snp != snps.end(); snp++)
	{
		int run = 0;
		std::map<string, long double> fammap;
		for (auto fam = (*snp).begin() ; fam != (*snp).end(); fam++)
		{
			vector<GLPROB> glpv;
			int famsize = 0;
			int homs = 0;
			int hets = 0;
			int unks = 0;
			for (auto ind = (*fam).begin(); ind!=(*fam).end(); ind++)
			{
				GLPROB glp;
				glp.setelem(*ind, ufm->GLs, "UNINF", phredFLAG);
				glpv.emplace_back(glp);
				famsize++;
				ulsnpid = ind->snpid;
				if (ind->sex == "1")
					homs++;
				else if (ind->sex == "2")
					hets++;
				else if (ind->sex == "0")
					unks++;
//				cout <<  "snpid:" << ind->snpid << " famid:" << ind->familyid << " individ:" <<ind->individualid << endl;
			}

			run++;
//			cout <<  run << ind->snpid << ind->familyid << ind->individualid << endl;
			string famkeyA = "A" + itoa(famsize) + "A";
			string famkeyX = "X" + itoa(famsize) + "X" + itoa(homs) + "X" + itoa(hets) + "X" + itoa(unks) + "X";
			if (fammap.size() != 0)
			{
				if (fammap.count(famkeyA)==0)
				{
					ulikelihood = cart_product(glpv, Penetrance, "A");
					fammap.insert(std::make_pair(famkeyA, ulikelihood));
					//ulfile << left << setw(10) << ulsnpid <<"\t"<< setw(10) << famkeyA << "\t" <<right << setw(15) << setprecision(6) << ulikelihood << endl;
				}
				if (fammap.count(famkeyX)==0)
				{
					ulikelihood = cart_product(glpv, Penetrance, "X");
					fammap.insert(std::make_pair(famkeyX, ulikelihood));
					//ulfile << left << setw(10) << "\t" << ulsnpid << setw(10) << famkeyX << "\t" << right << setw(15) << setprecision(6) << ulikelihood << endl;					
				}
			}
			else
			{
				ulikelihood = cart_product(glpv, Penetrance, "A");
				fammap.insert(std::make_pair(famkeyA, ulikelihood));
				//ulfile << left << setw(10) << ulsnpid <<"\t" <<  setw(10) << famkeyA <<"\t"<< right << setw(15) << setprecision(6)  << ulikelihood << endl;
				ulikelihood = cart_product(glpv, Penetrance, "X");
				fammap.insert(std::make_pair(famkeyX, ulikelihood));
				//ulfile << left << setw(10) << ulsnpid <<"\t"<<  setw(10) << famkeyX << "\t"<<right << setw(15) << setprecision(6)  << ulikelihood << endl; 				
			}
		}
		ufm++;
		UninfLikelihoods.emplace_back(fammap);
	}
	//ulfile.close();
//	cout << "Done!!!" << endl;




//	computing likelihoods  (33*(4**8)+51*(2**8)+16)


//	cout << "Computing likelihoods" << endl;

//	ofstream plfile ("/home/chaitanya/Desktop/nancy/pedigreelikelihoods.txt", ios::out);
//	ofstream plsfile ("/home/chaitanya/Desktop/nancy/snpScores.txt", ios::out  );
//	if (func_call == 1)
//	{
//		ofstream plfile ("/home/chaitanya/Desktop/nancy/pedigreelikelihoods.txt", ios::out);
//		plfile << "SNPID" << "\t" << "FAMID" << "\t" << "LIKELIHOOD" << "\t" << "uninf-LIKELIHOOD" << "\t" << "RATIO" << "\t" << "CPU-time" << endl;
//		ofstream plsfile ("/home/chaitanya/Desktop/nancy/snpScores.txt", ios::out  );
//		plsfile << "SNPID" << "\t" << "SCORE" << endl;
//	}
//	else
//	{
////		plfile.close();
////		plsfile.close();
		ofstream plfile (filename+".pedigreelikelihoods", ios::out|ios::app  );
		ofstream plsfile (filename+".snpScores", ios::out|ios::app  );
//	}


	
	vector < vector <long double>> InfLikelihoods;
	auto vecp = UninfLikelihoods.begin();
//	auto fm=FMGLs.begin();
	auto fm=popGLs.begin();
	long double inflik_A;
	long double inflik_X;
	string plfamid;
	string plsnpid;
	//int snp_count = 0;
	//long int sec;
	//time_t seconds;
	//long double pedlrtA = 0;
	//long double pedlrtX = 0;
	//long double llr= 0;
	//long double pval = 0;
	
//	ofstream pedlrtfile (filename+"_pedigreeLRTstatistic.txt", ios::out|ios::app );


	for (auto snp = snps.begin() ; snp != snps.end(); snp++)
	{
		auto famchrompos = (*snp).begin();
		auto indchrompos = (*famchrompos).begin();
		int pedigree = 0;
		long double pscoreA = 0;
		long double pscoreX = 0;
		long double pedLAut = 0;
		long double pedLSex = 0;		
//		pedlrtfile << left << setw(10) << plsnpid;
		for (auto fam = (*snp).begin() ; fam != (*snp).end(); fam++)
		{
			vector<GLPROB> glpv;
			int famsize = 0;
			int hets = 0;
			int homs = 0;
			int unks = 0;
			//sec = clock();
			//seconds = time (NULL);
			for (auto ind = (*fam).begin(); ind!=(*fam).end(); ind++)
			{
				GLPROB glp;
				glp.setelem(*ind, fm->GLs, "INF", phredFLAG);
				glpv.emplace_back(glp);
				famsize++;
				plfamid = ind->familyid;
				plsnpid = ind->snpid;
				if (ind->sex == "1")
					homs++;
				else if (ind->sex == "2")
					hets++;
				else if (ind->sex == "0")
					unks++;				
			}
			string famkeyA = "A" + itoa(famsize) + "A";
			string famkeyX = "X" + itoa(famsize) + "X" + itoa(homs) + "X" + itoa(hets) + "X" + itoa(unks) + "X"	;		
			inflik_A = cart_product(glpv, Penetrance, "A");
			inflik_X = cart_product(glpv, Penetrance, "X");
//			plfile << left << setw(10) << plsnpid << setw(10) << plfamid << setw(15) << setprecision(6) << inflik_A << setw(15) << setprecision(6) <<  inflik_X << setw(15) << setprecision(6) << (*vecp)[famkeyA] << setw(15) << setprecision(6) << (*vecp)[famkeyX] << left << setw(15) << setprecision(6) << inflik_A/(*vecp)[famkeyA] << left << setw(15) << setprecision(6) << inflik_X/(*vecp)[famkeyX]<< right << setw(15) << setprecision(10) << ((long int) clock()-sec) << endl; //time(NULL) - seconds << endl;
			plfile <<  indchrompos->chromosome << "\t" <<  indchrompos->position << "\t" <<  plsnpid << "\t" <<  plfamid << "\t" <<  inflik_A << "\t" <<  inflik_X << "\t" <<  (*vecp)[famkeyA] << "\t" <<  (*vecp)[famkeyX] << "\t" << inflik_A/(*vecp)[famkeyA] << "\t" <<  inflik_X/(*vecp)[famkeyX]<< endl; //<<right << setw(15) << setprecision(10) <<  //((long int) clock()-sec) << "\n"; //time(NULL) - seconds << endl;
			//pedlrtA = log(inflik_A); ///(*vecp)[famkeyA]);
			//pedlrtX = log(inflik_X); ///(*vecp)[famkeyX]);
			pscoreA += log(inflik_A/(*vecp)[famkeyA]);
			pscoreX += log(inflik_X/(*vecp)[famkeyX]);
			pedLAut += log(inflik_A);
			pedLSex += log(inflik_X);			
			//llr = -2*(pedlrtA-pedlrtX);
//			pval = bchisqr(1, llr);
//			pedlrtfile << setw(15) << setprecision(6) << pedlrtA << setw(15) << setprecision(6) << pedlrtX << setw(15) << setprecision(6) << llr << setw(15) << setprecision(6) << pval;
//			pedlrtfile << left << setw(10) << plsnpid << "\t" << setw(10) <<  plfamid << "\t" << setw(15) << setprecision(6) << pedlrtA << "\t" << setw(15) << setprecision(6) << pedlrtX << "\t" << setw(15) << setprecision(6) << llr << "\t" << setw(15) << setprecision(6) << pval << "\n";
			pedigree++;
		}
		fm++;
		vecp++;
		//double LLRstatistic = -2*(pedLAut-pedLSex);
//		double alpha = 0.067;
		long double normL = (alpha * exp(pedLSex) + (1-alpha) * exp(pedLAut));
		long double posterior_Prob_Sex_linkage = (alpha * exp(pedLSex ))/normL;
//		long double denomratio = alpha/(alpha+(1-alpha)*exp(pedLAut - pedLSex));
//		plsfile << left  << plsnpid << right << setw(15) << setprecision(6) << pscoreA << right << setw(15) << setprecision(6) << pscoreX << right << setw(15) << setprecision(6) << LLRstatistic << right << setw(15) << setprecision(6) <<  pedigree << right << setw(15) << setprecision(6) <<  bchisqr(pedigree, LLRstatistic) << endl;
//		plsfile << plsnpid << "\t" << setw(15) << setprecision(6) << pscoreA << "\t" << setw(15) << setprecision(6) << pscoreX << "\t"  << setw(15) << setprecision(6) << pedLAut << "\t"  << setw(15) << setprecision(6) << pedLSex << "\t"<< setw(15) << setprecision(6) << LLRstatistic << "\t" << setw(15) << setprecision(6) <<  pedigree << "\t" << setw(15) << setprecision(6) <<  bchisqr(pedigree, LLRstatistic) << "\n";
		plsfile << indchrompos->chromosome << "\t" << indchrompos->position << "\t" <<  plsnpid << "\t" << pscoreA << "\t"  << pscoreX << "\t"  <<  exp(pedLAut) << "\t"  <<  exp(pedLSex) << "\t"<<  posterior_Prob_Sex_linkage   << endl;
//		cout << snp_count++ << endl;
//		pedlrtfile << "\n";
	}	
	plfile.close();
	plsfile.close();	
//	pedlrtfile.close();
//	cout << "Done!!!" << endl;

}

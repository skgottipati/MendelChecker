#include "fileread.h"
#include "computeLikelihood.h"


std::string getFileName(const std::string strPath)
{
	//size_t iLastSeparator = 0;
//	cout << strPath.substr(iLastSeparator = strPath.find_last_of("/")) != std::string::npos ? iLastSeparator + 1 : 0 << endl;
//	cout << strPath.size()  << "\t" << strPath.find_last_of(".") << endl;
//	return strPath.substr((iLastSeparator = strPath.find_last_of("/")) != std::string::npos ? iLastSeparator + 1 : 0, strPath.size() - strPath.find_last_of("."));
	return strPath.substr(0, strPath.find_last_of("."));
}


bool map_compare (std::map<std::string, std::string> const &lhs, std::map<std::string, std::string> const &rhs) {
    // No predicate needed because there is operator== for pairs already.
	return (lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin()));
}


//void fileread3(char* argv[], double alpha){
void fileread(string fname, int bufsize, string phredFLAG, double alpha, string unfFLAG, string outfilename){

	ifstream file (fname, ios::in|ios::ate);
	std::ifstream inFile(fname);
	//string filename = argv[1];
//	cout << filename << endl;
	string filename = fname;
//	cout << filename << "\t alpha: " << alpha << endl;
  	//long numlines =  std::count(std::istreambuf_iterator<char>(inFile), std::istreambuf_iterator<char>(), '\n');
	//int running_size = 0;
	//std::streamsize step_size = 300;
//	cout << "Number of lines : " << numlines << endl;
	//long size;

	//int bufsize = 0;
	int snp_count = 0;

	vector< vector <vector <LINE>>> snps;

	int pedlength;
	int snpsperloop=2;
	//cout << "numsnps:" << snpsperloop << endl;	
	//snps.reserve(snpsperloop);

	//auto it = snps.begin();	
//	cout <<  "snps capacity " << snps.capacity() << endl;
	//int i = -1;
	string chrom = "-1";
	int pos = -1;
	string mom = "-1";
	string dad = "-1";
	string famid = "-1";
	vector< vector< LINE>> founders;
//	vector<vector <vector <LINE>>> missing;
	vector< LINE> pl;
	pl.reserve(100);
	vector< vector <LINE>> pp;
	vector< LINE> q;
	//int r = 0;


//	defining genotypes

	std::vector<string> genotypes;
	genotypes.reserve(10);
	genotypes.insert(genotypes.end(),  {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"});



//	defining sex & autosome indicators

	std::vector<string> chroms;
	chroms.reserve(4);
	chroms.insert(chroms.end(),  {"A", "X1", "X2", "X0"});



//	generating transition probability map
	
//	cout << "Generating transition probabilities" << endl;

	std::unordered_map<std::string, double> Penetrance;

	for (auto chrom = chroms.begin(); chrom != chroms.end(); chrom++)
	{
		for (auto hom = genotypes.begin(); hom != genotypes.end(); hom++)
		{
			for (auto het = genotypes.begin(); het != genotypes.end(); het++)
			{
				for (auto off = genotypes.begin(); off != genotypes.end(); off++)
				{
					string hom_first = hom->substr(0,1);
					string hom_second = hom->substr(1,1);
					string het_first = het->substr(0,1);
					string het_second = het->substr(1,1);
					string off_first = off->substr(0,1);
					string off_second = off->substr(1,1);

					if (*chrom == "A")
					{	
						double countMatches = 0;
						if ((*off == hom_first+het_first) || (*off == het_first+hom_first))
							countMatches++;
						if ((*off == hom_first+het_second) || (*off == het_second+hom_first))
							countMatches++;
						if ((*off == hom_second+het_first) || (*off == het_first+hom_second))
							countMatches++;
						if ((*off == hom_second+het_second) || (*off == het_second+hom_second))
							countMatches++;	
						Penetrance.insert(std::make_pair(*hom+*het+*off+*chrom, countMatches/4 ));
					}
	//				cout << *mo + "," + *fa+ "," + *off << "\t" <<  countMatches/4	<< endl;
					if (*chrom == "X1")
					{	
						double countMatches = 0;
						if (het_first == het_second)
						{
							if ((*off == hom_first+het_first) || (*off == het_first+hom_first))
								countMatches++;
							if ((*off == hom_second+het_first) || (*off == het_first+hom_second))
								countMatches++;
						}
						Penetrance.insert(std::make_pair(*hom+*het+*off+*chrom, countMatches/2 ));

      //if ($GENOHOM{$k} eq $HetParentGametes[0].$HomParentGametes[0] || $GENOHOM{$k} eq $HomParentGametes[0].$HetParentGametes[0]){$countMatches++;}
      //if ($GENOHOM{$k} eq $HetParentGametes[0].$HomParentGametes[1] || $GENOHOM{$k} eq $HomParentGametes[1].$HetParentGametes[0]){$countMatches++;}

      //determine the probability of the offspring genotype conditional on the parental genotypes & mutation rate
      //if ($countMatches > 0){$transition=($countMatches/2-$mutation);}
					}
					if (*chrom == "X2")
					{	
						double countMatches = 0;
						if (het_first == het_second && off_first == off_second)
						{
							if (off_first == hom_first) 
								countMatches++;
							if (off_first == hom_second)
								countMatches++;
						}
						Penetrance.insert(std::make_pair(*hom+*het+*off+*chrom, countMatches/2 ));
      //if ($GENOHET{$k} eq $HomParentGametes[0]){$countMatches++;}
      //if ($GENOHET{$k} eq $HomParentGametes[1]){$countMatches++;}

      //determine the probability of the offspring genotype conditional on the parental genotypes & mutation rate
      //if ($countMatches > 0){$transition=($countMatches/2-$mutation);}
					}
					if (*chrom == "X0")
					{	
						double countMatches = (Penetrance[*hom+*het+*off+"X1"]+Penetrance[*hom+*het+*off+"X2"])/2;
						Penetrance.insert(std::make_pair(*hom+*het+*off+*chrom, countMatches ));
					}
				}
			}
		}
	}

	//ofstream trprfile ("transitionPr.txt", ios::out );
	//trprfile << "transition_key_(HomGT+HetGT+OffGT+A|X{0-2})" << "\t" << "Prob" << endl;
	//for (auto pos=Penetrance.begin(); pos!=Penetrance.end(); pos++)
	//	trprfile << pos->first << "\t" << pos->second << endl;
	//trprfile.close();


	//ofstream mfglfile (filename+"_meanFounderGL.txt", ios::out );
	//mfglfile << "SNPID" << "\t" << "FOUNDERSwData" << "\t" ;
	//for (auto gt = genotypes.begin(); gt != genotypes.end(); gt++)
	//{
	//	mfglfile << *gt << "\t";
	//}
	//mfglfile << "SUM_GL" << endl;
	//mfglfile.close();

	//ofstream popGTfile (filename+"_populationGL.txt", ios::out );
	//popGTfile << "SNPID" << "\t" << "FOUNDERSwData" << "\t";
	//for (auto gt = genotypes.begin(); gt != genotypes.end(); gt++)
	//{
	//	popGTfile << *gt << "\t";
	//}
	//popGTfile <<  "SUM_GL" << endl;
	//popGTfile.close();
	//
	//ofstream ulfile (filename+"_uninformativelikelihoods.txt", ios::out );
	//ulfile << "SNPID" << "\t" << "OFFSPRINGS" << "\t"  << "uninfL" << endl;
	//ulfile.close();

	ofstream plfile (outfilename+".pedigreelikelihoods", ios::out);
	plfile << "CHR" << "\t" << "POS" << "\t" << "SNP" << "\t" << "FAMID" << "\t" << "AutoL" << "\t" << "SexL" << "\t"<< "AutoUninfL" << "\t"<< "SexUninfL" << "\t" << "AutoRATIO" << "\t" << "SexRATIO" << endl;
	ofstream plsfile (outfilename+".snpScores", ios::out  );
//	plsfile << "SNPID" << "\t" << "AutoSCORE" <<  "\t" << "SexSCORE" << "\t" << "AutoPedL" << "\t" << "SexPedL" << "\t" <<"LRT" << "\t" << "dof" << "\t" << "Pvalue" << endl;
	plsfile << "CHR" << "\t" << "POS" << "\t" << "SNP" << "\t" << "AutoSCORE" <<  "\t" << "SexSCORE" << "\t" << "AutoPedL" << "\t" << "SexPedL" << "\t" <<"PP_sex" << endl;
	plfile.close();
	plsfile.close();
	//ofstream pedlrtfile (filename+"_pedigreeLRTstatistic.txt", ios::out);
	//pedlrtfile << "SNPID" << "\t" << "FAMID" << "\t" << "AutoLogL" << "\t" << "SexLogL" << "\t"<< "LRT" << "\t"<< "p-value" << endl;	
	//pedlrtfile.close();
//
//
//
//
	//cout << "here" << endl;
	std::map<std::string, std::string> first_pedigree;
	std::map<std::string, std::string> iter_pedigree;
	int firstpedflag = 0;
	if (file.is_open())
	{
		file.seekg (0, ios::beg);
//		cout << "file is open" << endl;
		while(1)
		{
			string line;
			if (!getline(file, line))
			{
				if (firstpedflag != 0)
				{
					try
					{
						if (!map_compare(first_pedigree, iter_pedigree))
						{
							//for (auto ped = first_pedigree.begin(); ped!= first_pedigree.end(); ped++)
							//	cout << ped->first << " : " << ped->second << endl;
							//cout << "\n\n" << endl;
							//for (auto ped = iter_pedigree.begin(); ped!= iter_pedigree.end(); ped++)
							//	cout << ped->first << " : " << ped->second << endl;								
							//cout << "first" << endl;
							throw "Pedigrees across SNP's do not match.";
						}
					}
					catch(const char* Message)
					{
						cerr << "Error: " << Message << endl;
						exit (EXIT_FAILURE);
					}							
				}
				else
				{					
					verify_pedigrees(first_pedigree);
					pedlength = first_pedigree.size();
					snpsperloop = (1024/(double) pedlength)*1024*1024*((double) bufsize/1000);
					snps.reserve(snpsperloop);
				}
				firstpedflag +=1;
				iter_pedigree.clear();					
				//cout << line << "break" << endl;
				pp.emplace_back(pl);
				snps.emplace_back(pp);
				founders.emplace_back(q);
				if (snps.empty())
				{
					cerr << "No SNP's found" << endl;
					exit (EXIT_FAILURE);
				}				
				new_compute_likelihood(snps, founders, Penetrance, outfilename, alpha, unfFLAG, phredFLAG);
				break;
			}
			else if (line.empty())
			{
				continue;
			}
			else
			{
//				getline(file,line);
//				cout << line << endl;
				LINE l;
				l.setelem(line);
				if (line.substr(0,1) != "#")
				{
					vector<string> fields = split(line, '\t');
					stringstream ss1,ss2;
					ss1 << fields[3] << ":" << fields[4];
					ss2 << fields[3] << "\t" << fields[4] << "\t" << fields[5] << "\t" << fields[6] << "\t" << fields[7];
					
					//cout << firstpedflag << endl;
					//cout << first_pedigree.size() << ":" << iter_pedigree.size() << endl;
					if (l.chromosome == chrom && l.position == pos)
					{
						if (firstpedflag != 0)
						{
							//cout << "first" << endl;
							iter_pedigree.insert(std::make_pair(ss1.str(), ss2.str()));
						}
						else
						{
							first_pedigree.insert(std::make_pair(ss1.str(), ss2.str()));
						}
						
	//					cout << "snps capacity "<< snps.size() << " " << line << endl;
						if (l.familyid == famid)
						{					
							pl.emplace_back(l);
						}
						else
						{
							pp.emplace_back(pl);
							pl.clear();
							pl.emplace_back(l);
							famid = l.familyid;
						}
						if (l.mom == "0" && l.dad == "0")
						{
							q.emplace_back(l);
						}
						bufsize += line.size();
						continue;
					}
					else
					{
						if (pos != -1)
						{
							pp.emplace_back(pl);
							snps.emplace_back(pp);
							snp_count++;
							pl.clear();
							pp.clear();
							founders.emplace_back(q);
							q.clear();
							
							if (firstpedflag != 0)
							{
								//cout << "hi" << endl;
								try
								{
									if (!map_compare(first_pedigree, iter_pedigree))
									{
										//cout << "second" << endl;
										throw "Pedigrees across SNP's do not match.";
									}
								}
								catch(const char* Message)
								{
									cerr << "Error: " << Message << endl;
									exit (EXIT_FAILURE);
								}							
							}
							else
							{
								//for (auto ped=first_pedigree.begin(); ped!=first_pedigree.end(); ped++)
								//	cout << ped->second <<endl;
								//cout << endl;
								verify_pedigrees(first_pedigree);
								pedlength = first_pedigree.size();
								snpsperloop = (1024/(double) pedlength)*1024*1024*((double) bufsize/1000);
								snps.reserve(snpsperloop);								
							}
							firstpedflag +=1;
							iter_pedigree.clear();					

						}
//						cout << "snps capacity "<< snps.size() << " " << line << endl;
						if (l.mom == "0" && l.dad == "0")
						{
							q.emplace_back(l);

						}

//							cout << l.chromosome << " " << l.position << " " << l.snpid << endl;
						pl.emplace_back(l);
						famid = l.familyid;


						chrom = l.chromosome;
						pos = l.position;

						if (firstpedflag != 0)
						{
							//cout << "first" << endl;
							iter_pedigree.insert(std::make_pair(ss1.str(), ss2.str()));
						}
						else
						{
							first_pedigree.insert(std::make_pair(ss1.str(), ss2.str()));
						}
						
						bufsize = 0;
						bufsize += line.size();
					}
					if (snps.size() == (unsigned int) snpsperloop)
					{
//						cout << "snp count " << snpsperloop << endl;
//						pp.emplace_back(pl);
//						snps.emplace_back(pp);
//						founders.emplace_back(q);


//						for (vector< vector< vector<int> >>::size_type i = 0; i < snps.size(); i++){
//							cout << i <<  "\t" << snps[i].size() << "\t";		
//							for (vector< vector<int> >::size_type j = 0; j < snps[i].size(); j++){
//								cout << snps[i][j].size() << "\t";
//								for (vector<int>::size_type k=0; k < snps[i][j].size(); k++){
//									cout << snps[i][j][k].individualid << ":" << snps[i][j][k].snpid<< "\t";
//								}
//							}
//							cout <<  endl;
//						}
						snp_count = 0;
						new_compute_likelihood(snps, founders, Penetrance, outfilename, alpha, unfFLAG, phredFLAG);
						pp.clear();
						snps.clear();
						founders.clear();
		//				pl.clear();
		//				q.clear();
					}
				}
				else
					bufsize += line.size();
			}
		}
		file.close();
	}
	else
	{
		cout << "Warning: cannot open the file!" << endl;
	}
//
////						for (vector< vector< vector<int> >>::size_type i = 0; i < snps.size(); i++){
////							cout << i <<  "\t" << snps[i].size() << "\t";		
////							for (vector< vector<int> >::size_type j = 0; j < snps[i].size(); j++){
////								cout << snps[i][j].size() << "\t";
////								for (vector<int>::size_type k=0; k < snps[i][j].size(); k++){
////									cout << snps[i][j][k].individualid << ":" << snps[i][j][k].snpid << "\t";
////								}
////							}
////							cout <<  endl;
////						}
//
//
//
//
//

}

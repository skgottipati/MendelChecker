#include "genoped.h"


vector<string> splitdelim(string line){
    
    vector<string > fields;
    std::istringstream linestream (line, istringstream::in);
    string w;
    while(1){
        getline(linestream, w, '\t');
        fields.emplace_back(w);
    }
    return fields;
}


std::map<std::string, std::string> pedigree_reader(string pedfilename){
    
    //std::ifstream inFile(pedfilename);
    //long numlines =  std::count(std::istreambuf_iterator<char>(inFile), std::istreambuf_iterator<char>(), '\n');  
    ifstream file (pedfilename, ios::in|ios::ate);
    std::map<std::string, std::string> pedigree_str;
    if (file.is_open()){
    	file.seekg (0, ios::beg);
        cout << "Reading pedigree file." << endl;
        string line;
        while(getline(file, line))
        {
                if (line.substr(0,1) != "#"){
                    vector<string> fields = split(line, '\t');
                    std::stringstream ss2;
                    ss2 << fields[3] << "\t" << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[4];
                    std::string s2 = ss2.str();
                    std::stringstream ss1;
                    ss1 << fields[3] << ":" << fields[0];
                    std::string s1 = ss1.str();
                    pedigree_str.insert(std::make_pair(s1, s2));
                    //cout << s1 << "\t" << s2 << endl;
                }
        }
        file.close();
    }
    else
    {
        cout << "Warning: cannot open pedigree file!" << endl;
        file.close();
        exit (EXIT_FAILURE);
    }
    return pedigree_str;
}


vector<string> get_individuals_from_pedigreemap(const std::map<std::string, std::string>* pedigree){
    vector<string> fams;
    vector<string> samples;
    vector<string> fields;
    for (auto ped=pedigree->begin(); ped!=pedigree->end(); ped++){
        fields = split(ped->first, ":");
        fams.emplace_back(fields.at(0));
        samples.emplace_back(fields.at(1));
    }
    return samples;
}



// F(j/k) = (k*(k+1)/2)+j
std::unordered_map<std::string, int> generateGTmap(vector<string> alleles){
    int j = 0;
    int k = 0;
    std::unordered_map<std::string, int> GTmap;
    for (auto ind1=alleles.begin(); ind1!=alleles.end(); ind1++){
        k = 0;
        for (auto ind2=alleles.begin(); ind2!=alleles.end(); ind2++){
            std::stringstream ss;
            if (*ind1 > *ind2)
                ss << *ind2 << *ind1;
            else
                ss << *ind1 << *ind2;
            std::string s = ss.str();
            GTmap.insert(std::make_pair(s, (k*(k+1)/2)+j));
            k = k + 1;
        }
        j = j + 1;
    }
    return GTmap;
}




std::unordered_map<std::string, double> computePenetrance(){
    
//	defining genotypes

	std::vector<string> genotypes;
	genotypes.reserve(10);
	genotypes.insert(genotypes.end(),  {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"});


//	defining sex & autosome indicators

	std::vector<string> chroms;
	chroms.reserve(4);
	chroms.insert(chroms.end(),  {"A", "X1", "X2", "X0"});    
    
//	generating transition probability map
    
    cout << "Generating transition probabilities." << endl;

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
    return Penetrance;
}


std::string get_FileName(const std::string strPath)
{
	//size_t iLastSeparator = 0;
//	cout << strPath.substr(iLastSeparator = strPath.find_last_of("/")) != std::string::npos ? iLastSeparator + 1 : 0 << endl;
	cout << strPath.size()  << "\t" << strPath.find_last_of(".") << endl;
//	return strPath.substr((iLastSeparator = strPath.find_last_of("/")) != std::string::npos ? iLastSeparator + 1 : 0, strPath.size() - strPath.find_last_of("."));
	return strPath.substr(0, strPath.find_last_of("."));
}




void read_geno(string genofield, int numsnps, string vcfname, const std::map<std::string, std::string>* pedigree, std::unordered_map<std::string, double> Penetrance, double alpha, string unfFLAG){
    ifstream file (vcfname, ios::in|ios::ate);
    std::unordered_map<std::string, std::string> pedigree_str;
    std::vector<string> genotypes;
    genotypes.reserve(10);
    genotypes.insert(genotypes.end(),  {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"});    
    int snpcount = 0;
    vector<string> samples;
    vector< vector <vector <LINE>>> snps;
    snps.reserve(numsnps);
    vector< vector< LINE>> founders;
    vector< LINE> pl;
    vector< vector <LINE>> pp;
    vector< LINE> q;
    string famid;
    
    string filename = vcfname;
    //ofstream mfglfile (filename+"_meanFounderGL.txt", ios::out );
    //mfglfile << "SNPID" << "\t" << "FOUNDERSwData" << "\t" ;
    //for (auto gt = genotypes.begin(); gt != genotypes.end(); gt++)
    //{
    //        mfglfile << *gt << "\t";
    //}
    //mfglfile << "SUM_GL" << endl;
    //mfglfile.close();

    //ofstream popGTfile (filename+"_populationGL.txt", ios::out );
    //popGTfile << "SNPID" << "\t" << "FOUNDERSwData" << "\t";
    //for (auto gt = genotypes.begin(); gt != genotypes.end(); gt++)
    //{
    //        popGTfile << *gt << "\t";
    //}
    //popGTfile <<  "SUM_GL" << endl;
    //popGTfile.close();
    //
    //ofstream ulfile (filename+"_uninformativelikelihoods.txt", ios::out );
    //ulfile << "SNPID" << "\t" << "OFFSPRINGS" << "\t"  << "uninfL" << endl;
    //ulfile.close();

    ofstream plfile (filename+"_pedigreelikelihoods.txt", ios::out);
    plfile << "SNPID" << "\t" << "FAMID" << "\t" << "AutoL" << "\t" << "SexL" << "\t"<< "AutoUninfL" << "\t"<< "SexUninfL" << "\t" << "AutoRATIO" << "\t" << "SexRATIO"<< "\t" << "CPUtime" << endl;
    ofstream plsfile (filename+"_snpScores.txt", ios::out  );
//	plsfile << "SNPID" << "\t" << "AutoSCORE" <<  "\t" << "SexSCORE" << "\t" << "AutoPedL" << "\t" << "SexPedL" << "\t" <<"LRT" << "\t" << "dof" << "\t" << "Pvalue" << endl;
    plsfile << "SNP" << "\t" << "AutoSCORE" <<  "\t" << "SexSCORE" << "\t" << "AutoPedL" << "\t" << "SexPedL" << "\t" <<"PP_sex" << endl;
    plfile.close();
    plsfile.close();
        
        
    if (file.is_open()){
    	file.seekg (0, ios::beg);
        cout << "Reading VCF file." << endl;
        string line;
        while(1){
            if(!getline(file, line)){
//                cout << line << "break" << endl; 						
                //pp.emplace_back(pl);
                //snps.emplace_back(pp);
                //founders.emplace_back(q);
                new_compute_likelihood(snps, founders, Penetrance, vcfname, alpha, unfFLAG);
                break;                
            }
            else {
                if (line.substr(0,2) == "##"){
                    line.clear();
                    continue;
                }
                else if (line.substr(0,1) == "#"){
                    vector<string> fields = split(line, '\t');
                    samples.reserve(fields.size()-9);
                    samples.assign(fields.begin()+9, fields.end());
                    cout << samples.size() << " samples were found in the VCF file." << endl;
                    for (auto sam = samples.begin(); sam!= samples.end(); sam++ )
                        cout << *sam << "\t";
                    cout << endl;
                    line.clear();
                    fields.clear();
                }
                else
                {
                    //cout << line << endl;
                    vector<string> fields = split(line, '\t');
                    string chrom = fields.at(0);
                    string pos = fields.at(1);
                    string snpid = fields.at(2);
                    string ref = fields.at(3);
                    vector<string> alt = split(fields.at(4), ",");
                    vector <string> alleles;
                    alleles.emplace_back(ref);
                    alleles.resize(alt.size()+1);
                    std::copy(alt.begin(), alt.end(), alleles.begin()+1);
                    vector<string> format = split(fields.at(8), ':');
                    int GTindex = -1;
                    int PLindex = -1;
                    int runind = 0;
                    for (auto ind=format.begin(); ind != format.end(); ind++)
                    {
                        if (*ind == "GT")
                            GTindex = runind;
                        if (*ind == genofield)
                            PLindex = runind;
                        runind = runind + 1;
                    }
                    if (GTindex ==-1 || PLindex == -1)
                    {
                        cout << "Genotypes or genotype likelihood do not exist in the format" << endl;
                        continue;
                    }
                    std::unordered_map<std::string, int> GTs = generateGTmap(alleles);
                    std::unordered_map<std::string, std::string> sample_GLs;
                    vector<string> sample_fields;
                    vector<string> GLvec;
                    std::map<string,string> GL;
                    int samind = 0;
                    //cout << "fields_size:" << fields.size() << endl;
                    for (auto ind=fields.begin()+9; ind!=fields.end(); ind++)
                    {
                        sample_fields = split(*ind, ':');
                        stringstream ss;
                        if (sample_fields.at(0) == "./." ){
                            //cout << " was here" << endl;
                            for (int i=0; i<10; i++)
                                GL.insert(std::make_pair(genotypes[i], "-1" ));
                        }
                        else {
                            GLvec = split(sample_fields.at(PLindex), ",");
                            for (int i=0; i<10; i++)
                                GL.insert(std::make_pair(genotypes[i], "" ));
                            for (auto it=GTs.begin(); it!=GTs.end(); it++)
                            {
                                GL.at(it->first) = GLvec.at(it->second);
//                                cout << it->first << ":" << it->second << ":" << GLvec.at(it->second) << endl;
                            }
                            GLvec.clear();
                        }
                        for (auto it=GL.begin(); it!=GL.end(); it++)
                            ss << it->second << "\t";
                        string s = ss.str();
                        //cout << s << endl;
                        sample_GLs.insert(std::make_pair(samples.at(samind), s));
                        samind = samind + 1;
                        GL.clear();
                        sample_fields.clear();
                    }
                    fields.clear();
                    snpcount = snpcount + 1;
                    //for (auto it=sample_GLs.begin(); it!=sample_GLs.end(); it++)
                    //    cout << it->first << "\t" << it->second << "\t" << ": " <<sample_GLs.size() << endl;
                    int pedcount = 0;
                    for (auto ped=pedigree->begin(); ped!=pedigree->end(); ped++){
//                        cout << ped->first << endl;
                        vector<string> pedkey = split(ped->first, ":");
                        std::stringstream sspedsnp;
                        sspedsnp << chrom << "\t" << pos << "\t" << snpid << "\t" << ped->second << "\t" << sample_GLs.at(pedkey.at(1));
//                        cout << sspedsnp.str() << endl;
                        string line = sspedsnp.str();
                        LINE l;
                        l.setelem(line);
//                        cout << l.familyid << "\t" << l.chromosome << "\t" << l.position << endl;
                        if (pedcount == 0)
                            famid = l.familyid;
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
                        pedcount = pedcount + 1;
                        //bufsize += line.size();
                        //continue;                    
                    }
                    pp.emplace_back(pl);
                    snps.emplace_back(pp);
                    founders.emplace_back(q);
                    pp.clear();
                    pl.clear();
                    q.clear();
                    //cout << snpcount << "\t" << numsnps << "\t" << PLindex << endl;
                    if (snpcount == numsnps){
                        //cout << "was here" << endl;
                        new_compute_likelihood(snps, founders, Penetrance, vcfname, alpha, unfFLAG);
                        pp.clear();
                        snps.clear();
                        founders.clear();
                        line.clear();
                        snpcount = 0;
                    }
                    line.clear();
                }
            }
        }
        //if (file.eof())
        //    return -1;
        //else
        //    return file.tellg();
    }
    else
    {
        cout << "Warning: cannot open the file!" << endl;
        exit (EXIT_FAILURE);
        //return 0;
    }

    //for (vector< vector< vector<int> >>::size_type i = 0; i < snps.size(); i++){
    //    cout << i <<  "\t" << snps[i].size() << "\t";		
    //    for (vector< vector<int> >::size_type j = 0; j < snps[i].size(); j++){
    //        cout << snps[i][j].size() << "\t";
    //        for (vector<int>::size_type k=0; k < snps[i][j].size(); k++){
    //            cout << snps[i][j][k].individualid << ":" << snps[i][j][k].snpid << "\t";
    //        }
    //    }
    //    cout <<  endl;
    //}

    snps.clear();
    founders.clear();
    
    file.close();
}


//    
//    vector< vector <vector <LINE>>> snps;
//
//    int snpsperloop = 100;
//    snps.reserve(snpsperloop);
//
//    auto it = snps.begin();	
//    cout <<  "snps capacity " << snps.capacity() << endl;
//    int i = -1;
//    string chrom = "-1";
//    int pos = -1;
//    string mom = "-1";
//    string dad = "-1";
//    string famid = "-1";
//    vector< vector< LINE>> founders;
////	vector<vector <vector <LINE>>> missing;
//    vector< LINE> pl;
//    pl.reserve(100);
//    vector< vector <LINE>> pp;
//    vector< LINE> q;
////  Chr     Pos     SNPId   FamId   IndvId   DadId     MomId   Sex    
//    for (auto ped=pedigree.begin(); ped!=pedigree.end(); ped++){
//        vector<string> pedkey = split(ped->first, ":");        
//        std::stringstream ss;
//        ss << chrom << "\t" << pos << "\t" << snpid << "\t" << ped->second << sample_GLs.at(pedkey.at(1));        
//        string line = ss.str();
//        LINE l;
//	l.setelem(line);
//        if (l.familyid == famid)
//        {                                   
//                pl.emplace_back(l);
//        }
//        else
//        {
//                pp.emplace_back(pl);
//                pl.clear();
//                pl.emplace_back(l);
//                famid = l.familyid;
//        }
//        if (l.mom == "0" && l.dad == "0")
//        {
//                q.emplace_back(l);
//        }
//        bufsize += line.size();
//        continue;       
//    }





//
//
//    			string line;
//			if (!getline(file, line))
//			{
//				cout << line << "break" << endl; 						
//				pp.emplace_back(pl);
//				snps.emplace_back(pp);
//				founders.emplace_back(q);
//				new_compute_likelihood(snps, founders, Penetrance, filename, alpha, unfFLAG);
//				break;
//			}
//			else
//			{
////				getline(file,line);
////				cout << line << endl;
//				LINE l;
//				l.setelem(line);
//				if (line.substr(0,1) != "#")
//				{					
//					if (l.chromosome == chrom && l.position == pos)
//					{
//	//					cout << "snps capacity "<< snps.size() << " " << line << endl;
//
//						if (l.familyid == famid)
//						{					
//							pl.emplace_back(l);
//						}
//						else
//						{
//							pp.emplace_back(pl);
//							pl.clear();
//							pl.emplace_back(l);
//							famid = l.familyid;
//						}
//						if (l.mom == "0" && l.dad == "0")
//						{
//							q.emplace_back(l);
//						}
//						bufsize += line.size();
//						continue;
//					}
//					else
//					{
//						if (pos != -1)
//						{
//							pp.emplace_back(pl);
//							snps.emplace_back(pp);
//							snp_count++;
//							pl.clear();
//							pp.clear();
//							founders.emplace_back(q);
//							q.clear();
//
//						}
////						cout << "snps capacity "<< snps.size() << " " << line << endl;
//						if (l.mom == "0" && l.dad == "0")
//						{
//							q.emplace_back(l);
//
//						}
//
////							cout << l.chromosome << " " << l.position << " " << l.snpid << endl;
//						pl.emplace_back(l);
//						famid = l.familyid;
//
//
//						chrom = l.chromosome;
//						pos = l.position;
//
//
//						bufsize = 0;
//						bufsize += line.size();
//					}
//					if (snps.size() == snpsperloop)
//					{
//						cout << "snp count " << snpsperloop << endl;
////						pp.emplace_back(pl);
////						snps.emplace_back(pp);
////						founders.emplace_back(q);
//
//
////						for (vector< vector< vector<int> >>::size_type i = 0; i < snps.size(); i++){
////							cout << i <<  "\t" << snps[i].size() << "\t";		
////							for (vector< vector<int> >::size_type j = 0; j < snps[i].size(); j++){
////								cout << snps[i][j].size() << "\t";
////								for (vector<int>::size_type k=0; k < snps[i][j].size(); k++){
////									cout << snps[i][j][k].individualid << ":" << snps[i][j][k].snpid<< "\t";
////								}
////							}
////							cout <<  endl;
////						}
//
//
//						snp_count = 0;
//						new_compute_likelihood(snps, founders, Penetrance, filename, alpha, unfFLAG);
//						pp.clear();
//						snps.clear();
//						founders.clear();
//		//				pl.clear();
//		//				q.clear();
//					}
//				}
//				else
//					bufsize += line.size();					
//			}
//
//		}
//}



//
//void parse(string& line, bool parseSamples) {
//
//    // clean up potentially variable data structures
//    info.clear();
//    infoFlags.clear();
//    format.clear();
//    alt.clear();
//    alleles.clear();
//
//    // #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE1 .. SAMPLEN]
//    vector<string> fields = split(line, '\t');
//    if (fields.size() < 7) {
//        cerr << "broken VCF record (less than 7 fields)" << endl
//             << line << endl;
//        exit(1);
//    }
//
//    sequenceName = fields.at(0);
//    char* end; // dummy variable for strtoll
//    position = strtoll(fields.at(1).c_str(), &end, 10);
//    id = fields.at(2);
//    ref = fields.at(3);
//    alt = split(fields.at(4), ","); // a comma-separated list of alternate alleles
//
//    // make a list of all (ref + alts) alleles, allele[0] = ref, alleles[1:] = alts
//    // add the ref allele ([0]), resize for the alt alleles, and then add the alt alleles
//    alleles.push_back(ref);
//    alleles.resize(alt.size()+1);
//    std::copy(alt.begin(), alt.end(), alleles.begin()+1);
//
//    // set up reverse lookup of allele index
//    altAlleleIndexes.clear();
//    int n = 0;
//    for (vector<string>::iterator a = alt.begin();
//            a != alt.end(); ++a, ++n) {
//        altAlleleIndexes[*a] = n;
//    }
//
//    convert(fields.at(5), quality);
//    filter = fields.at(6);
//    if (fields.size() > 7) {
//        vector<string> infofields = split(fields.at(7), ';');
//        for (vector<string>::iterator f = infofields.begin(); f != infofields.end(); ++f) {
//            if (*f == ".") {
//                continue;
//            }
//            vector<string> kv = split(*f, '=');
//            if (kv.size() == 2) {
//                split(kv.at(1), ',', info[kv.at(0)]);
//            } else if (kv.size() == 1) {
//                infoFlags[kv.at(0)] = true;
//            }
//        }
//    }
//    // check if we have samples specified
//    // and that we are supposed to parse them
//    if (parseSamples && fields.size() > 8) {
//        format = split(fields.at(8), ':');
//        // if the format changed, we have to rebuild the samples
//        if (fields.at(8) != lastFormat) {
//            samples.clear();
//            lastFormat = fields.at(8);
//        }
//        vector<string>::iterator sampleName = sampleNames.begin();
//        vector<string>::iterator sample = fields.begin() + 9;
//        for (; sample != fields.end() && sampleName != sampleNames.end(); ++sample, ++sampleName) {
//            string& name = *sampleName;
//            if (*sample == "." || *sample == "./.") {
//                samples.erase(name);
//                continue;
//            }
//            vector<string> samplefields = split(*sample, ':');
//            vector<string>::iterator i = samplefields.begin();
//            if (samplefields.size() != format.size()) {
//                // ignore this case... malformed (or 'null') sample specs are caught above
//                // /*
//                // cerr << "inconsistent number of fields for sample " << name << endl
//                // << "format is " << join(format, ":") << endl
//                // << "sample is " << *sample << endl;
//                // exit(1);
//                // *
//            }
//            else {
//                for (vector<string>::iterator f = format.begin(); f != format.end(); ++f) {
//                    samples[name][*f] = split(*i, ','); ++i;
//                }
//            }
//        }
//        if (sampleName != sampleNames.end()) {
//            cerr << "error: more sample names in header than sample fields" << endl;
//            cerr << "samples: " << join(sampleNames, " ") << endl;
//            cerr << "line: " << line << endl;
//            exit(1);
//        }
//        if (sample != fields.end()) {
//            cerr << "error: more sample fields than samples listed in header" << endl;
//            cerr << "samples: " << join(sampleNames, " ") << endl;
//            cerr << "line: " << line << endl;
//            exit(1);
//        }
//    } else if (!parseSamples) {
//        originalLine = line;
//    }
//
//    //return true; // we should be catching exceptions...
//}


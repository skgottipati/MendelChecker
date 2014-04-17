#include "verify_snps.h"


//class LINE {
//
//	public:
//	std::map<std::string,int> GLs;	
//	void setelem(string);	
//	string chromosome;
//	int position;
//	string snpid;
//	string familyid;
//	string individualid;
//	string mom;
//	string dad;
//	string sex;	
//};

void verify_nuclear_family(std::vector<std::string> fam){
    std::map<std::string, std::string> parents;
    std::map<std::string, std::map<std::string, std::string>> offsprings;
    vector<string> fields;
    std::map<std::string, std::string> offspring_fields;
    string mom = "mom";
    string dad = "dad";
    string sex = "sex";
    string parsex = "";
    string parind = "";
    //cout << fam.size() << endl;
    for (auto it = fam.begin(); it!= fam.end(); it++){
        fields = split(*it, "\t");
        //cout << fields.at(1) << endl;
        if (fields.at(0) =="")
            throw "All individuals must belong to a family";
        if (fields.at(2) == "0" && fields.at(3) == "0")
        {    
            parents.insert(std::make_pair(fields.at(1), fields.at(4)));
            if (parsex != "" && parsex == fields.at(4))
                throw "Both parents are of the same sex";
            if (parind != "" && parind == fields.at(1))
                throw "Both mother and father have the same id";
            if (fields.at(4) != "1" && fields.at(4) != "2")
            {
                //cout << fields.at(4) << endl;
                throw "All parents must have sex";
            }
            if (parsex == "")
                parsex = fields.at(4);
            if (parind == "")
                parind = fields.at(1);                
        }
        else
        {
            offspring_fields.insert(make_pair(mom, fields.at(3)));
            offspring_fields.insert(make_pair(dad, fields.at(2)));
            offspring_fields.insert(make_pair(sex, fields.at(4)));
            offsprings.insert(make_pair(fields.at(1), offspring_fields));
        }
    }
    //cout << "parents-size: " << parents.size() << "\t offsprings-size: " << offsprings.size() << endl;
    if (parents.size()<2)
        throw "One or both parents are missing";
    else if (parents.size()>2)
        throw "More than 2 parents exist in this family";
    
    parents.clear();
    offsprings.clear();
}




void verify_pedigrees(std::map<std::string, std::string> pedigree){
    string famid = "";
    vector<std::string> fam;
    vector<string> famidped;
    unsigned int pedind = 1;
    //cout << pedigree.size() << endl;
    cout << "Running verification of pedigree file." << endl;
    
    for (auto it=pedigree.begin(); it!=pedigree.end(); it++){
        famidped = split(it->first, ":");
        if (famid == famidped.at(0))
        {
            fam.emplace_back(it->second);
            //cout << "was heres" << endl;
            if (pedind == pedigree.size())
            {
                //cout << "was herer" << endl;
                try
                {
                    verify_nuclear_family(fam);
                }
                catch(const char* Message)
                {
                    cerr << "Error: " << Message << endl;
                    exit(EXIT_FAILURE);
                }
            }   
        }
        else
        {
            if (famid != "")
            {
                cout << "was herer" << endl;
                try
                {
                    //cout << "Running verification of pedigree file." << endl;
                    verify_nuclear_family(fam);
                    //cout << "Pedigree file passed all filters, however there potentially could be other problems with input files." << endl;
                }
                catch(const char* Message)
                {
                    cout << "Error: " << Message << endl;
                    exit(EXIT_FAILURE);
                }
                famid = famidped.at(0);
                fam.clear();
                fam.emplace_back(it->second);
            }
            else if (famid == "")
            {
                fam.emplace_back(it->second);
                famid = famidped.at(0);
            }
        }
        famidped.clear();
        pedind = pedind + 1;
    }
    cout << "Pedigree file passed all filters, however there potentially could be other format problems with input files." << endl;
}





//void verify_snps(vector< vector <vector <LINE>>> snps, vector< vector< LINE>> founders)
//{
//    for (auto snp = snps.begin() ; snp != snps.end(); snp++)
//    {
//        for (auto fam = (*snp).begin() ; fam != (*snp).end(); fam++)
//        {
//            vector<string > dads;
//            vector<string > moms;
//            vector<string > famids;
//            for (auto ind = (*fam).begin(); ind!=(*fam).end(); ind++)
//            {
//                
//                #####General notes
//                #pedigree fields are FAM ID FA MO SEX
//
//                
//                #all parents must have a sex
//                if (ind->mom == "0" && ind->dad =="0")
//                    if ind->sex == "2"
//                        throw "All parents must have a sex";
//                
//                #all offspring must have exactly two parents
//                if (ind->mom != "0" && ind->dad != "0")
//                {
//                    if (ind->familyid +ind->mom == ind->familyid +ind->dad)
//                        throw "Mom and dad cannot be the same individual";
//                    moms.emplace_back(ind->mom);
//                    dads.emplace_back(ind->dad);
//                }
//                #all individuals must have a family
//                if (ind->familyid == "0")
//                    throw "All individuals must have a family";
//
//                #duplicate IDs are allowed, so the unique identifier in a pedigree is FAM and ID
//                #FAM and ID accomodates splitting multiple generations into nuclear families
//                ##### pseudocode for problem with pedigree
//
//                # some definitions
//                # FAMID = concatenate text in FAM and ID fields (required for checks to accomodate multigenerational families that have been split up 'by hand')
//                # FAMFA = concatenate FAM and FA
//                # FAMMO = concatenate FAM and MO
//                # FAMID(FAMMO) SEX is the value of the variable SEX listed on the line where the variable ID matches the value MO
//
//                # set problem flag, if this doesn't change from false, the pedigree is legal
//                # problem = false
//    
//                #check for duplicate individuals
//                #foreach FAMID, if FAMID exists in the list of all FAMIDs more than once THEN problem = true & output line
//                famids.emplace_back(ind->familyid + ind->individualid);
//                
//                
//                
//                
//                
//                #check for parent violations
//                #foreach FAMID, if FAMFA is defined AND FAMMO is defined {
//                #individual is an offspring
//                    #if FAMFA == FAMMO THEN problem = true & output line # there are identical parents
//                    #if FAMID != FAMID(MO) FAM or FAM != FAMID(FAMFA) FAM THEN problem = true & output line
//                        #the FAM of offspring must match the FAM of parent
//                #checks for FAMFA
//                    #if FAMFA does not exist as an FAMID THEN problem = true & output line #parent must exist elsewhere in pedigree
//                    #for FAMID(FAMFA) #the line corresponding to FAMFA, i.e. lookup the line for which the value of FAMFA (collected from the line) is the FAMID (from the list of all individuals) 
//                        #if FAMID(FAMFA) SEX is not male THEN problem = true & output line #sex of parent must match sex of FAMID
//                        #if FAMID(FAMFA) is defined OR FAMID(FAMMO) is defined THEN problem = true & output line #the parents of any offspring cannot be offspring themselves, i.e. multigenerational families
//                #checks for MO
//                    #if MO does not exist as an FAMID THEN problem = true & output line #parent must exist elsewhere in pedigree
//                    #for FAMID(FAMMO) #the line corresponding to FAMMO, i.e. lookup the line for which the value of FAMMO (collected from the line) is the FAMID (from the list of all individuals) 
//                        #if FAMID(FAMMO) SEX is not female THEN problem = true & output line #sex of parent must match sex of FAMID
//                        #if FAMID(FAMFA) is defined OR FAMID(FAMMO) is defined THEN problem = true & output line #the parents of any offspring cannot be offspring themselves, i.e. multigenerational families
//                #} ELSE {
//                    #if FAMFA is undefined OR FAMMO is undefined THEN problem = true & output line #only one parent defined, the individual cannot be a legal parent or offspring
//                #}
//#if every FAMID passes this set of checks, it is a legal pedigree
//
//            }
//            if (moms.size() != 1 || dads.size() != 1)
//                throw "All offsprings do not belong to the same family.";
//                
//        }
//            
//    }
//
//
//
//
//}

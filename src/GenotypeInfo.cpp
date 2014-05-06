#include "GenotypeInfo.h"


void LINE::setelem(string line){
	std::istringstream linestream (line, istringstream::in);
	string w;
	getline (linestream, w, '\t');
	chromosome = w;
	getline (linestream, w, '\t');
	position = atoi(w.c_str());
	getline (linestream, w, '\t');
	snpid = w;
	getline (linestream, w, '\t');
	familyid = w;
	getline (linestream, w, '\t');
	individualid = w;
	getline (linestream, w, '\t');
	dad = w;
	getline (linestream, w, '\t');
	mom = w;
	getline (linestream, w, '\t');
	sex = w;
	int i;
	std::vector<string> genotypes;
	genotypes.reserve(10);
	genotypes.insert(genotypes.end(),  {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"});
	std::map<string,int> GL;
//	cout << GL.size() << endl;
	for (i=0; i<10; i++){
		getline (linestream, w, '\t');
//		cout << genotypes[i] << ":" << atof(w.c_str()) << "\t";
		if (w!="")
			GL.insert(std::make_pair(genotypes[i], atoi(w.c_str()) ));
		else
			continue;
//			GL.insert(std::make_pair(genotypes[i], 0 ));
//		GL[genotypes[i]]= atoi(w.c_str());
	}
//	cout << endl;
//	for (auto pos=GL.begin(); pos!=GL.end(); pos++){
//		cout << pos->first << ":" << pos->second << "\t";
//	}
//	cout << endl;
//	return GL;
//	cout << GL.size() << endl;
	GLs = GL;
}

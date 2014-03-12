#include "genoped.h"


void pedigree_reader(string pedfilename){
    
    std::ifstream inFile(fname);
    long numlines =  std::count(std::istreambuf_iterator<char>(inFile), std::istreambuf_iterator<char>(), '\n');  
    ifstream file (pedfilename, ios::in|ios::ate);
    std::unordered_map<std::string, std::string> pedigree_str;
    if (file.is_open()){
    	file.seekg (0, ios::beg);
        cout << "file is open" << endl;
        while(1)
        {
            string line;
            if (getline(file, line)){
                if (line.substr(0,1) != "#"){
                    vector<string> fields = split(line, '\t');
                    std::stringstream ss2;
                    ss2 << fields[3] << "\t" << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[3];
                    std::string s2 = ss2.str();
                    std::stringstream ss1;
                    ss1 << fields[3] << ":" << fields[0];
                    std::string s1 = ss1.str();
                    pedigree_str.insert(std::make_pair(ss1, ss2));
                }
            }
        }
    }
    else
    {
        cout << "Warning: cannot open the file!" << endl;
    }
    file.close();
}


void parse(string& line, bool parseSamples) {

    // clean up potentially variable data structures
    info.clear();
    infoFlags.clear();
    format.clear();
    alt.clear();
    alleles.clear();

    // #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE1 .. SAMPLEN]
    vector<string> fields = split(line, '\t');
    if (fields.size() < 7) {
        cerr << "broken VCF record (less than 7 fields)" << endl
             << line << endl;
        exit(1);
    }

    sequenceName = fields.at(0);
    char* end; // dummy variable for strtoll
    position = strtoll(fields.at(1).c_str(), &end, 10);
    id = fields.at(2);
    ref = fields.at(3);
    alt = split(fields.at(4), ","); // a comma-separated list of alternate alleles

    // make a list of all (ref + alts) alleles, allele[0] = ref, alleles[1:] = alts
    // add the ref allele ([0]), resize for the alt alleles, and then add the alt alleles
    alleles.push_back(ref);
    alleles.resize(alt.size()+1);
    std::copy(alt.begin(), alt.end(), alleles.begin()+1);

    // set up reverse lookup of allele index
    altAlleleIndexes.clear();
    int n = 0;
    for (vector<string>::iterator a = alt.begin();
            a != alt.end(); ++a, ++n) {
        altAlleleIndexes[*a] = n;
    }

    convert(fields.at(5), quality);
    filter = fields.at(6);
    if (fields.size() > 7) {
        vector<string> infofields = split(fields.at(7), ';');
        for (vector<string>::iterator f = infofields.begin(); f != infofields.end(); ++f) {
            if (*f == ".") {
                continue;
            }
            vector<string> kv = split(*f, '=');
            if (kv.size() == 2) {
                split(kv.at(1), ',', info[kv.at(0)]);
            } else if (kv.size() == 1) {
                infoFlags[kv.at(0)] = true;
            }
        }
    }
    // check if we have samples specified
    // and that we are supposed to parse them
    if (parseSamples && fields.size() > 8) {
        format = split(fields.at(8), ':');
        // if the format changed, we have to rebuild the samples
        if (fields.at(8) != lastFormat) {
            samples.clear();
            lastFormat = fields.at(8);
        }
        vector<string>::iterator sampleName = sampleNames.begin();
        vector<string>::iterator sample = fields.begin() + 9;
        for (; sample != fields.end() && sampleName != sampleNames.end(); ++sample, ++sampleName) {
            string& name = *sampleName;
            if (*sample == "." || *sample == "./.") {
                samples.erase(name);
                continue;
            }
            vector<string> samplefields = split(*sample, ':');
            vector<string>::iterator i = samplefields.begin();
            if (samplefields.size() != format.size()) {
                // ignore this case... malformed (or 'null') sample specs are caught above
                // /*
                // cerr << "inconsistent number of fields for sample " << name << endl
                // << "format is " << join(format, ":") << endl
                // << "sample is " << *sample << endl;
                // exit(1);
                // *
            }
            else {
                for (vector<string>::iterator f = format.begin(); f != format.end(); ++f) {
                    samples[name][*f] = split(*i, ','); ++i;
                }
            }
        }
        if (sampleName != sampleNames.end()) {
            cerr << "error: more sample names in header than sample fields" << endl;
            cerr << "samples: " << join(sampleNames, " ") << endl;
            cerr << "line: " << line << endl;
            exit(1);
        }
        if (sample != fields.end()) {
            cerr << "error: more sample fields than samples listed in header" << endl;
            cerr << "samples: " << join(sampleNames, " ") << endl;
            cerr << "line: " << line << endl;
            exit(1);
        }
    } else if (!parseSamples) {
        originalLine = line;
    }

    //return true; // we should be catching exceptions...
}


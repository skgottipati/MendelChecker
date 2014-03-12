#ifndef __SPLIT_H
#define __SPLIT_H

// functions to split a string by a specific delimiter
#include <string>
#include <vector>
#include <sstream>
#include <string.h>


void pedigree_reader(string pedfilename);

void parse(string& line, bool parseSamples);

#endif
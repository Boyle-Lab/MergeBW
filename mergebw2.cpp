//
//  blacklist.cpp
//  blacklist generation for ENCODE
//
//  Author: Alan Boyle
//  Blacklist Copyright (c) 2018 Alan Boyle
//  This program comes with ABSOLUTELY NO WARRANTY.
//  This is free software, and you are welcome to redistribute it
//  under certain conditions.

extern "C" {
    #include "./libBigWig/bigWig.h"
}
#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <stdint.h>
#include <math.h>
#include <dirent.h>
#include <string>
#include <algorithm>
#include <iterator>
using namespace std;



int main(int argc, char* argv[])
{

	std::vector<float> binsInput (10000000, 0);
	string wigFile = "/nfs/boylelabnr_turbo/ENCODE/bigwig/DNase_regulome/ENCFF013IRY.bigWig";
	string chromosome = "chr1";

	char *fname = new char[wigFile.length() + 1];
	strcpy(fname, wigFile.c_str());
	bigWigFile_t *bwFile = bwOpen(fname, NULL, "r");

	char *chrom = new char[chromosome.length() + 1];
	strcpy(chrom, chromosome.c_str());


	if(bwFile == NULL){
		cerr << "Failed to open file: " << wigFile << endl;
		exit(1);
	}

      	bwOverlappingIntervals_t *ptr = bwGetValues(bwFile, chrom,
			static_cast<uint32_t>(1),
			static_cast<uint32_t>(100), 0);

        for(int k = 0; k < (int)(ptr->l); ++k){
		if(!isnan( ptr->value[k] )){
                	binsInput[k] = roundf(ptr->value[k] * 1000.0) / 1000.0;
			cout << binsInput[k];
                }
        }
	cout << binsInput[1000000] << endl;
}

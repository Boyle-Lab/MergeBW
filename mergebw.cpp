//
//  mergebw.cpp
//  tool for merging bigwig files
//
//  Author: Alan Boyle
//  Copyright (c) 2020 Alan Boyle
//  This program comes with ABSOLUTELY NO WARRANTY.
//  This is free software, and you are welcome to redistribute it
//  under certain conditions.

extern "C" {
    #include "./libBigWig/bigWig.h"
}

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <stdint.h>
#include <math.h>
#include <dirent.h>
#include <string>
#include <algorithm>
#include <iterator>
using namespace std;


class SignalData {
	public:
		std::vector<float> binsInput;

		// vector of vectors
		// chrom names and chrom lengths

		SignalData(std::string wigFile, std::string chromosome, int chromLength) {
			this->wigFile = wigFile;
			this->chromosome = chromosome;
			this->chromLength = chromLength;
		}

		SignalData() {
		}

		void setWigFile(std::string wigFile) {
			this->wigFile = wigFile;
		}

		std::string getWigFile() {
			return wigFile;
		}


// Method:
// 1: Read through every poisiton to create a ranked mean - for every cell type sort all positions and create a mean of the values at each rank
//  This rankedMean list can probably stay in memory
// need rankedMean and rankedMeanTie for each bin in each chromosome
// 2: Read through every cell type rankify all bins and replace with appropriate ranked mean/mean tie
//  Either output as Quantile normalized or this can be used to identify max, min, etc at each position

		// Read in the bigWig file for all chromosomes into memory
		void readBigWig()
		{
			int chromStart = 0;
			int chromEnd = chromLength;
			std::vector<float> binsTemp (chromEnd, 0);

			// Note that libBigWig uses char* for input so we have to convert
		        char *fname = new char[wigFile.length() + 1];
		        strcpy(fname, wigFile.c_str());
		        bigWigFile_t *bwFile = bwOpen(fname, NULL, "r");


		        char *chrom = new char[chromosome.length() + 1];
		        strcpy(chrom, chromosome.c_str());

			// If there is no file attached die
		        if(bwFile == NULL){
                		cerr << "Failed to open file: " << wigFile << endl;
		                exit(1);
        		}

//			binsTemp.resize() //# chroms, size of each chrom... no

		        bwOverlappingIntervals_t *ptr = bwGetValues(bwFile, chrom,
                				        static_cast<uint32_t>(chromStart),
                        				static_cast<uint32_t>(chromEnd),
							1);

			for(int k = 0; k < (int)(ptr->l); ++k){
                		if(!isnan( ptr->value[k] )){
                    			binsTemp[k] = roundf(ptr->value[k] * 1000.0) / 1000.0;
		                }
			}
			//cout << bwFile->cl->chrom[0] << " " << bwFile->cl->len[0] << endl;
		}

	private:
		std::string wigFile;
		std::string chromosome;
		int chromLength;

};

int getdir(std::string dirname, std::vector<string> & files, std::string filetype)
{
	DIR *dir;
	int pos;
	struct dirent *ent;
	dir = opendir(dirname.c_str());
	if (dir != NULL) {
	  	while ((ent = readdir (dir)) != NULL) {
			pos = strlen(ent->d_name) - 7;
			if (! strcmp(&ent->d_name[pos], filetype.c_str())) {
				//printf("%s\n", ent->d_name); //DEBUG
				files.push_back(string(ent->d_name));
			}
	  	}
  		closedir (dir);
	} else {
 	 	/* could not open directory */
		cerr << "Unable to read input files!" << endl;
		exit(1);
	}
}

int main(int argc, char* argv[])
{
	std::vector<string> fileList;
	std::vector<SignalData> inputData;
	std::string wigFile;
	std::string chromosome;
	int chromLength;


	// Will need to read in a file of chromosome sizes and names
	chromosome = "chr1";
	chromLength = 249250621;

	//Parameters
	std::string inputFolder = "/nfs/boylelabnr_turbo/ENCODE/bigwig/DNase_regulome/"; // make this a commandline option - ends with /!
	std::string tempFolder = ".";

	//Read in bigWig files for processing
	getdir(inputFolder, fileList, ".bigWig"); // may need to also consider bw files

	for(int i = 0; i < fileList.size(); i++) {
		wigFile = inputFolder + fileList[i];
		inputData.push_back(SignalData(wigFile, chromosome, chromLength));
		cerr << "Processing " << wigFile << endl;
		inputData.back().readBigWig();
	}

	//Options for output: raw or quantile normalized: max, min, quantile
	//Quantile Normalize

}

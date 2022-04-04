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
#include <tuple>
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
			this->binsInput.clear();

			// Note that libBigWig uses char* for input so we have to convert
		        char *fname = new char[wigFile.length() + 1];
		        strcpy(fname, wigFile.c_str());
		        bigWigFile_t *bwFile = bwOpen(fname, NULL, "r");

			// convert string to char*
		        char *chrom = new char[chromosome.length() + 1];
		        strcpy(chrom, chromosome.c_str());

			// If there is no file attached die
		        if(bwFile == NULL){
                		cerr << "Failed to open file: " << wigFile << endl;
		                exit(1);
        		}

		        bwOverlappingIntervals_t *ptr = bwGetValues(bwFile, chrom,
                				        static_cast<uint32_t>(chromStart),
                        				static_cast<uint32_t>(chromEnd),
							0);

			for(int k = 0; k < (int)(ptr->l); ++k){
                		if(!isnan( ptr->value[k] )){
                    			this->binsInput.push_back(roundf(ptr->value[k] * 1000.0) / 1000.0);
		                }
			}
//			cout << bwFile->cl->chrom[0] << " " << bwFile->cl->len[0] << endl;
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

double quantile(std::vector<double> v, double q) {
	sort(v.begin(), v.end());
	double h = ((v.size() - 1) * q) + 1;
	if(h >= v.size()) {
		h = v.size() - 1;
	}
	return v[floor(h)] + ((h - floor(h)) * (v[floor(h) + 1] - v[floor(h)]));
}

int quantile(std::vector<int> v, double q) {
	sort(v.begin(), v.end());
	double h = (((double)v.size() - 1) * q) + 1;
	return floor((double)v[floor(h)] + ((h - floor(h)) * ((double)v[floor(h) + 1] - (double)v[floor(h)])));
}

// Function to find rank
void rankify(std::vector<double>& A) {
	int n = A.size();
	std::vector<double> R(n,0);
	std::vector< std::tuple<double, int> > T;
	int r = 1;

	// Create array of tuples storing value and index
	for(int j = 0; j < n; j++) {
		T.push_back(std::make_tuple(A[j], j));
	}

	// Sort tubples by data value
	std::sort(begin(T), end(T), [](auto const &t1, auto const &t2) {
        	return get<0>(t1) < get<0>(t2); // or use a custom compare function
	});

	int i = 0;
	int index, j;
	while(i < n) {
		j = i;

		// Get elements of same rank
		while(j < n && std::get<0>(T[j]) == std::get<0>(T[j+1])) {
			j++;
		}

		int m = j - i + 1;

		for(j = 0; j < m; j++) {
			// For each equal element use .5
			index = std::get<1>(T[i+j]);
			R[index] = r + (m-1)*0.5;
		}

		// Increment rank and index
		r+=m;
		i+=m;
	}

	A.swap(R);
}

// convert vectors to vector of SignalData
void quantileNormalize(std::vector<SignalData>& data) {
	int cellCount = data.size();
	int binCount = data[0].binsInput.size();

	//First calculate rank means
	std::vector<double> rankedMean(binCount,0);
	for(int cellID = 0; cellID < cellCount; cellID++) {
		std::vector<double> x(binCount,0);
		for(int binID = 0; binID < binCount; binID++) {
			x[binID] = data[cellID].binsInput[binID];
		}

		sort(x.begin(), x.end());

		for(int binID = 0; binID < binCount; binID++) {
			rankedMean[binID] += x[binID];
		}
	}
	for(int binID = 0; binID < binCount; binID++) {
		rankedMean[binID] /= (double)cellCount;
	}

	//calculate half value for ties
	std::vector<double> rankedMeanTie(binCount-1,0);
	for(int binID = 0; binID < (binCount-1); binID++) {
		rankedMeanTie[binID] = ((rankedMean[binID]+rankedMean[binID+1])/2);
	}

	//Iterate through each cell line
	for(int s = 0; s < cellCount; s++) {
		std::vector<double> bins(binCount,0);
		for(int p = 0; p < binCount; p++) {
			bins[p] = data[s].binsInput[p];
		}
		rankify(bins);

		std::vector<double> binsQuantileNormalized(binCount, 0);
		for(int p = 0; p < binCount; p++) {
			if(std::fmod(bins[p],1) != 0) {
				binsQuantileNormalized[p] = rankedMeanTie[(int)floor(bins[p])-1];
			} else {
				binsQuantileNormalized[p] = rankedMean[(int)(bins[p]-1)];
			}

			data[s].binsInput[p] = binsQuantileNormalized[p];
		}
	}

}

int main(int argc, char* argv[])
{
	std::vector<string> fileList;
	std::vector<SignalData> inputData;
	std::string wigFile;
	std::string chromosome;
	int chromLength;
	int segmentSize;
	std::ifstream chromStream;
	std::string chromName;
	std::string fileName;
	int numChroms  = 0;
	std::vector<int> chromLocation; 

	if (argc < 2) {
		cerr << "Not enough arguments" << endl;
		return 0;
	}
	int readIndex; 
	std::string inputFolder = argv[1];
	std::string tempFolder = "."; 
	std::vector<string> chromList;
	fileName = argv[2];
	int chromNums;

	cout << "Input segment size: " << endl;
	cin >> segmentSize;
	cout << "Input number of chromosomes: ";
	cin >> chromNums;
	cout << "Input starting chromosome index: ";
	cin >> readIndex;

	chromStream.open(inputFolder + fileName);
	if (!chromStream.is_open()) {
		cerr << "Error opening text file" << endl;
		exit(1);
	}

// read in file that has chromosome names and lengths
	while (chromStream >> chromName >> chromLength) {
		numChroms++;
		chromLocation.push_back(chromLength);
		chromList.push_back(chromName);
	}
	chromStream.close(); 

	//Read in bigWig files for processing
	getdir(inputFolder, fileList, ".bigWig"); // may need to also consider bw files

	// Now we need to do this for every chromosome

	// And then for every sub-segment - need to be careful at the end
	std::vector<vector<SignalData>> inputDataVector;
	// For each file perform calculation on the chromosome and segment
	//  we should be able to keep the results in memory and write it all out at the end
	for (int j = readIndex; j < readIndex + chromNums; j++) {
		chromosome = chromList.at(j);
		// read in by segment for each file
		for (int k = 0; k < chromLocation.at(j); k += segmentSize) {
			int l = 0; 
			for (int i = 0; i < fileList.size(); i++) {
				wigFile = inputFolder + fileList[i];
				inputData.push_back(SignalData(wigFile, chromosome, segmentSize + k));
				inputData.back().readBigWig();
			}
			inputDataVector.push_back(inputData);
			inputData.clear();
		}
	}

	//Options for output: raw or quantile normalized: max, min, quantile

	// now quantileNormalize for each segment of each chrom of each file
	// each element of the vector consists of a vector of inputData for each segment for each file
	for (int i = 0; i < inputDataVector.size(); i++) {
		quantileNormalize(inputDataVector[i]);
	}
	
	// do output by segment as well
	// raw output

}

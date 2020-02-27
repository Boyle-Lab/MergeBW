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

		SignalData(std::string wigFile) {
			this->wigFile = wigFile;
		}

		SignalData() {
		}

		void setWigFile(std::string wigFile) {
			this->wigFile = wigFile;
		}

		std::string getWigFile() {
			return wigFile;
		}


		void getSignalBins(std::string chromosome, int chromStart, int chromEnd)
		{
			this->binsInput.clear();
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

			// Read in entire chromosome
		        bwOverlappingIntervals_t *ptr = bwGetValues(bwFile, chrom,
                				        static_cast<uint32_t>(chromStart),
                        				static_cast<uint32_t>(chromEnd),
							1);

			for(int k = 0; k < (int)(ptr->l); ++k){
                		if(!isnan( ptr->value[k] )){
                    			binsTemp[k] = roundf(ptr->value[k] * 1000.0) / 1000.0;
		                }
			}
		}

	private:
		std::string wigFile;


};

int getdir(std::string dirname, std::vector<string> & files, std::string filetype)
{
	DIR *dir;
	int pos;
	struct dirent *ent;
	dir = opendir(dirname.c_str());
	if (dir != NULL) {
	  	while ((ent = readdir (dir)) != NULL) {
			pos = strlen(ent->d_name) - 4;
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
	std::vector<std::tuple<double, int>> T;
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

void quantileNormalize(std::vector<std::vector<double>>& data) {
	int cellCount = data.size();
	int binCount = data[0].size();

	//First calculate rank means
	std::vector<double> rankedMean(binCount,0);
	for(int cellID = 0; cellID < cellCount; cellID++) {
		std::vector<double> x(binCount,0);
		for(int binID = 0; binID < binCount; binID++) {
			x[binID] = data[cellID][binID];
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
			bins[p] = data[s][p];
		}
		rankify(bins);

		std::vector<double> binsQuantileNormalized(binCount, 0);
		for(int p = 0; p < binCount; p++) {
			if(std::fmod(bins[p],1) != 0) {
				binsQuantileNormalized[p] = rankedMeanTie[(int)floor(bins[p])-1];
			} else {
				binsQuantileNormalized[p] = rankedMean[(int)(bins[p]-1)];
			}

			data[s][p] = binsQuantileNormalized[p];
		}
	}

}



std::vector<double> getAbnormalRegions(std::vector<SequenceData> inputData, std::vector<int> binsMap, int type, bool normalize) {
	// Types:
	//  1 - Read
	//		- normalize = Reads / mapability
	//  2 - Multimapping
	//		- normalize = multimapping reads / total reads
	//  3 - Spike (not used)

	std::vector<double> normForQuantile;
	std::vector<double> result(binsMap.size());
	double normTemp;
	double quantileVal = 0;
	std::vector<std::vector<double>> data(inputData.size(), std::vector<double>(binsMap.size()));

	for(int j = 0; j < inputData.size(); j++) { //process by column

		// Generalize this a little
		if(type == 1) { //reads
			inputData[j].binsTemp = inputData[j].binsInput;
		} else if (type == 2) { //multimapping
			inputData[j].binsTemp = inputData[j].binsMultimapping;
		} else if (type == 3) { //spikes
			inputData[j].binsTemp = inputData[j].binsSpikes;
		}

		for(int i = 0; i < binsMap.size(); i++) { //track all rows in this column
			if(normalize) {
				if (type == 2) { //multimapping
					normTemp = (double)inputData[j].binsTemp[i] / (double)inputData[j].binsInput[i];
				} else {
					normTemp = (double)inputData[j].binsTemp[i] / (double)binsMap[i];
				}
			} else {
				normTemp = (double)inputData[j].binsTemp[i] / (double)inputData[j].totalReads * (double)1000000;
			}

			// Save it if we didn't divide by zero
			if(!isnan(normTemp) && !isinf(normTemp) && normTemp > 0) {
				data[j][i] = normTemp;
			} else {
				data[j][i] = 0.0;
			}
		}
	}

	quantileNormalize(data);

	//Now we collapse rows
	std::vector<double> means(inputData.size());
	for(int i = 0; i < binsMap.size(); i++) { // over each row
		for(int j = 0; j < inputData.size(); j++) { //over each column
			means[j] = data[j][i];
		}
		result[i] = quantile(means, 0.5); // This is median signal
	}

	return result;
}

int main(int argc, char* argv[])
{
	std::vector<string> fileList;
	std::vector<SequenceData> inputData;
	std::string wigFile;
	std::string chromosome;
	int chromLength;

	std::string mappabilityFile;
	std::string bamIndexFile;
	std::string refName;
	int binCount;
	std::vector<double> readNormList;
	std::vector<double> multiList;
/*
	//Parameters
	int binSize = 1; // make this a commandline option
	std::string inputFolder = "input/"; // make this a commandline option - ends with /!

	//Read in bigWig files
	getdir(inputFolder, fileList, ".bigwig"); // may need to also consider bw files
	for(int i = 0; i < fileList.size(); i++) {
		wigFile = inputFolder + inputFileList[i];
		inputData.push_back(SignalData(wigFile));
		cerr << "Processing " << wigFile << endl;
		inputData.back().getSignalBins(binSize, "chr1", 1, 1000000);
	}
*/
        wigFile = "/nfs/boylelabnr_turbo/ENCODE/bigwig/DNase_regulome/ENCFF013IRY.bigWig";
	chromosome = "chr1";
	chromLength = 249250621;

	inputData.push_back(SignalData(wigFile));
	inputData.back().getSignalBins("chr1", 1, chromLength);

	//Options for output: raw or quantile normalized: max, min, quantile
	//Quantile Normalize

}

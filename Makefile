prefix=$(shell pwd)

BAMTOOLS_INCLUDE_DIR=$(prefix)/bamtools/install/include/bamtools
BAMTOOLS_LIB_DIR=$(prefix)/bamtools/install/lib/bamtools

all:
	g++ -std=c++14 -o MergeBW mergebw.cpp -I$(BAMTOOLS_INCLUDE_DIR) -L$(BAMTOOLS_LIB_DIR) -lbamtools -lz -Wl,-rpath,$(BAMTOOLS_LIB_DIR)

debug:
	g++ -std=c++14 -g -o MergeBW mergebw.cpp -I$(BAMTOOLS_INCLUDE_DIR) -L$(BAMTOOLS_LIB_DIR) -lbamtools -lz -Wl,-rpath,$(BAMTOOLS_LIB_DIR)


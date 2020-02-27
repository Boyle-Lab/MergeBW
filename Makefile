prefix=$(shell pwd)

LIBBIGWIG=$(prefix)/libBigWig

all:
	g++ -std=c++14 -o MergeBW mergebw.cpp -L$(LIBBIGWIG) -lBigWig -Wl,-rpath,$(LIBBIGWIG)
#	g++ -std=c++14 -o MergeBW mergebw2.cpp -I$(LIBBIGWIG_INCLUDE_DIR) -L$(LIBBIGWIG_LIB_DIR) -llibbigwig -lz -Wl,-rpath,$(LIBBIGWIG_LIB_DIR)

debug:
	g++ -std=c++14 -g -o MergeBW mergebw.cpp -I$(BAMTOOLS_INCLUDE_DIR) -L$(BAMTOOLS_LIB_DIR) -lbamtools -lz -Wl,-rpath,$(BAMTOOLS_LIB_DIR)


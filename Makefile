prefix=$(shell pwd)

LIBBIGWIG=$(prefix)/libBigWig

all:
#	g++ -std=c++14 -o MergeBW mergebw.cpp -L$(LIBBIGWIG) -lBigWig -Wl,-rpath,$(LIBBIGWIG)
	g++ -std=c++14 -o MergeBW mergebw2.cpp -L$(LIBBIGWIG) -lBigWig -Wl,-rpath,$(LIBBIGWIG)


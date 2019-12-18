all : 
	cd src && $(MAKE) parallel-world CXX=g++ CXXFLAGS="-fopenmp -g -lgsl -lgslcblas -lm"

intel : 
	cd src && $(MAKE) parallel-world CXX=icc CXXFLAGS="-qopenmp -g -std=c++17 -lgsl -lgslcblas -lm"

.PHONY : doc clean

doc : 
	cd src && $(MAKE) doc

clean : 
	rm -f bin/parallel-world
	rm -f bin/generated-files/*
	rm -rf doc/*

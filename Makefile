ifeq (,$(findstring intercept,$(CXX)))
CXX=g++
endif
#ifneq ($(CC),intercept-c)
ifeq (,$(findstring intercept,$(CC)))
CC=gcc
endif
############################################################
#change these directories to reflect your boost installation
BOOSTLIBDIR=/usr/local/lib 
BOOSTHDRDIR=/usr/local/include
############################################################
BOOST_LIB=-L $(BOOSTLIBDIR) -lmpc -lmpfr -lgmp
CXXFLAGS= -Wall -std=c++17 -O3 -I $(BOOSTHDRDIR) 
HEADERS=./quartic.hpp ./pvector.hpp ./rpoly.hpp ./cpoly.hpp
LDFLAGS=-lm -llapack -lblas $(BOOST_LIB) 
all: statanalysis poly_real poly_mp poly_cmplx timingtest

poly_real: poly_real.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o poly_real poly_real.cpp  

poly_cmplx: poly_cmplx.cpp $(HEADERS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o poly_cmplx poly_cmplx.cpp  

poly_mp: poly_mp.cpp $(HEADERS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o poly_mp poly_mp.cpp  

statanalysis: statanalysis.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o statanalysis statanalysis.cpp  

timingtgest: timingtest.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o timingtest timingtest.cpp  

clean:
	rm -f timingtest poly_real poly_cmplx poly_mp statanalysis *.o

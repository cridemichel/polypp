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
CXXFLAGS= -Wall -std=c++17 -O3 -I $(BOOSTHDRDIR) $(LDFLAGS) 
HEADERS=./quartic.hpp ./pvector.hpp 
LDFLAGS=-lm -llapack -lblas $(BOOST_LIB) 
all: quartic quartic_mp quartic_cmplx accuracytest statanalysis

quartic: quartic.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o quartic quartic.cpp  

quartic_mp: quartic_mp.cpp $(HEADERS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o quartic_mp quartic_mp.cpp  

quartic_cmplx: quartic_cmplx.cpp $(HEADERS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o quartic_cmplx quartic_cmplx.cpp  

accuracytest: accuracytest.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o accuracytest accuracytest.cpp  

statanalysis: statanalysis.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o statanalysis statanalysis.cpp  

clean:
	rm -f quartic quartic_mp quartic_cmplx accuracytest statanalysis *.o

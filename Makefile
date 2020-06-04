ifeq (,$(findstring intercept,$(CXX)))
CXX=g++-9
endif
############################################################
#change this directories to reflect your boost installation
ifeq ($(BOOSTDIR),)
BOOSTDIR=/usr/local/
endif
############################################################
BOOSTLIBDIR=$(BOOSTDIR)/lib 
BOOSTHDRDIR=$(BOOSTDIR)/include # if PARALLEL is set to 1 parallelization through openmp is enabled,
# but you have to use gnu gcc for this.
PARALLEL=1
BOOST_LIB=-L $(BOOSTLIBDIR) -lmpc -lmpfr -lgmp -lgmpxx
CXXFLAGS= -Wall -std=c++17 -O3 -I $(BOOSTHDRDIR) 
HEADERS=quartic.hpp pvector.hpp rpoly.hpp cpoly.hpp
ifeq ($(PARALLEL),1)
  PARFLA=-fopenmp
else
  PARGLA=
endif 
LDFLAGS=-lm -llapack -lblas $(BOOST_LIB) $(PARFLA)

all: statanalysis poly_real poly_mp poly_cmplx timingtest

poly_real: poly_real.cpp $(HEADERS)
	$(CXX) poly_real.cpp $(CXXFLAGS) -o poly_real 

poly_cmplx: poly_cmplx.cpp $(HEADERS) 
	$(CXX) poly_cmplx.cpp $(CXXFLAGS) $(LDFLAGS) -o poly_cmplx   

poly_mp: poly_mp.cpp $(HEADERS) 
	$(CXX) poly_mp.cpp $(CXXFLAGS) $(LDFLAGS) -o poly_mp   

statanalysis: statanalysis.cpp $(HEADERS)
	$(CXX) statanalysis.cpp $(CXXFLAGS) $(LDFLAGS) -o statanalysis 

timingtgest: timingtest.cpp $(HEADERS)
	$(CXX) timingtest.cpp $(CXXFLAGS) $(LDFLAGS) -o timingtest 

clean:
	rm -f timingtest poly_real poly_cmplx poly_mp statanalysis *.o

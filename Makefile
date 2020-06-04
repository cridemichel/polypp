ifeq (,$(findstring intercept,$(CXX)))
CXX=g++-9
endif
############################################################
#change these directories to reflect your boost installation
BOOSTLIBDIR=/usr/local/lib 
BOOSTHDRDIR=/usr/local/include
############################################################
# if PARALLEL is set to 1 parallelization through openmp is enabled,
# but you have to use gnu gcc for this.
PARALLEL=0
BOOST_LIB=-L $(BOOSTLIBDIR) -lmpc -lmpfr -lgmp
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

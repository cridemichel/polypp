# check if homebrew is installed 
ifeq ($(shell which brew),)
   $(error Please install homebrew first!)
endif

HBPACK=$(shell brew ls --versions gmp gcc@9 boost) 
ifeq ($(shell echo $(HBPACK)),)
   $(warning Please install gmp, boost and gcc@9 homebrew packages through the command:)
   $(warning > brew install gmp boost gcc@9)
   $(error Aborting...)
endif
ifeq ($(shell echo $(HBPACK)|grep gmp),)
   $(warning Please install gmp homebrew packages through the command:)
   $(warning > brew install gmp)
   $(error Aborting...)
endif
ifeq ($(shell echo $(HBPACK)|grep gcc | grep 9),)
   $(warning Please install gcc homebrew packages through the command:)
   $(warning > brew install gcc@9)
   $(error Aborting...)
endif
ifeq ($(shell echo $(HBPACK)|grep boost),)
   $(warning Please install boost homebrew packages through the command:)
   $(warning > brew install boost)
   $(error Aborting...)
endif

ifeq ($(HBDIR),)
HBDIR=$(shell brew --prefix)
endif
ifeq (,$(findstring intercept,$(CXX)))
  CXXHB=$(HBDIR)/bin/g++-9
  #check if g++-9 exists
  ifneq ("$(wildcard $(HBDIR))","")
    CXX=$(CXXHB)
  else
    CXX=g++
  endif
endif
HBLIBDIR=$(HBDIR)/lib 
HBHDRDIR=$(HBDIR)/include # if PARALLEL is set to 1 parallelization through openmp is enabled,
# but you have to use gnu gcc for this.
PARALLEL=1
LIBS=-L $(HBLIBDIR) -lmpc -lmpfr -lgmp -lgmpxx
CXXFLAGS= -Wall -std=c++17 -O3 -I $(HBHDRDIR) 
HEADERS=quartic.hpp pvector.hpp rpoly.hpp cpoly.hpp
ifeq ($(PARALLEL),1)
  PARFLA=-fopenmp
else
  PARGLA=
endif 
LDFLAGS=-lm $(LIBS) $(PARFLA)

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

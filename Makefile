BREWINST=.brew_packages_installed_
ifneq ($(MAKECMDGOALS),poly_real)
ifeq ("$(wildcard $(BREWINST))","")
ifeq ($(shell command -v brew 2> /dev/null),)
  $(error Please install homebrew first!)
endif
# check if homebrew is installed 
HBPACK=$(shell brew ls --versions gmp gcc boost) 
ifeq ($(shell echo $(HBPACK)),)
  $(warning Please install gmp, boost and gcc homebrew packages through the command:)
  $(warning > brew install gmp boost gcc)
  $(error Aborting...)
endif
ifeq ($(shell echo $(HBPACK)|grep gmp),)
  $(warning Please install gmp homebrew packages through the command:)
  $(warning > brew install gmp)
  $(error Aborting...)
  GMPPAK=0
else	
  GMPPAK=1
endif
ifeq ($(shell echo $(HBPACK)|grep gcc),)
  $(warning Please install gcc homebrew packages through the command:)
  $(warning > brew install gcc)
  $(error Aborting...)
  GCCPAK=0
else	
  $(shell CC=$[$CC+1])
  GCCPAK=1
endif
ifeq ($(shell echo $(HBPACK)|grep boost),)
  $(warning Please install boost homebrew packages through the command:)
  $(warning > brew install boost)
  $(error Aborting...)
  BOOSTPAK=0
else
  BOOSTPAK=1
endif
ifeq ($(GMPPAK),1)
  ifeq ($(GCCPAK),1)
    ifeq ($(BOOSTPAK),1) 
      $(shell touch $(BREWINST))
    endif
  endif
endif
endif
endif
ifneq ($(shell command -v brew 2> /dev/null),)
ifeq ($(HBDIR),)
  HBDIR=$(shell brew --prefix)
endif
endif
ifeq (,$(findstring intercept,$(CXX)))
  CXXHB=g++-13
  #check if g++-9 exists
  ifneq ("$(wildcard $(HBDIR))","")
    CXX=$(CXXHB)
  else
    CXX=g++
  endif
endif
ifneq ($(HBDIR),)
HBLIBS=-L $(HBDIR)/lib -lmpc -lmpfr -lgmp -lgmpxx
HBHDRS=-I $(HBDIR)/include
else
HBLIBS=
HBHDRS=
endif
# if PARALLEL is set to 1 parallelization through openmp is enabled,
# but you have to use gnu gcc for this.
ifeq ($(PARALLEL),)
PARALLEL=0
endif
LIBS=$(HBLIBS) 
CXXFLAGS= -Wall -std=c++17 -O3 
CXXFLAGSMP=$(CXXFLAGS) $(HBHDRS)
HEADERS=quartic.hpp pvector.hpp rpoly.hpp cpoly.hpp cpolyvp.hpp rpolyvp.hpp
ifeq ($(PARALLEL),1)
  PARFLA=-fopenmp
else
  PARGLA=
endif 
LDFLAGS=-lm $(LIBS) $(PARFLA)

all: statanalysis poly_real poly_mp poly_cmplx timingtest accuracytest timingtest_vp accuracytest_vp

poly_real: poly_real.cpp $(HEADERS)
	$(CXX) poly_real.cpp $(CXXFLAGS) -o poly_real 

poly_cmplx: poly_cmplx.cpp $(HEADERS) 
	$(CXX) poly_cmplx.cpp $(CXXFLAGSMP) $(LDFLAGS) -o poly_cmplx   

poly_mp: poly_mp.cpp $(HEADERS) 
	$(CXX) poly_mp.cpp $(CXXFLAGSMP) $(LDFLAGS) -o poly_mp   

statanalysis: statanalysis.cpp $(HEADERS)
	$(CXX) statanalysis.cpp $(CXXFLAGSMP) $(LDFLAGS) -o statanalysis 

timingtest: timingtest.cpp $(HEADERS)
	$(CXX) timingtest.cpp $(CXXFLAGSMP) $(LDFLAGS) -o timingtest 

accuracytest: accuracytest.cpp $(HEADERS)
	$(CXX) accuracytest.cpp $(CXXFLAGSMP) $(LDFLAGS) -o accuracytest 

timingtest_vp: timingtest_vp.cpp $(HEADERS)
	$(CXX) timingtest_vp.cpp $(CXXFLAGSMP) $(LDFLAGS) -o timingtest_vp 

accuracytest_vp: accuracytest_vp.cpp $(HEADERS)
	$(CXX) accuracytest_vp.cpp $(CXXFLAGSMP) $(LDFLAGS) -o accuracytest_vp 

clean:
	rm -f timingtest_vp accuracytest_vp timingtest poly_real poly_cmplx poly_mp statanalysis accuracytest *.o $(BREWINST)

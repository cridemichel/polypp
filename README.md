**polypp**
===========

This project provides two C++ classes (*rpoly.hpp* and *cpoly.hpp*) for solving complex and real polynomials of any degree in multiprecision. Quadratic and cubic solver are based on Numerical Recipe book [1].
Quartic solvers are implemented according to Ref. [2].
Roots of polynomials of higher degree are found by using Ehrlich–Aberth algorithm as discussed in Refs. [3-5].
Stopping criterion and accurate calculation of correction term in Aberth method 
have been implemented as suggested in [6].

**References**

```bibliography
[1] W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, Numerical Recipes - The Art of Scientific
    Computing (3rd ed.), Cambridge University Press, Cambridge, UK (2007).
[2] A. G. Orellana and C. De Michele, ACM Transactions on Mathematical Software, 46, No. 2, Article 20 (2020),
    doi: https://doi.org/10.1145/3386241.
[3] D. A. Bini, Numerical Algorithms 13, 179-200 (1996).
[4] D. A. Bini and G. Fiorentino, Numerical Algorithms 23, 127–173 (2000).
[5] D. A. Bini and L. Robol, J. Comp. and Appl. Mat. 272, 276-292 (2014).
[6] T. R. Cameron, Numerical Algorithms, 82, 1065–1084 (2019), doi: https://doi.org/10.1007/s11075-018-0641-9
[7] J. L. Aurentz et al., Siam J. Matrix Anal. Appl., 36, 3 942–973 (2015).
```

In addition to header files you will find some .cpp files with examples on how to use this class, a tool for a statistical analysis of solver accuracy (statanalysis) and another one for testing its performance (timingtest).
Multiprecision is implemented through boost multiprecision libraries (https://www.boost.org/doc/libs/1_73_0/libs/multiprecision/doc/html/index.html) and you need to have both boost and gmp (https://gmplib.org/) installed.
Boost and gmp are conveniently provided by *boost* and *gmp* homebrew packages (https://brew.sh/). 
Note that homebrew not only supports Mac OSX but also Linux and Windows (see https://docs.brew.sh/Homebrew-on-Linux).
To build the sources you need a c++ compiler which complies with *c++17*. It is strongly recommended to install g++ from homebrew. The package is called gcc and it will provide the executable g++-9.

The class itself can be used without boost and gmp and it does not require gcc from homebrew. 
This means that the *poly_real* example can be straightforwardly compiled by the command:
```shell
make poly_real
```
without installing anything else.

If you installed boost, gmp and g++ through homebrew (https://brew.sh), the Makefile should work out of the box for compiling all tools
and examples, i.e. by issuing the command:

```bash
make all
```
the following executables will be created: 

**poly_real**: example of usage of rpoly.hpp class for finding the roots of real polynomials.

**poly_mp**:  example of usage of rpoly.hpp to find the roots of real polynomials in multiple precision (using boost).

**poly_cmplx**: example of usage of cpoly.hpp to find the roots of complex polynomials in  multiple precision (using boost).

**statanalysis**: it performs some statistical analyses inspired by the ones which can be found in Ref. [2]. 
The syntax is the following (where '>' is the shell prompt string):

```bash
statanalysis <trials> <output> <sample> <degree>
```

where:

*trials*: it is the number of roots (samples A-E) or coefficients (sample F) to generate

*output*: every <output> trials save the the probability distribution function P(eps_rel) 
	in the file named P_of_eps_rel.dat and the cumulative distribution function F(eps_rel) 
	in the file named F_of_eps_rel.dat. 

*sample*: it is an integer between 0 and 4 which specifies the sample to generate 
	(cf. Table 4 in Ref. [2] with 0={sample A}, 1={sample B}, 2={sample C}, 3={sample D}, 4={sample E} and
	5={sample F}).

*degree*: degree of the polynomials.

**timingtest**:

this program performs a timing test of the polynomial solver. The syntax is the following (where '>' is the shell prompt string):

time timingtest <trials>

where

<trials> is the number of quartic equations to solve. 
	
**accuracytest**: 

this program carries out several accuracy tests which can be found in Ref. [7]. The syntax is the following (where '>' is the shell prompt string):

> accuracytest <test_number>

where

<test_number> is a number between 1 and 18 which identified the test to perform (see the sources to recognize the different tests)

**Parallelization**

Parallel calculation of roots is implemeneted through OpenMP (https://www.openmp.org/). 
To enable the parallelization gnu gcc must be used and in the Makefile replace
```makefile
PARALLEL=0
```
with 
```makefile
PARALLEL=1
```
then issue the commands: 
```shell
> make clean
> make
```
and to execute the timing test in parallel, do:
```shell
> export OMP_NUM_THREADS=4
> time ./timingtest 1
```
where you can play with the number of threads by changing the value of the environment variable OMP_NUM_THREADES 
and the degree of the polynomial used for the timing test by changing the macro NDEG in the file timingtest.cpp. 

For polynomials of degree less than 5 the present solver is much much faster than tools like MPSolve (ver. 3.1.8, https://numpi.dm.unipi.it/software/mpsolve), while it has a comparable efficiency for higher degree polynomials.
Differently from MPSolve, the present algorithm does not perform any cluster analysis to boost the convergence 
in case of multiple roots and you can only set the working precision but not the precision of the roots.

Enjoy!

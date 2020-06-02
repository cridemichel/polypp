**polypp**
===========

C++ class for solving complex and real polynomials of any degree in multiprecision.
Quadratic and cubic solver are based on Numerical Recipe book [1].
Quartic solvers are implemented according to Ref. [2].
Roots of polynomials of higher degree are found by using Aberth algorithm as discussed in Refs. [3-5].
Stopping criterion and accurate calculation of correction term in Aberth method 
have been implemented as suggested in [6].

**References**

```bibliography
[1] W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery. 2007. Numerical Recipes - The Art of Scientific
    Computing (3rd ed.). Cambridge University Press, Cambridge, UK.
[2] A. G. Orellana and C. De Michele, ACM Transactions on Mathematical Software, 46, No. 2, Article 20 (2020),
    doi: https://doi.org/10.1145/3386241.
[3] D. A. Bini, Numerical Algorithms 13, 179-200 (1996).
[4] D. A. Bini and G. Fiorentino, Numerical Algorithms 23, 127–173 (2000).
[5] D. A. Bini et al. Numerical Algorithms 34, 217–227 (2003).
[6] T. R. Cameron, Numerical Algorithms, 82, 1065–1084 (2019), doi: https://doi.org/10.1007/s11075-018-0641-9
```
In addition to header files you will find some .cpp files with examples on how to use this class.
Multiprecision is implemented through boost multiprecision libraries (https://www.boost.org/doc/libs/1_73_0/libs/multiprecision/doc/html/index.html).
Before compiling change the directories related to boost in the Makefile (BOOSTLIBDIR and BOOSTHDRDIR).
If you installed boost through homebrew (https://brew.sh) the Makefile should work out of the box.

By issuing the command:

```bash
make all
```

you obtain the following executables:

**poly_real**: example of usage of rpoly.hpp class for finding the roots of real polynomials.

**poly_mp**:  example of usage of rpoly.hpp to find the roots of real polynomials in multiple precision (using boost).

**poly_cmplx**: example of usage of cpoly.hpp to find the roots of complex polynomials in  multiple precision (using boost).

**statanalysis**: it performs some statistical analyses inspired by the ones which can be found in Ref. [2]. 
The syntax is the following (where '>' is the shell prompt string):

```bash
> statanalysis <trials> <output> <sample> <degree>
```

where:

*trials*: it is the number of roots (samples A-E) or coefficients (sample F) to generate

*output*: every <output> trials save the the probability distribution function P(eps_rel) 
	in the file named P_of_eps_rel.dat and the cumulative distribution function F(eps_rel) 
	in the file named F_of_eps_rel.dat. 

*sample*: it is an integer between 0 and 4 which specifies the sample to generate 
	(see Table 4 in Ref. [2] with 0={sample A}, 1={sample B}, 2={sample C}, 3={sample D}, 4={sample E} and
	5={sample F}).

*degree*: degree of the polynomials.


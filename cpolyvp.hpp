#ifndef _CPOLYVP_
#define _CPOLYVP_
//#define BINI_CONV_CRIT
/* 
 * NOTES:
 *
 * Find_roots uses by default aberth method (implicit deflation based on newton-raphson method),
 * as discussed in [1] and [2]
 *
 * Aberth method can be parallelized without issues
 
*
 * stopping criterion and accurate calculation of correction term in Abrth method have been implemented as suggested in [4]
 * By using the method set_polish(true), polishing by maehly's method is enabled (by default is disabled since it is slower)
 *
 *
 * set_output_prec() set the precision to use for checking convergence
 *
 * This class can be used also using multiprecision types from boost library (e.g. mpfr and mpc )
 *
 * References
 * [1] D. A. Bini, Numerical Algorithms 13, 179-200 (1996).
 * [2] D. A. Bini and G. Fiorentino, Numerical Algorithms 23, 127–173 (2000).
 * [3] D. A. Bini et al. Numerical Algorithms 34, 217–227 (2003). 
 * [4] T. R. Cameron, Numerical Algorithms, 82, 1065–1084 (2019), doi: https://doi.org/10.1007/s11075-018-0641-9 
 * */
#include "./cpoly.hpp"
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<cmath>
#include <algorithm> 
#include <limits>
#include <cstdlib>
#include <vector>
#include <array>
// To enable parallelizazion use gnu g++ and use the flag -fopenmp, i.e.
// g++ -fopenmp ... 
#if defined(_OPENMP)
#include<omp.h>
#endif
#define Sqr(x) ((x)*(x))
//#define BINI_CONV_CRIT // Bini stopping criterion is slightly less accurate and slightly slower than Cameron one.
#define USE_CONVEX_HULL // <- faster (from  T. R. Cameron, Numerical Algorithms, 82, 1065–1084 (2019)
 
#define USE_ABERTH_REAL //<--- faster
#if defined(_OPENMP)
#define USE_ROLD
#endif
using namespace std;

template <class cmplx, class ntype, class dcmplx> 
class cpolyvp_base_dynamic 
{
public:
  
};

template <class cmplx, class ntype> 
class cpolyvp: public numeric_limits<ntype> {
  const ntype pigr=acos(ntype(-1.0));
  const cmplx I = cmplx(0.0,1.0);
  cpoly<cmplx,-1,ntype> pol;
  ntype eps05, meps, maxf, maxf2, maxf3, minf, scalfact, cubic_rescal_fact;
  unsigned input_precision, output_precision; 
  int maxdigits, n;
  ntype goaleps;
  ntype Kconv;
  bool gpolish;
  bool use_dbl_iniguess;
  bool guess_provided, calc_err_bound;
  pvector<cmplx> coeff, roots;
  bool *found;
 
  void set_coeff(pvector<ntype> v)
    {
      set_input_precision(v[0].precision());
      for (int i=0; i <= n; i++)
        coeff[i].assign(cmplx(v[i],0.0),input_precision);
    }

 void deallocate(void)
    {
      coeff.deallocate();
      delete[] found;
    }
  void allocate(int nc)
    {
      n=nc;
      coeff.allocate(n+1);
      found = new bool[n];
    }
  vector<cmplx> rg;
  vector<int> k;
  struct scoped_precision
   {
      unsigned p;
      scoped_precision(unsigned new_p) : p(ntype::default_precision())
      {
         ntype::default_precision(new_p);
         cmplx::default_precision(new_p);
      }
      ~scoped_precision()
      {
         ntype::default_precision(p);
         cmplx::default_precision(p);
      }
   };
  void set_precision(unsigned p)
    {
      ntype::default_precision(p);
      cmplx::default_precision(p);
    }

public:
  void set_coeff(pvector<cmplx>& v)
    {
      set_input_precision(v[0].precision());
       for (int i=0; i <= n; i++)
        coeff[i].assign(v[i],input_precision);
    }
  void set_input_precision(unsigned p)
    {
      input_precision=p;
    }
  
  int degree()
    {
      return n; 
    }
  
  // find roots by default uses aberth method which is faster than laguerre implicit method
  void find_roots(pvector<cmplx>& roots)
    {
      int prec=output_precision+5;
      set_precision(prec);
      pvector<cmplx> cvp(n+1);
      cpoly<cmplx,-1,ntype,complex<long double>, long double> rp(n);
      pvector<cmplx> roini(n);
      ntype errb, maxerr=0, EPS=pow(ntype(2.0),-ntype(output_precision)*log(10.0)/log(2.0));
      int j;
      for (j=0; j < n+1; j++)
        cvp[j].assign(coeff[j], cvp[j].precision());
      rp.set_coeff(cvp);
      rp.find_roots(roini);
      int nf=0;
      //cout << setprecision(200) << "EPS=" << EPS << "\n";
      for (j=0; j < n; j++)
        {
          errb=rp.calcerrb(roini[j]);
          if (j==0 || errb > maxerr)
            maxerr = errb;
          if (errb <= EPS*abs(roini[j]))
            {
              nf++;
              found[j] = true;
            }
          else
            found[j] = false;
        }
      //cout << setprecision(200) << "maxerr= " << maxerr << "\n";
      if (nf == n)
        {
          for (j=0; j < n; j++)
            roots[j].assign(roini[j], output_precision);
          return;
        }
      else
        {
          for (j=0; j < n; j++)
            roots[j].assign(roini[j], roots[j].precision());
        }
      prec *=2;
      for (int ip=0; ip < 8; ip++)
        {
          set_precision(prec);    
          pvector<cmplx> cvp(n+1);
          cpoly<cmplx,-1,ntype,cmplx,ntype> rp(n);
          ntype errb;
          for (j=0; j < n+1; j++)
            cvp[j].assign(coeff[j], cvp[j].precision());

          pvector<cmplx> ro(n);
          for (j=0; j < n; j++)
            ro[j].assign(roots[j], ro[j].precision());
          rp.use_this_guess(ro);
          rp.set_prec_reached(found);
          rp.set_coeff(cvp);
          rp.find_roots(ro);
          //rp.aberth(roots);
          ntype maxerr;
          nf=0;
          for (j=0; j < n; j++)
            {
              if (found[j])
                {
                  nf++;
                  continue;
                }
              errb=rp.calcerrb(ro[j]);
              if (j==0 || errb > maxerr)
                maxerr = errb;
              if (errb <= EPS*abs(ro[j]))
                {
                  nf++;
                  found[j] = true;
                }
              else
                found[j] = false;
            }

          //cout << setprecision(200) << "2) maxerr= " << maxerr <<  "nf=" << nf << "\n";
          prec *= 2;
          if (nf == n)
            {
              for (j=0; j < n; j++)
                roots[j].assign(ro[j], output_precision);

              break;
            }
          else
            {
              for (j=0; j < n; j++)
                roots[j].assign(ro[j], roots[j].precision());
            }
        }
    }

  void set_output_precision(unsigned op)
    {
      output_precision = op;
    } 

  void init_const(void)
    {
      //cout << "numeric digits=" << maxdigits << " meps=" << meps << "\n";
      input_precision=200;
    }
  cpolyvp(int nc): coeff(nc+1), roots(nc)
  {
    init_const();
    found = new bool[nc];
    n=nc;
  }
  ~cpolyvp()
    {
      delete[] found;
    }
   cpolyvp() = delete;
};
#endif

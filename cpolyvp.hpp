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
template <class cmplx, class ntype>
class azero
{
public:
  cmplx z;
  ntype r;
  vector<int> bonds;
};

template <class cmplx, class ntype, class dcmplx=complex<long double>, class dntype=long double> 
class cpolyvp: public numeric_limits<ntype> {
  const ntype pigr=acos(ntype(-1.0));
  const cmplx I = cmplx(0.0,1.0);
  cpoly<cmplx,-1,ntype> pol;
  ntype eps05, meps, maxf, maxf2, maxf3, minf, scalfact, cubic_rescal_fact;
  unsigned input_precision, output_precision; 
  int maxdigits, n;
  unsigned current_precision;
  ntype goaleps;
  ntype Kconv;
  bool gpolish;
  bool use_dbl_iniguess;
  bool guess_provided, calc_err_bound;
  pvector<cmplx> coeff, roots;
  bool *found;
 
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
  unsigned initial_precision;
  double prec_fact;
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
      current_precision=p;
      ntype::default_precision(p);
      cmplx::default_precision(p);
    }
#if 0
  vector<azero<cmplx,ntype>> particles;
  vector<ntype> errbarr;
  void set_particles(pvector<cmplx>& ro)
    {
      particles.resize(n);
      for (int i=0; i < n; i++)
        {
          particles[i].z = ro[i];
          particles[i].r = errbarr[i];
        }
    }
  void find_bonds(void)
    {
      int i, j; 
      for (i=0; i < n; i++)
        particles[i].bonds.clear();
      /* we need to employ linked cell lists to find bonds
       * otherwise can become very slow for large n */
      for (i=0; i < n; i++)
        {
          for (j=i+1; j < n; j++)
            {
              auto sig = particles[i].r + particles[j].r;
              if (abs(particles[i].z - particles[j].z) <= sig)
                {
                  particles[i].bonds.add(j);
                  particles[j].bonds.add(i);
                }
            }
        }
    }

  void find_clusters(pvector<cmplx>& ro)
    {

    }
  void cluster_analysis(pvector<cmplx>& ro)
    {
      set_particles(ro);
      find_bonds();
      find_clusters(ro);
    }
#endif
public:

  void show(void)
    {
      show(NULL);
    }

  void show(const char* str)
    {
      set_precision(input_precision);
      cpoly<cmplx,-1,ntype> rp(n);
      rp.set_show_digits(input_precision);
      rp.set_coeff(coeff);
      rp.show(str);
    } 
  void set_coeff(pvector<ntype>& v)
    {
      set_input_precision(v[0].precision());
       for (int i=0; i <= n; i++)
        coeff[i].assign(cmplx(v[i],0.0),input_precision);
    }
  void set_coeff(pvector<cmplx>& v)
    {
      set_input_precision(v[0].precision());
       for (int i=0; i <= n; i++)
        coeff[i].assign(v[i],input_precision);
    }
  void set_initial_precision(unsigned p)
    {
      initial_precision=p;
    }
  void set_input_precision(unsigned p)
    {
      input_precision=p;
    }
  
  int degree()
    {
      return n; 
    }
  unsigned auto_precision(void)
    {
      return output_precision+12;
    }
  // find roots by default uses aberth method which is faster than laguerre implicit method
  void find_roots(pvector<cmplx>& roots)
    {
      set_precision(output_precision+5);
      ntype errb, maxerr=0, EPS=pow(ntype(2.0),-ntype(output_precision)*log(10.0)/log(2.0));
      //ntype maxrelerr=0, relerr;
      cmplx roo;
      unsigned prec=initial_precision<=0?auto_precision():initial_precision;
      set_precision(prec);
      pvector<dcmplx> cvp(n+1);
      cpoly<dcmplx,-1,dntype,dcmplx,dntype> rp(n);
      pvector<dcmplx> roini(n);
      int j;
      if constexpr (!(is_same<dcmplx,complex<float>>::value && 
                    is_same<dntype,float>::value) &&
                    !(is_same<dcmplx,complex<double>>::value &&
                    is_same<dntype, double>::value) &&
                    !(is_same<dcmplx,complex<long double>>::value &&
                    is_same<dntype,long double>::value))
        {
          cout << "dcmplx and dntype must be either float, double or long double\n";
          exit(1);
        }

      for (j=0; j < n+1; j++)
        cvp[j]=dcmplx(coeff[j]);
      rp.iniguess_slow();
      rp.set_coeff(cvp);
      rp.find_roots(roini);
      int nf=0;
#if 0
      errbarr.resize(n);
      particles.resize(n);
#endif
      //cout << setprecision(200) << "EPS=" << EPS << "\n";
      for (j=0; j < n; j++)
        {
          errb.assign(ntype(rp.calcerrb(roini[j])), errb.precision());
          //errbarr[j].assign(errb,current_precision); 
#if 0
          if (roinid[j]==dcmplx(0,0))
            relerr = errb;
          else
            relerr = errb/ntype(abs(roinid[j]));
#endif
          if (j==0 || errb > maxerr)
            maxerr = errb;
          //if (j==0 || relerr> maxrelerr)
            //maxrelerr = relerr;

          roo.assign(cmplx(roini[j]), roo.precision());
          if (errb <= EPS*abs(roo))
            {
              nf++;
              found[j] = true;
            }
          else
            {
              //cout << "1)root #" << j << " not accurate enough (errb=" << errb << ")\n";
              found[j] = false;
            }
        }
      //cout << setprecision(200) << "maxerr= " << maxerr << "\n";
      if (nf == n)
        {
          for (j=0; j < n; j++)
            roots[j].assign(cmplx(roini[j]), roots[j].precision());
          return;
        }
      else
        {
          for (j=0; j < n; j++)
            roots[j].assign(cmplx(roini[j]), roots[j].precision());
        }
      //prec = (unsigned)(double(prec)*1.1*abs(log10(EPS)/log10(maxrelerr)));
      //cout << "INIPREC=" << prec << "\n";
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
              errb.assign(rp.calcerrb(ro[j]),errb.precision());
              roo.assign(ro[j],roo.precision());
              if (j==0 || errb > maxerr)
                maxerr = errb;
              if (errb <= EPS*abs(roo))
                {
                  nf++;
                  found[j] = true;
                }
              else
                {
                  //cout << "2)root #" << j << " not accurate enough (errb=" << errb << ")\n";
                  found[j] = false;
                }
            }
          //cout << setprecision(200) << "2) maxerr= " << maxerr <<  "nf=" << nf << "\n";
          prec = (int)(double(prec)*prec_fact);
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
          //cout << "newprec2=" << prec << " iter=" << ip << "\n";
        }
    }

  void set_output_precision(unsigned op)
    {
      output_precision = op;
    } 

  void init_const(void)
    {
      //cout << "numeric digits=" << maxdigits << " meps=" << meps << "\n";
      input_precision=16;
      initial_precision=0; // 0 means "auto"
      prec_fact=2.0;
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
   cpolyvp()=delete;
};
#endif

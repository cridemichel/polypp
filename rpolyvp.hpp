#ifndef _RPOLYVP_
#define _RPOLYVP_
#include "./pvector.hpp"
#include "./cpolyvp.hpp"
#include "./quartic.hpp"
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<cmath>
#include <algorithm> 
#include <limits>
#include <cstdlib>
#include <vector>
#include <array>
#define Sqr(x) ((x)*(x))
//#define FAST_MATH
using namespace std;
template <class ntype, class cmplx, class dcmplx=complex<long double>, class dntype=long double> 
class rpolyvp: public numeric_limits<ntype> {
  int n, maxdigits;
  cpolyvp<cmplx,ntype,dcmplx,dntype> polvp;
  pvector<ntype,-1> coeff, cmon;
  unsigned initial_precision;
  double prec_fact;
  unsigned input_precision, output_precision; 

public:
  
  void set_coeff(pvector<ntype> v)
    {
      set_input_precision(v[0].precision());
      for (int i=0; i <= n; i++)
        coeff[i].assign(v[i],input_precision);
      cmon[n].assign(1.0,input_precision);
      for (int i=n-1; i >=0; i--)
        {
          cmon[i].assign(coeff[i]/coeff[n], input_precision);
        }

    }

 void deallocate(void)
    {
      coeff.deallocate();
      cmon.deallocate();
      polvp.deallocate();
    }
  void allocate(int nc)
    {
      n=nc;
      coeff.allocate(n+1);
      cmon.allocate(n+1);
      polvp.allocate(n);
    }

  int degree()
    {
      return n; 
    }
  
  void set_show_digits(int p)
    {
      maxdigits=p;
    }

  void show(const char* str)
    {
      int i;
      if (str!=NULL)
	cout <<  str;
      for (i=n; i >= 0; i--)
	{
	  if (coeff[i] > 0)
    	    {
	      if (i < n)
		cout << "+";
	    }
	  else
	    { 
	      cout << "-";
	    }
	  if (i==0)
	    cout << setprecision(maxdigits) << abs(coeff[i]);
	  else if (i > 0 && abs(coeff[i]) != 1.0)
	    cout << setprecision(maxdigits) << abs(coeff[i])<< "*";
	 
	  if ( i > 1)
	    {
	      cout << "x^" << i;
	    }
	  else if (i==1)
	    cout << "x";
	}
      cout << "\n";
    }
   cmplx evalpoly(cmplx x)
    {
      // evaluate polynomail via Horner's formula 
      cmplx bn=cmplx(0.0);
      for (int i=n; i >= 0; i--)
        {
          bn = cmon[i] + bn*x;
        }
      return bn;
    }
   cmplx evaldpoly(cmplx x)
    {
      // evaluate polynomail via Horner's formula 
      cmplx bn=0.0;
      for (int i=n-1; i >= 0; i--)
        {
          bn = (i+1)*cmon[i+1] + bn*x;
        }
      return bn;
    }
   cmplx evalddpoly(cmplx x)
    {
      // evaluate second derivative of polynomail via Horner's formula 
      cmplx bn=0.0;
      if (n == 1)
        return 0;
      for (int i=n-2; i >= 0; i--)
        {
          bn = (i+2)*(i+1)*cmon[i+2] + bn*x;
        }
      return bn;
    }
   ntype evalpoly(ntype x)
    {
      // evaluate polynomail via Horner's formula 
      ntype bn=0.0;
      for (int i=n; i >= 0; i--)
        {
          bn = cmon[i] + bn*x;
        }
      return bn;
    }

  ntype evaldpoly(ntype x)
    {
      // evaluate first derivative of polynomail via Horner's formula 
      ntype bn=0.0;

      for (int i=n-1; i >= 0; i--)
        {
          bn = (i+1)*cmon[i+1] + bn*x;
        }
      return bn;
    }
  ntype evalddpoly(ntype x)
    {
      // evaluate second derivative of polynomail via Horner's formula 
      ntype bn=0.0;
      if (n == 1)
        return 0;
      for (int i=n-2; i >= 0; i--)
        {
          bn = (i+2)*(i+1)*cmon[i+2] + bn*x;
        }
      return bn;
    }
  void set_precision(unsigned p)
    {
      ntype::default_precision(p);
      cmplx::default_precision(p);
    }

  void set_initial_precision(unsigned p)
    {
      initial_precision=p;
    }
  void set_input_precision(unsigned p)
    {
      input_precision=p;
    }

  void set_output_precision(unsigned op)
    {
      output_precision = op;
    } 

  void find_roots(pvector<cmplx>& roots)
    {
      polvp.set_input_precision(input_precision);
      if (initial_precision <= 0)
        initial_precision=polvp.auto_precision();
      polvp.set_initial_precision(initial_precision);
      polvp.set_output_precision(output_precision);
      polvp.set_coeff(coeff);
      polvp.find_roots(roots);
    }     
   // get machine precision for "ntype" type (ntype can float, double, long double)
  ntype epsilon()
    {
      return numeric_limits<ntype>::epsilon(); 
    }
  ntype getmax()
    {
      return numeric_limits<ntype>::max();
    }
   void init_const(void)
    {
      input_precision=16;
      initial_precision=0; // 0 means "auto"
      prec_fact=2.0;
    }
  ~rpolyvp()
    {
    }
  rpolyvp() = delete;

  rpolyvp(int nc): polvp(nc)
    {
      n=nc;
      cmon.allocate(nc+1);
      coeff.allocate(nc+1);
      init_const(); 
    }
};
#endif

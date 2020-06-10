#include <stdlib.h>
#include <stdio.h>
#include "./cpolyvp.hpp"
#define WP 0
#define WPD 20
#define WPBS 200
#define MPC_MP
#ifdef CPP_MP
#include <boost/multiprecision/cpp_bin_float.hpp> 
#include <boost/multiprecision/cpp_complex.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using numty = number<cpp_bin_float<WP>>;
using cmplx = cpp_complex<WP>;
using dntype= number<cpp_bin_float<WPD>>;
using dcmplx= cpp_complex<WPD>;
#elif defined(GMP_MP)
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using numty=mpf_float;
using cmplx=number<complex_adaptor<gmp_float<0>>>;
using dntype=number<gmp_float<WPD>>;
using dcmplx=number<complex_adaptor<gmp_float<WPD>>>;
#ifdef BACKSTAB
using bsdbl=number<gmp_float<WPBS>>;
using bscmplx=number<complex_adaptor<gmp_float<WPBS>>>;
#endif
#elif defined(MPC_MP)
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using boost::multiprecision::mpfr_float;
using numty=mpfr_float;
using cmplx=mpc_complex;
using dntype=mpfr_float;
using dcmplx=mpc_complex;
#ifdef BACKSTAB
using bsdbl=number<mpfr_float_backend<WPBS>>;
using bscmplx=number<mpc_complex_backend<BS>>;
#endif
#else
#ifdef BACKSTAB
using bsdbl=number<mpfr_float_backend<WPBS>>;
using bscmplx=number<mpc_complex_backend<WPBS>>;
#endif
using numty=double;
using cmplx=complex<numty>;
using dntype=numty;
using dcmplx=cmplx;
#endif
//#include<complex>
#ifndef NDEG
#define NDEG 1000
#endif

double gauss(void)
{
  double  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  double sum, r, r2;
  int i;

  sum = 0.0;

  for(i=0; i < 12; i++)
    {
      sum = sum + drand48();
    }
  
  r  = ( sum - 6.0 ) / 4.0;
  r2 = r * r;

  return  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r;

}

using namespace std;
int main(int argc, char* argv[])
{
  numty::default_precision(100);
  cmplx::default_precision(100);
  pvector<cmplx> c(NDEG+1);
  pvector<cmplx> roots(NDEG);
  cpolyvp<cmplx,numty> rp(NDEG);
  int j, maxiter;
#ifdef _OPENMP
  cout << "# thread=" << omp_get_max_threads() << "\n";
#endif
  srand48(4242);
  numty sig=1.0;
  if (argc>=2)
    {
      maxiter = atoi(argv[1]);
    }
  else
    maxiter = 1000000;
  c[NDEG]=1.0;
  rp.set_output_precision(100); 

  /* initial precision should be around input_precision (which is automatically
   * set according to the precision of the coefficients) plus 15 for optimale performance */
  //rp.set_initial_precision(115); // if initial precision is not provided, it is automatically estimated

  for (int i=0; i < maxiter; i++)
    {
      //cout << "iter #" << i << "\n";

      for (j=0; j < NDEG; j++)
        c[j]=cmplx(2.0*sig*(drand48()-0.5),0.0);
      rp.set_coeff(c);
      rp.find_roots(roots);
    }

  return 0;
}

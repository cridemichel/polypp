#include <stdlib.h>
#include <stdio.h>
#include "./cpolyvp.hpp"
#define WP 0
#define WPD 20
#define WPBS 200
#define MPC_MP
//#define GMP_MP
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
using numty=number<gmp_float<WP>>;
using cmplx=number<complex_adaptor<gmp_float<WP>>>;
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
using boost::multiprecision::mpfr_float;
int main(int argc, char* argv[])
{
  mpfr_float::default_precision(200);
  mpc_complex::default_precision(200);
  pvector<cmplx> c(NDEG+1);
  pvector<cmplx> roots(NDEG);
  cpolyvp<cmplx,numty> rp(NDEG);
  int j, maxiter;
#ifdef _OPENMP
  cout << "# thread=" << omp_get_max_threads() << "\n";
#endif
  srand48(time(0));
  numty sig=1.0;
  if (argc>=2)
    {
      maxiter = atoi(argv[1]);
    }
  else
    maxiter = 1000000;
  c[NDEG]=1.0;
  rp.set_output_precision(32); 
  for (int i=0; i < maxiter; i++)
    {
      //cout << "iter #" << i << "\n";

      for (j=0; j < NDEG; j++)
        c[j]=cmplx(2.0*sig*(drand48()-0.5),0.0);
      rp.set_coeff(c);
      rp.find_roots(roots);
#if 0
      mpfr_float::default_precision(prec);
      mpc_complex::default_precision(prec);
      pvector<cmplx> cvp(NDEG+1);
      cpoly<cmplx,-1,numty> rp(NDEG);
      pvector<cmplx,-1> roini(NDEG);


      for (j=0; j < NDEG+1; j++)
        cvp[j].assign(c[j], cvp[j].precision());
      rp.set_coeff(cvp);
      rp.find_roots(roini);
      for (j=0; j < NDEG; j++)
        roots[j].assign(roini[j], roots[j].precision());
      cout << "precision=" << roots[0].precision() << "\n";
      prec *=2;
      //roots.show("roots=");
      for (int ip=0; ip < 2; ip++)
        {
          cout << "prec=" << prec << "\n";
          mpfr_float::default_precision(prec);
          mpc_complex::default_precision(prec);
          pvector<cmplx> cvp(NDEG+1);
          cpoly<cmplx,-1,numty,cmplx,numty> rp(NDEG);
          
          for (j=0; j < NDEG+1; j++)
            cvp[j].assign(c[j], cvp[j].precision());
     
          pvector<cmplx,-1> ro(NDEG);
          for (j=0; j < NDEG; j++)
            ro[j].assign(roots[j], ro[j].precision());
          rp.use_this_guess(ro);
          rp.set_coeff(cvp);
          rp.find_roots(ro);
          //rp.aberth(roots);

          prec *= 2;
          for (j=0; j < NDEG; j++)
            roots[j].assign(ro[j], roots[j].precision());
#if 0
          cout << "roprec=" << ro[0].precision() << "\n";
          cout << "cpvprec="<< cvp[0].precision() << "\n"; 
          cout << "rootsprec=" << roots[0].precision() << "\n";
          cout << "coeffprec=" << c[0].precision() << "\n";
#endif
        }
#endif
    }

  return 0;
}

#include "./rpoly.hpp"
#include<stdlib.h>
#include<stdio.h>
#include "./rpoly.hpp"
#include <complex>
#define WP 50
// N.B. you can use either CPP, GMP or MPC backend by
// defining CPP_MP, GMP_MP or MPC_MP
#define NDEG 4
#define MPC_MP
#ifdef CPP_MP
#include <boost/multiprecision/cpp_bin_float.hpp> 
#include <boost/multiprecision/cpp_complex.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using mpreal = number<cpp_bin_float<WP>>;
using mpcmplx = cpp_complex<WP>;
#elif defined(GMP_MP)
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using mpreal=number<gmp_float<WP>>;
using mpcmplx=number<complex_adaptor<gmp_float<WP>>>;
#elif defined(MPC_MP)
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
//we set 100 digits working precision!
using mpreal=number<mpfr_float_backend<WP>>;
using mpcmplx=number<mpc_complex_backend<WP>>;
#endif
int main(void)
{
  rpoly<mpreal,NDEG,mpcmplx> P;
  pvector<mpreal,NDEG+1> c;
  pvector<mpcmplx,NDEG> r;

  for (int i=0; i <= NDEG; i++)
    c[i] = mpreal("1.0");
  P.set_coeff(c);
  P.find_roots(r);
  //r.show("roots");
  int cc=0;
  for (auto& r0: r)
    {
      cout << setprecision(WP) << "root #" << cc <<  "=" << r0 << "\n";
      cout << setprecision(WP) << "p(#" << cc << ")=" << P.evalpoly(r0) << "\n\n";
      cc++;
    }
  return 0;
}

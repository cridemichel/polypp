#include <stdlib.h>
#include <stdio.h>
//#define CPOLY
//#define STATIC
#ifdef CPOLY
#include "./cpoly.hpp"
#else
#include "./rpoly.hpp"
#endif
#define WP 100
#define WPD 20
#define WPBS 200
//#define MPC_MP
#define GMP_MP
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
using numty=number<mpfr_float_backend<WP>>;
using cmplx=number<mpc_complex_backend<WP>>;
using dntype=number<mpfr_float_backend<WPD>>;
using dcmplx=number<mpc_complex_backend<WPD>>;
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

#ifdef BACKSTAB
#ifdef CPOLY
void calc_coeff(bscmplx co[], bscmplx er[])
{
  bscmplx c[NDEG+1], alpha;

  int ii, jj;

  for (ii=0; ii < NDEG; ii++)
    {
      c[ii]  = 0.0;
    }
  c[NDEG]=1.0;
  ii=0;
  
  while (ii < NDEG)
    { 
      alpha = -er[ii];
      for (jj=ii; jj >= 0; jj--)
        {         
          //do jj=ii,1,-1
          if (jj==0)
            c[jj] = c[jj] + alpha;
          else
            c[jj] = c[jj] + alpha*c[jj-1];
        }
      ii=ii+1;
    }
  for (ii=0; ii < NDEG; ii++)
     co[ii] = c[NDEG-ii-1];
  co[NDEG]=1.0;
}
#else
void calc_coeff(bsdbl co[], bscmplx er[])
{
  bscmplx c[NDEG+1], alpha;

  int ii, jj;

  for (ii=0; ii < NDEG; ii++)
    {
      c[ii]  = 0.0;
    }
  c[NDEG]=1.0;
  ii=0;
  
  while (ii < NDEG)
    { 
      alpha = -er[ii];
      for (jj=ii; jj >= 0; jj--)
        {         
          //do jj=ii,1,-1
          if (jj==0)
            c[jj] = c[jj] + alpha;
          else
            c[jj] = c[jj] + alpha*c[jj-1];
        }
      ii=ii+1;
    }
  for (ii=0; ii < NDEG; ii++)
     co[ii] = real(c[NDEG-ii-1]);
  co[NDEG]=1.0;
}
#endif
#ifdef STATIC
#define NN NDEG
#else
#define NN -1
#endif
#ifdef CPOLY
bsdbl calc_backward_err(pvector<cmplx,NN>& roots,  pvector<cmplx,NN>& c)
{
  int i;
  vector<cmplx> rv;
  bscmplx *er = new bscmplx[NDEG];
  bscmplx *cbs = new bscmplx[NDEG+1];
  rv.resize(NDEG);
  for (i=0; i < NDEG; i++)
    {
      rv[i] = cmplx(roots[i]);
    }
  // sort is needed to obtain the right backward error
  std::sort(rv.begin(),rv.end(),[&](cmplx a, cmplx b)->bool {return abs(a) < abs(b);}); 

  for (i=0; i < NDEG; i++)
    {
      er[i] = bscmplx(rv[i]);
    }  
  calc_coeff(cbs, er);
  bsdbl err, errmax;
  for (i=0; i < NDEG+1; i++)
    {
      if (c[i]==cmplx(0.0))
        err=abs(cbs[i] - bscmplx(c[i]));
      else
        err=abs((cbs[i] - bscmplx(c[i]))/bscmplx(c[i]));
      if (err > 100)
        cout << "i=" << i << " cbs=" << cbs[i] << " c=" << c[i] << "\n";
      if (i==0 || err > errmax)
        errmax=err;
   }
  return errmax;
}
#else
bsdbl calc_backward_err(pvector<cmplx,NN> roots,  pvector<numty,NN> c)
{
  int i;
  bscmplx *er = new bscmplx[NDEG];
  bsdbl *cbs = new bsdbl[NDEG+1];
  vector<cmplx> rv;
   
  bsdbl err, errmax;
  rv.resize(NDEG);
  for (i=0; i < NDEG; i++)
    {
      rv[i] = cmplx(roots[i]);
    }
  std::sort(rv.begin(),rv.end(),[&](cmplx a, cmplx b)->bool {return abs(a) < abs(b);}); 
  for (i=0; i < NDEG; i++)
    {
      er[i] = bscmplx(rv[i]); 
    }
  calc_coeff(cbs, er);
  for (i=0; i < NDEG+1; i++)
    {
      if (abs(c[i])==0)
        err=abs(cbs[i] - bsdbl(c[i]));
      else
        err=abs((cbs[i] - bsdbl(c[i]))/bsdbl(c[i]));
      if (i==0 || err > errmax)
        errmax=err;
   }
  return errmax;
}
#endif
#endif
using namespace std;
int main(int argc, char* argv[])
{
#ifdef BACKSTAB
  bsdbl berr, berrmax;
  int poltype, m;
#endif
  int j, maxiter;
#ifdef STATIC
#ifdef CPOLY
  cpoly<cmplx,NDEG,numty> rp;
  pvector<cmplx,NDEG+1> c;
  pvector<cmplx,NDEG> roots;
#else
  rpoly<numty,NDEG,cmplx> rp;
  pvector<numty,NDEG+1> c;
  pvector<cmplx,NDEG> roots;
#endif
#else
#ifdef CPOLY
  cpoly<cmplx,-1,numty> rp(NDEG);
  pvector<cmplx,-1> c(NDEG+1);
  pvector<cmplx,-1> roots(NDEG);
#else
  rpoly<numty,-1,cmplx> rp(NDEG);
  pvector<numty,-1> c(NDEG+1);
  pvector<cmplx,-1> roots(NDEG);
#endif
#endif
#ifdef _OPENMP
  cout << "# thread=" << omp_get_max_threads() << "\n";
#endif
  //srand48(time(0));
  srand48(4242);  
  numty sig=1.0;
  if (argc>=2)
    {
      maxiter = atoi(argv[1]);
    }
  else
    maxiter = 1000000;
#ifdef BACKSTAB
  if (argc >= 3)
    poltype = atoi(argv[2]);
  else
    poltype=0;
#endif
  c[NDEG]=1.0;

#ifdef BACKSTAB
  if (poltype !=0 )
    maxiter=1;
#endif
  for (int i=0; i < maxiter; i++)
    {
      //cout << "iter #" << i << "\n";

#ifdef BACKSTAB
#ifdef CPOLY
      for (j=0; j < NDEG; j++)
        c[j]=cmplx(2.0*sig*(drand48()-0.5),0.0);
#else
      if (poltype==0)
        {
          for (j=0; j < NDEG; j++)
            c[j]=sig*(drand48()-0.5);
        }
      else if (poltype==1)
        {
          if (NDEG%2==1)
            {
              cout << "You have to set an even polynomial degree for this test\n";
              exit(1);
            }
          for (j=0; j < NDEG; j++)
            c[j] = 0.0;
          c[0] = 1.0;
          c[NDEG/2] = NDEG/(NDEG+1.0)+(NDEG+1.0)/NDEG;
          c[NDEG]=1.0;
        }
      else if (poltype==2)
        {
          if (NDEG%2==1)
            {
              cout << "You have to set an even polynomial degree for this test\n";
              exit(1);
            }
          m=NDEG/2;
          c[m]=m+1.0;
          for (j=0; j < m-1; j++)
            c[j] = m+j;
          for (j=2*m; j >= m+1; j--)
            c[j] = 3*m-j;
          for (j=0; j <= 2*m; j++)
            c[j] /= m;
        }
      else if (poltype==3)
        {
          m=NDEG-1;
          numty lambda=0.9;
          for (j=0; j <= NDEG; j++)
            c[j] = 0;
          c[m+1] = 1.0-lambda;
          c[m] = -(lambda+1.0);
          c[1] = lambda+1.0;
          c[0] = -(1.0-lambda); 
          for (j=1; j <= NDEG; j++)
            c[j] /= c[NDEG];
          c[NDEG]=1.0;
        }
      else if (poltype==4)
        {
          m=NDEG-1;
          numty lambda=0.99;
          for (j=0; j <= NDEG; j++)
            c[j] = 0;
          c[m+1] = 1.0-lambda;
          c[m] = -(lambda+1.0);
          c[1] = lambda+1.0;
          c[0] = -(1.0-lambda);
          for (j=1; j <= NDEG; j++)
            c[j] /= c[NDEG];
          c[NDEG]=1.0;
        }
#endif
#else
#ifdef CPOLY
      for (j=0; j < NDEG; j++)
        c[j]=cmplx(2.0*sig*(drand48()-0.5),0.0);
#else

      for (j=0; j < NDEG; j++)
        c[j]=sig*(drand48()-0.5);
#endif
#endif
      rp.set_coeff(c);
      rp.find_roots(roots);
      //rp.aberth(roots);
#ifdef BACKSTAB
      berr=calc_backward_err(roots, c);
      if (i==0 || berr > berrmax)
        berrmax=berr;
#endif
    }
#ifdef BACKSTAB
  cout << "Max backward error=" << berrmax << "\n";
#endif

  return 0;
}

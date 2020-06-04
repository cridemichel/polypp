#ifndef _CPOLY_
#define _CPOLY_
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
#include "pvector.hpp"
#include "quartic.hpp"
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

template <class cmplx, class ntype, class dcmplx, int N> 
class cpoly_base_static 
{
public:
  int n;
  constexpr static int dynamic = false;
  pvector<cmplx, N+1> coeff;
  pvector<ntype, N+1> acmon;
  pvector<cmplx, N+1> cmon;
  pvector<cmplx, N-1> abscmon;
  pvector<ntype, N+1> alpha;
  quartic<ntype,cmplx,false> quar;
  bool found[N];
#ifdef USE_ROLD
  pvector<cmplx, N> rold;
#endif
  pvector<dcmplx,N> droots;
  void set_coeff(pvector<ntype,N+1> v)
    {
      for (int i=0; i <= N; i++)
        coeff[i] = cmplx(v[i],0.0);
      cmon[n]=1.0;
      for (int i=n-1; i >=0; i--)
        {
          cmon[i]=coeff[i]/coeff[n];
        }
    }
 
  void set_coeff(pvector<cmplx,N+1> v)
    {
      coeff = v;
      cmon[n]=1.0;
      for (int i=n-1; i >=0; i--)
        {
          cmon[i]=coeff[i]/coeff[n];
        }
    }
  cpoly_base_static()
    {
      n=N;
    }
  cpoly_base_static(int nc): cpoly_base_static()
  {
    n=nc;
  }
};
template <class cmplx, class ntype, class dcmplx, int N> 
class cpoly_base_dynamic 
{
public:
  int n;
  constexpr static int dynamic = true;
  pvector<cmplx> coeff;
  pvector<cmplx> cmon;
  pvector<ntype> acmon;
  pvector<ntype> alpha;
#ifdef USE_ROLD
  pvector<cmplx> rold;
#endif
  pvector<dcmplx> droots;
  quartic<ntype,cmplx,true> quar;
  bool *found;

  cpoly_base_dynamic() = default;

  void set_coeff(pvector<ntype,-1> v)
    {
      for (int i=0; i <= n; i++)
        coeff[i] = cmplx(v[i],0.0);
      cmon[n]=1.0;
      for (int i=n-1; i >=0; i--)
        {
          cmon[i]=coeff[i]/coeff[n];
        }
    }

  void set_coeff(pvector<cmplx,-1> v)
    {
      coeff = v;
      cmon[n]=1.0;
      for (int i=n-1; i >=0; i--)
        {
          cmon[i]=coeff[i]/coeff[n];
        }
    }
  void use_vec(int nc, cmplx* coeffv, cmplx* cmonv,
               ntype* acmonv, ntype *alphav)
    {
      coeff.use_vec(nc+1,coeffv);
      cmon.use_vec(nc+1,cmonv);
      acmon.use_vec(nc+1,acmonv);
      alpha.use_vec(nc+1,alphav);
      n=nc;
    }

  cpoly_base_dynamic(int nc): coeff(nc+1), cmon(nc+1), acmon(nc+1), alpha(nc+1), droots(nc)
  {
#ifdef USE_ROLD
    rold.allocate(nc);
#endif
    found = new bool[nc];
    n=nc;
  }
  ~cpoly_base_dynamic()
    {
      delete[] found;
    }
  void deallocate(void)
    {
      coeff.deallocate();
      cmon.deallocate();
      acmon.deallocate();
      alpha.deallocate();
#ifdef USE_ROLD
      rold.deallocate();
#endif
      droots.deallocate();
      delete[] found;
    }
  void allocate(int nc)
    {
      n=nc;
      coeff.allocate(n+1);
      cmon.allocate(n+1);
      acmon.allocate(n+1);
      alpha.allocate(n+1);
#ifdef USE_ROLD
      rold.allocate(n);
#endif
      droots.allocate(n);
      found = new bool[n];
    }
};

template <class cmplx, class ntype, class dcmplx, int N> using cpolybase = 
typename std::conditional<(N>0), cpoly_base_static <cmplx, ntype, dcmplx, N>,
         cpoly_base_dynamic <cmplx, ntype, dcmplx, N>>::type;

template <class cmplx, int N=-1, class ntype=double, class dcmplx=complex<long double>, class dntype=long double> 
class cpoly: public numeric_limits<ntype>, public cpolybase<cmplx,ntype,dcmplx,N> 
{
  using cpolybase<cmplx,ntype,dcmplx,N>::n;
  using cpolybase<cmplx,ntype,dcmplx,N>::coeff;
  using cpolybase<cmplx,ntype,dcmplx,N>::cmon;
  using cpolybase<cmplx,ntype,dcmplx,N>::acmon;
  using cpolybase<cmplx,ntype,dcmplx,N>::alpha;
  using cpolybase<cmplx,ntype,dcmplx,N>::droots;
  using cpolybase<cmplx,ntype,dcmplx,N>::quar;
  using cpolybase<cmplx,ntype,dcmplx,N>::found;
  const ntype pigr=acos(ntype(-1.0));
  const cmplx I = cmplx(0.0,1.0);
 
#ifdef USE_ROLD
  using cpolybase<cmplx,ntype,dcmplx,N>::rold;
#endif
  template <class vtype>
    using pvecNm1 = typename std::conditional<(N>0),  pvector<vtype, N-1>,
          pvector<vtype, -1>>::type;
  template <class vtype>
    using pvecNp1 = typename std::conditional<(N>0),  pvector<vtype, N+1>,
          pvector<vtype, -1>>::type;
  using cpolyNm1 = typename std::conditional<(N>0),  cpoly<cmplx, N-1,ntype>,
        cpoly<cmplx, -1,ntype>>::type;
  template <class vtype>
    using pdvecNp1 = typename std::conditional<(N>0),  pvector<vtype, N+1>,
          pvector<vtype, -1>>::type;
  const int maxiter_polish=8;
  int imaxarg1,imaxarg2;
  ntype eps05, meps, maxf, maxf2, maxf3, minf, scalfact, cubic_rescal_fact;
  int maxdigits;
  ntype goaleps;
  ntype Kconv;
  bool gpolish;
  bool use_dbl_iniguess;
  bool guess_provided;

#ifdef USE_CONVEX_HULL
  using point=pair<ntype,ntype>; 

  typedef ntype coord_t;   // coordinate type
  typedef ntype coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2

  struct Point {
    coord_t x;
    coord_t y;
    bool operator <(const Point &p) const {
      return x < p.x || (x == p.x && y < p.y);
    }
  };
  // 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
  // Returns a positive value, if OAB makes a clockwise turn,
  // negative for counter-clockwise turn, and zero if the points are collinear.
  vector<Point> pts, convh;

  coord2_t cross(const Point &O, const Point &A, const Point &B)
    {
      return ((A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x));
    }

  // Returns a list of points on the convex hull in counter-clockwise order.
  // Note: the last point in the returned list is the same as the first one.
  // Method: Andrew’s monotone chain algorithm O(n) in this case since we do not need sorting
  vector<Point> upper_convex_hull(vector<Point> P)
    {
      int n = P.size(), k = 0;
      if (n <= 3) 
        return P;
      vector<Point> H(n+1);
      // points are already ordered in the present case!
      //sort(P.begin(), P.end());
      // we need just upper hull for Bini's algorithm
      // Build upper hull
      for (int i = n-1; i >= 0; i--) {
        while (k >= 2 && cross(H[k-1], H[k], P[i]) < 0) 
          {
            k--;
          }
        k++;
        H[k] = P[i];
      }
      H.resize(k+1);
      return H;
    }
#else
  vector<ntype> vk, uk;
#endif
  vector<cmplx> rg;
  vector<int> k;

public:
  void use_this_guess(pvector<dcmplx,N>& rg)
    {
      droots=rg;
      guess_provided=true;
    }
  void no_guess(void)
    {
      guess_provided=false;
    }
  void show(void)
    {
      show(NULL);
    }
  void show(const char* str)
    {
      int i;
      bool re, im;
      if (str!=NULL)
        cout <<  str;
      if (maxdigits <=0)
        maxdigits=30;
      for (i=n; i >= 0; i--)
        {
          re=false;
          im=false;
          if (real(coeff[i]) > 0.0)
            {
              if (i < n)
                cout << "+";
              if (imag(coeff[i])!=0)
                cout << "(";
              if (abs(real(coeff[i]))!=1.0 || imag(coeff[i])!=0)
                {
                  cout << setprecision(maxdigits) << abs(real(coeff[i]));
                }
              re=true;
            }
          else if (real(coeff[i]) < 0)
            { 
              cout << "-";
              cout << setprecision(maxdigits) << abs(real(coeff[i]));
              re=true;
            }
          else
            re=false;

          if (imag(coeff[i]) > 0.0)
            {
              if (i < n && re==true)
                cout << "+" << abs(imag(coeff[i]));
              im=true;
            }
          else if (imag(coeff[i]) < 0)
            { 
              cout << "-" << abs(imag(coeff[i]));
              im=true;
            }
          else
            im=false;
          if (im==true && re==true)
            cout << ")";
          if (im==true || re==true)
            {
              if ((im==true || abs(real(coeff[i]))!=1.0) && i > 0) 
                cout << "*";
              if ( i > 1)
                {
                  cout << "x^" << i;
                }
              else if (i==1)
                cout << "x";
            }
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
  // quadratic equation
  void solve_quadratic(pvector<cmplx,N>&sol)
    {
      cmplx acx,bcx,zx1,zx2,cdiskr,zxmax,zxmin;
      acx = coeff[1]/coeff[2];
      bcx = coeff[0]/coeff[2];
      cdiskr=sqrt(acx*acx-ntype(4.0)*bcx);
      zx1 = -ntype(0.5)*(acx+cdiskr);
      zx2 = -ntype(0.5)*(acx-cdiskr);
      if (abs(zx1) > abs(zx2))
        zxmax = zx1;
      else
        zxmax = zx2;
      if (zxmax==cmplx(0.0,0.0))
        zxmin=0;
      else
        zxmin = bcx/zxmax;

      sol[0] = zxmin;
      sol[1] = zxmax;
    }
  
  bool nr_aberth_real_rev(cmplx &r0, pvector<cmplx,N> &roots, int iac)
    {
      int j;
#ifndef BINI_CONV_CRIT
      ntype err;
#endif
      cmplx p1pc;
      const ntype EPS=goaleps;
      // root polishing by NR
      ntype tt[2], xcR, xcI, aR, aI, p1p[2];
#ifdef BINI_CONV_CRIT
      ntype pa[2], s;
#else
      ntype absp;
#endif
      ntype abx, p10, p0, x[2], p[2], p1[2], dx[2], invden;
      //cout << "xc=" << xc << "\n";
      xcR=real(r0);
      xcI=imag(r0);
      invden = 1.0/(xcR*xcR+xcI*xcI);
      x[0] = xcR*invden;// 1/x
      x[1] = -xcI*invden;
      p[0]=real(cmon[0]);
      p[1]=imag(cmon[0]);
      p1[0]=0.0;
      p1[1]=0.0;
      for (j=1;j<=n;j++) {
       p10 = p1[0];
       p1[0] = x[0]*p1[0] - x[1]*p1[1] + p[0];
       p1[1] = x[0]*p1[1] + x[1]*p10 + p[1];
       p0 = p[0];
       p[0] = x[0]*p[0] - x[1]*p[1] + real(cmon[j]);
       p[1] = x[0]*p[1] + x[1]*p0 + imag(cmon[j]);
     }

#ifdef BINI_CONV_CRIT
     s=acmon[n];
     abx = sqrt(xcR*xcR+xcI*xcI);
     pa[0]=real(cmon[n]);
     pa[1]=imag(cmon[n]);
     for (j=n-1; j >=0; j--) 
       {
         s=abx*s+acmon[j];
         p0 = pa[0];
         pa[0] = xcR*pa[0] - xcI*pa[1] + real(cmon[j]);
         pa[1] = xcR*pa[1] + xcI*p0 + imag(cmon[j]);
       }
     if (abs(cmplx(pa[0],pa[1])) <= 2.0*EPS*(4.0*ntype(n)+1)*s) // stopping criterion of bini 
       {
         return true;
       }
#else
     err=alpha[0];//abs(cmon[0])*(Kconv*m+1);
     abx=1.0/sqrt(xcR*xcR+xcI*xcI);
     for (j=1;j<=n;j++) {
       err=abx*err+alpha[j];//abs(cmon[j])*(Kconv*j+1);
     }
     absp=abs(cmplx(p[0],p[1]));
     if (absp <= EPS*err || (EPS*err <= minf && absp < minf))
       {
         return true;
       }
#endif 

      p1pc=cmplx(p1[0],p1[1])/cmplx(p[0],p[1]);
      p1p[0]=real(p1pc);
      p1p[1]=imag(p1pc);
      tt[0] = ntype(n)-x[0]*p1p[0]+x[1]*p1p[1];
      tt[1] = -x[0]*p1p[1]-x[1]*p1p[0];
      p1p[0] = x[0]*tt[0] - x[1]*tt[1];
      p1p[1] = x[1]*tt[0] + x[0]*tt[1];

      for (j=0; j < n; j++)
        {
          if (j==iac)
            continue;
          aR = xcR-real(roots[j]);
          aI = xcI-imag(roots[j]);
          invden = 1.0/(aR*aR+aI*aI);
          p1p[0] -= aR*invden;
          p1p[1] -= -aI*invden;
        }
      invden=1.0/(p1p[0]*p1p[0]+p1p[1]*p1p[1]);
      dx[0]=p1p[0]*invden;
      dx[1]=-p1p[1]*invden;
      r0 = cmplx(xcR-dx[0],xcI-dx[1]);
      return false;
    }
  bool nr_aberth_real(cmplx &r0, pvector<cmplx,N> &roots, int iac)
    {
      int j;
#ifndef BINI_CONV_CRIT
      ntype err;
#endif
      cmplx p1pc;
      const ntype EPS=goaleps;
      // root polishing by NR
      ntype aR, aI, p1p[2];
      ntype abx, p10, p0, x[2], p[2], p1[2], dx[2], invden;
#ifdef BINI_CONV_CRIT
      ntype s;
#else
      ntype absp;
#endif
      x[0] = r0.real();
      x[1] = r0.imag();
     //its=iter;
      p[0]=real(cmon[n]);
      p[1]=imag(cmon[n]);
#ifndef BINI_CONV_CRIT
      err=alpha[n];
#endif
      p1[0]=0.0;
      p1[1]=0.0;
      abx=sqrt(x[0]*x[0]+x[1]*x[1]);
      for (j=n-1;j>=0;j--) {
        p10 = p1[0];
        p1[0] = x[0]*p1[0] - x[1]*p1[1] + p[0];
        p1[1] = x[0]*p1[1] + x[1]*p10 + p[1];
        p0 = p[0];
        p[0] = x[0]*p[0] - x[1]*p[1] + real(cmon[j]);
        p[1] = x[0]*p[1] + x[1]*p0 + imag(cmon[j]);
#ifndef BINI_CONV_CRIT
        err=abx*err+alpha[j];
#endif
     }
#ifdef BINI_CONV_CRIT
     s=acmon[n];
     for (j=n-1; j >=0; j--) 
       {
         s=abx*s+acmon[j];
       }

     if (abs(cmplx(p[0],p[1])) <= 2.0*EPS*(4.0*ntype(n)+1)*s) // stopping criterion of bini 
       {
         return true;
       }
#else
     absp=abs(cmplx(p[0],p[1]));
     if (absp <= EPS*err || (EPS*err <= minf && absp < minf))
       {
         return true;
       }
#endif 

      p1pc=cmplx(p1[0],p1[1])/cmplx(p[0],p[1]);
      p1p[0]=real(p1pc);
      p1p[1]=imag(p1pc);
      for (j=0; j < n; j++)
        {
          if (j==iac)
            continue;
          aR = x[0]-real(roots[j]);
          aI = x[1]-imag(roots[j]);
          invden = 1.0/(aR*aR+aI*aI);
          p1p[0] -= aR*invden;
          p1p[1] -= -aI*invden;
        }
      invden=1.0/(p1p[0]*p1p[0]+p1p[1]*p1p[1]);
      dx[0]=p1p[0]*invden;
      dx[1]=-p1p[1]*invden;
      x[0] -= dx[0];
      x[1] -= dx[1];
      r0 = cmplx(x[0],x[1]);
      return false;
    }   
  bool nr_aberth_cmplx(cmplx &r0, pvector<cmplx,N> &roots, int iac)
    {
      int i;
#ifdef POLISH_NR_REAL
      cmplx povp1;
#endif
      const ntype EPS=goaleps;
      // root polishing by NR
      cmplx p, p1, p1p;
      ntype absp, abx;
#ifndef BINI_CONV_CRIT
      ntype err;
#else
      ntype K=2.0*EPS*(4.0*ntype(n)+1.0);
      ntype s;
#endif
      abx=abs(r0);
#ifndef BINI_CONV_CRIT
      err=alpha[n];
#endif
      p=cmon[n];
      p1=0;
#ifndef BINI_CONV_CRIT
      err=alpha[n];//abs(cmon[m])*(Kconv*m+1);
#endif
      for(i=n-1;i>=0;i--) {
        p1=p1*r0+p;
        p=p*r0+cmon[i];
#ifndef BINI_CONV_CRIT
        err=abx*err+alpha[i];
#endif
      }
#ifdef BINI_CONV_CRIT
      s=acmon[n];
      for (i=n-1; i >=0; i--) 
        {
          s=abx*s+acmon[i];
        }
      absp=abs(p);
      if (absp <= K*s) // stopping criterion of bini 
        {
          return true;
        }
#else
      absp=abs(p);
      if (absp <= EPS*err || (EPS*err <= minf && absp < minf))
        {
          return true;
        }
#endif
      p1p=p1/p;
      for (i=0; i < n; i++)
        {
          if (i==iac)
            continue;
          p1p -= ntype(1.0)/(r0-roots[i]);
        }
      r0 -= ntype(1.0)/p1p;
      return false;
    }
  bool nr_aberth_cmplx_rev(cmplx &r0, pvector<cmplx,N> &roots, int iac)
    {
      int i;
      const ntype EPS=goaleps;
      // root polishing by NR
      cmplx pa, p, p1, p1p, x;
      ntype absp, abx;
#ifndef BINI_CONV_CRIT
      ntype err;
#else
      ntype K=2.0*EPS*(4.0*ntype(n)+1.0);
      ntype s;
#endif
      x=ntype(1.0)/r0;
      p=cmon[0];
      p1=0;
      for(i=1;i<=n;i++) {
        p1=p1*x+p;
        p=p*x+cmon[i];
      }
#ifdef BINI_CONV_CRIT
      s=acmon[n];
      pa=cmon[n];
      abx=abs(r0);
      for (i=n-1; i >=0; i--) 
        {
          s=abx*s+acmon[i];
          pa=r0*pa+cmon[i];
        }
      absp=abs(pa);
      if (absp <= K*s) // stopping criterion of bini 
        {
          return true;
        }
#else
      err=alpha[0];//abs(cmon[0])*(Kconv*m+1);
      abx=1.0/abs(r0);
      for (i=1;i<=n;i++) {
        err=abx*err+alpha[i];//abs(cmon[j])*(Kconv*j+1);
      }
      absp=abs(p);
      if (absp <= EPS*err|| (EPS*err<=minf && absp<=minf)) 
        {
          return true;
        }
#endif

      p1p=p1/p;
      p1p=x*(ntype(n)-x*p1p);
      for (i=0; i < n; i++)
        {
          if (i==iac)
            continue;
          p1p -= ntype(1.0)/(r0-roots[i]);
        }
      r0 -= ntype(1.0)/p1p;
      return false;
    }

  void refine_root_maehly(cmplx &r0, pvector<cmplx,N> &roots, int iac)
    {
      int iter,i;
      ntype err;
#ifdef POLISH_NR_REAL
      cmplx povp1;
#endif
      // root polishing by NR
      cmplx r0old, p, sum;
      cmplx p1, p1c;
      ntype errold;
      for (iter=0; ; iter++)
        {
          p=cmon[n]*r0+cmon[n-1];
          p1=cmon[n];
          for(i=n-2;i>=0;i--) {
            p1=p+p1*r0;
            p=cmon[i]+p*r0;
          }
          if (iter > 0)
            errold=err;

          err = abs(p);//abs(p.real())+abs(p.imag());
          if (err==0)
            break;

          if (iter > 0 && err >= errold)
            {
              r0=r0old;
              break;
            }

          if (p1==cmplx(0,0))
            {
              break;
            }
          if (iter==maxiter_polish)
            break;
          r0old=r0;
          sum=0;
          for (i=0; i < iac; i++)
            sum += ntype(1.0)/(r0-roots[i]);
          sum *= p;
          p1c = p1 - sum;
          r0 -= p/p1c;
          if (isnan(abs(r0)) || isinf(abs(r0)))
            {
              r0=r0old;
              break;
            }
        }
    }
  void refine_root(cmplx &r0)
    {
      int iter,i;
      ntype err;
      // root polishing by NR
#ifdef POLISH_NR_REAL
      // root polishing by NR
      ntype r0old[2], p[2], r0new[2];
      ntype p1[2], p0, p10;
      ntype errold, invnp1;
      cmplx povp1;
      r0new[0]=r0.real();
      r0new[1]=r0.imag();
      for (iter=0; ;iter++)
        {
          p[0]=real(cmon[n])*r0new[0]-imag(cmon[n])*r0new[1] + real(cmon[n-1]);
          p[1]=real(cmon[n])*r0new[1]+imag(cmon[n])*r0new[0] + imag(cmon[n-1]);    
          p1[0]=real(cmon[n]);
          p1[1]=imag(cmon[n]);
          for(i=n-2;i>=0;i--) {
            p10=p1[0];
            p1[0]=p[0]+p1[0]*r0new[0]-p1[1]*r0new[1];
            p1[1]=p[1]+p10*r0new[1]+p1[1]*r0new[0];
            p0=p[0];
            p[0]=real(cmon[i])+p[0]*r0new[0]-p[1]*r0new[1];
            p[1]=imag(cmon[i])+p0*r0new[1]+p[1]*r0new[0];
          }
          if (iter > 0)
            errold=err;
          err = abs(cmplx(p[0],p[1]));// abs(p[0])+abs(p[1]);

          if (err==0)
            {
              break;
            }
          if (iter > 0 && err >= errold)
            {
              r0new[0]=r0old[0];
              r0new[1]=r0old[1];
              break;
            }
          if (p1[0]==0 && p1[1]==0)
            {
              break;
            }
          if (iter == maxiter_polish)
            break;
          r0old[0]=r0new[0];
          r0old[1]=r0new[1];
          povp1 = cmplx(p[0],p[1])/cmplx(p1[0],p1[1]);
          r0new[0] -= povp1.real();
          r0new[1] -= povp1.imag();
          if (isnan(r0new[0]) || isnan(r0new[1]) ||isinf(r0new[0]) || isinf(r0new[1]))
            {
              r0new[0]=r0old[0];
              r0new[1]=r0old[1];
              break; 
            }
        }
      r0=cmplx(r0new[0], r0new[1]);

#else
      cmplx r0old, p;
      cmplx p1;
      ntype errold;
      for (iter=0; ; iter++)
        {
          p=cmon[n]*r0+cmon[n-1];
          p1=cmon[n];
          for(i=n-2;i>=0;i--) {
            p1=p+p1*r0;
            p=cmon[i]+p*r0;
          }
          if (iter > 0)
            errold=err;

          err = abs(p);//abs(p.real())+abs(p.imag());
          if (err==0)
            break;

          if (iter > 0 && err >= errold)
            {
              r0=r0old;
              break;
            }

          if (p1==cmplx(0,0))
            {
              break;
            }
          if (iter==maxiter_polish)
            break;

          r0old=r0;
          r0 -= p/p1;

          if (isnan(abs(r0)) || isinf(abs(r0)))
            {
              r0=r0old;
              break;
            }
        }
#endif
    }

  ntype calc_upper_bound_kal(void)
    {
      // Kalantari's formula as found in McNamee Pan Vol. 1
      int k;
      static const ntype K= 1.0/0.682338;
      ntype ximax, xi;
      cmplx cnsq, term;

      cnsq=cmon[n-1]*cmon[n-1];
      for (k=4; k <= n+3; k++)
        {
          term=0;
          if (n-k+3 >= 0)
            term += cnsq*cmon[n-k+3]-cmon[n-2]*cmon[n-k+3];
          if (n-k+2 >= 0)
            term +=-cmon[n-1]*cmon[n-k+2];
          if (n-k+1 >= 0)
            term += cmon[n-k+1];
          xi = pow(abs(term),1.0/(k-1));
          if (k==4 || xi > ximax)
            ximax=xi;
        }
      ximax=K*ximax;
      if (isinf(ximax)||isnan(ximax))
        {
          //cout << "I cannot calculate the upper bound...\n";
          //cmon.show("coeff");     
          //exit(-1);
          return pow(maxf/1.618034,1.0/n);
        }

      return ximax;
    }
  void initial_guess(pvector<cmplx,N>& roots)
    {
      ntype sigma;
#ifndef USE_CONVEX_HULL
      ntype ukc, vkc, maxt=0.0, t;
#endif
      int q, i, j;
#ifdef USE_CONVEX_HULL
      Point p;
      pts.resize(n+1);
      for (i = 0; i <= n; i++)
        {
          p.x=ntype(i);
          if (abs(cmon[i])!=0)
            {
              p.y=log(abs(cmon[i]));
            }
          else
            {
              p.y=-maxf/1.618034;
            }
          pts[i] = p;
        }
      convh=upper_convex_hull(pts);
      q=convh.size()-1;
      if (q<=2)
        {
          q=2;
          k.resize(q+1);
          k[1]=n;
          k[2]=0;
        }
      else
        {
          k.resize(q+1);
        }
      i=1; 
      for (i=1; i <= q; i++)
        {
          k[i] = ((int)convh[i].x);
        }
#endif
#ifndef USE_CONVEX_HULL
      int kk;
      vk.resize(n+1);
      uk.resize(n+1);
      for (i=1; i <= n; i++) 
        {
          t=abs(cmon[i]/cmon[0]);
          if (i==1 || t > maxt)
            maxt=t;
        }
      uk[0] = 1.0/(1.0+maxt);   

      for (i=0; i < n; i++) 
        {
          t=abs(cmon[i]/cmon[n]);
          if (i==0 || t > maxt)
            maxt=t;
        }
      vk[n] = 1.0 + maxt;

      for (kk=1; kk <= n; kk++)
        {
          for (i=0; i < kk; i++) 
            {
              ukc=pow(abs(cmon[i]/cmon[kk]),1.0/(kk-i));
              if (i==0 || ukc > uk[kk])
                uk[kk]=ukc;
            }
        }
      for (kk=0; kk < n; kk++)
        {
          bool first=true;
          for (i=kk+1; i <= n; i++) 
            {
              if (cmon[i]==cmplx(0.0,0.0)) 
                continue;
              vkc=pow(abs(cmon[i]/cmon[kk]),1.0/(kk-i));
              if (first || vkc < vk[kk])
                {
                  first=false;
                  vk[kk]=vkc;
                }
            }
        } 
      k.resize(n+2); 
      q=1;
      for (kk=0; kk <= n; kk++)
        {
          if (uk[kk] <= vk[kk])
            {
              k[q] = kk; 
              q++;
            }
        }  
      q--;
      //sigma=2.0*M_PI*(drand48()-0.5);
      sigma=0.7;
      rg.resize(n); 
      int cc=0, dk;
      for (i=1; i < q; i++)
        {
          dk = k[i+1]-k[i];
          for (j=1; j <= dk; j++)
            {
              rg[(k[i]+j-1)] = polar(ntype(uk[k[i+1]]),ntype(j*2.0*pigr/dk+2.0*pigr*ntype(i)/ntype(n)+sigma));
              cc++;
            }
        }
#else
      //sigma=2.0*M_PI*(drand48()-0.5);
      ntype ukip1;
      sigma=0.7;
      rg.resize(n); 
      int cc=0;

      for (i=q-1; i >= 1; i--)
        {
          ukip1 = pow(abs(cmon[k[i+1]]/cmon[k[i]]),1.0/(k[i]-k[i+1]));
          for (j=1; j <= k[i]-k[i+1]; j++)
            {
              rg[cc] = polar(ukip1,ntype(j*2.0*pigr/(k[i]-k[i+1])+2.0*pigr*ntype(i)/ntype(n)+sigma));
              cc++;
            }
        }
#endif
      //std::sort(rg.begin(),rg.end(), [&] (cmplx a, cmplx b)-> bool {return abs(a) > abs(b);});
      for (i=0; i < n; i++)
        {
          roots[i]=rg[i]; //dal più piccolo al più grande
        }
    }

  void aberth(pvector<cmplx,N>& roots, bool polish=false)
    {
      bool ret;
#ifdef _OPENMP
      volatile bool fine=false;
#else
      bool fine=false;
#endif
      int i, nf=0, iter, is=0;
      int itmax=1000, nold;
      ntype r;
      cmplx cn;
      if (coeff[n]==cmplx(0.0,0.0))
        {
          cout << "That's not an " << n << " degree polynomial!\n";
          return;
        }
      i=0;
      cn = coeff[n];
      // if first is coefficients are zero we x=0 root with multiplicity is
      nold=n;
      while (coeff[i]==cmplx(0.0,0.0))
        {
          roots[n-1]=0;
          n--;
          i++;
        }
      is=i;
      cmon[n]=1.0;
      for (i=n-1; i >=0; i--)
        {
          cmon[i]=coeff[i+is]/cn;
        }
   
      for (i=0; i <=n; i++)
        {
          acmon[i] = abs(cmon[i]);
          alpha[i] = acmon[i]*Kconv*(ntype(i)+ntype(1.0));
        }
      //absolute value of coefficients are used in Bini's stopping criterion
      //we use coeff vector to store them since it won't be used anymore from here on
      if (guess_provided)
        {
          for (i=0; i < n; i++)
            {
              roots[i]=cmplx(droots[i]);
            }
        }
      else if (use_dbl_iniguess)
        {
          if constexpr (!is_same<cmplx,dcmplx>::value 
                        && !is_same<cmplx,complex<float>>::value
                        && !is_same<cmplx,complex<double>>::value)
            {
              pvecNp1<dcmplx> dcoeff;
              if constexpr (N < 0)
                dcoeff.allocate(n+1);
              for (i=0; i < n+1; i++)
                dcoeff[i]=dcmplx(cmon[i]);
              cpoly<dcmplx,N,dntype,dcmplx,dntype> drs;
              if constexpr (N < 0)
                drs.allocate(n);
              drs.iniguess_slow();
              drs.set_coeff(dcoeff);
              drs.find_roots(droots);
              for (i=0; i < n; i++)
                roots[i]=cmplx(droots[i]); 
            }
          else
            initial_guess(roots);
        }
      else
        {
          initial_guess(roots);
        }
      for (i=0;i < n; i++)
        found[i]=false;
      for (iter=0; iter < itmax && !fine; iter++)
        {
#ifdef USE_ROLD
          rold=roots;
#endif

#if defined(_OPENMP)
#pragma omp parallel for private(r) private(ret) //schedule(static)
#endif
          for (i=0; i < n; i++)
            {
              // refine conjugate pairs only once!
              if (found[i]) 
                continue;
              r=abs(roots[i]);
#ifdef USE_ROLD
#ifdef USE_ABERTH_REAL
              if (r > 1)
                ret = nr_aberth_real_rev(roots[i], rold, i);
              else
                ret = nr_aberth_real(roots[i], rold, i);
#else
              if (r > 1)
                ret = nr_aberth_cmplx_rev(roots[i], rold, i);
              else
                ret = nr_aberth_cmplx(roots[i], rold, i);
#endif
#else
#ifdef USE_ABERTH_REAL
              if (r > 1)
                ret = nr_aberth_real_rev(roots[i], roots, i);
              else
                ret = nr_aberth_real(roots[i], roots, i);
#else
              if (r > 1)
                ret= nr_aberth_cmplx_rev(roots[i], roots, i);
              else
                ret = nr_aberth_cmplx(roots[i], roots, i);
#endif
#endif
              if (ret)
                {
                  found[i]=true;
                  // all processess should wait here before updating nf
#if defined(_OPENMP)
#pragma omp critical // "critical" prevents multiple access of code between { ... } 
                    {
                      nf++;
                      if (nf==n)
                        {
                          fine=true;
                          i=n;// this is like a break but thread-safe
                        }
                    }
#else
                  nf++;
                  if (nf==n)
                    {
                      fine=true;
                      i=n;// this is like a break but thread-safe
                    }
#endif
                }
            }
        }
      if (nf < n)
        {
          cout << "Found " << nf << " roots out of " << n << "\n";
        }

      // refine 
      if (polish==true)
        {
          for (i=0; i < n; i++)
            refine_root_maehly(roots[i], roots, i);
        }
      n=nold;
    }
   

  pvecNp1<ntype> get_coeff()
    {
      return coeff;
    }
  int degree()
    {
      return n; 
    }
  
  void set_polish(bool p)
    {
      gpolish=p;
    }
  void solve_cubic(pvector<cmplx,N>& sol)
    {
      /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
       * (see sec. 2.2 in the manuscript) */ 
      cmplx a,b,c;
      cmplx K, Q, R, Q3, R2, A, B, Ap, Am, ApB, AmB;
      ntype theta, Q3r, R2r, Ar, Br;
      const ntype sqrt3=sqrt(ntype(3.0))/ntype(2.0);
      a = coeff[2]/coeff[3];
      b = coeff[1]/coeff[3];
      c = coeff[0]/coeff[3]; 
      Q = (Sqr(a) - cmplx(3.0)*b)/cmplx(9.0);
      R = (cmplx(2.0)*Sqr(a)*a - cmplx(9.0)*a*b + cmplx(27.0)*c)/cmplx(54.0);
      if (Q.imag()==0 && R.imag()==0)
        {
          Q3r = Sqr(Q.real())*Q.real();
          R2r = Sqr(R.real());
          if (R2r < Q3r)
            {
              theta = acos(R.real()/sqrt(Q3r));
              sol[0] = -ntype(2.0)*sqrt(Q.real())*cos(theta/ntype(3.0))- a/cmplx(3.0);
              sol[1] = -ntype(2.0)*sqrt(Q.real())*cos((theta+ntype(2.0)*pigr)/ntype(3.0)) - a/cmplx(3.0);
              sol[2] = -ntype(2.0)*sqrt(Q.real())*cos((theta-ntype(2.0)*pigr)/ntype(3.0)) - a/cmplx(3.0);
            }
          else
            {
              Ar = -copysign((ntype)1.0,R.real())*cbrt(abs(R.real()) + sqrt(R2r - Q3r));
              if (Ar==0.0)
                Br=0.0;
              else
                Br = Q.real()/Ar;
              sol[0] = (Ar+Br) - a/cmplx(3.0);
              sol[1] = -ntype(0.5)*(Ar+Br)-a/cmplx(3.0)+I*sqrt3*(Ar-Br);
              sol[2] = -ntype(0.5)*(Ar+Br)-a/cmplx(3.0)-I*sqrt3*(Ar-Br);
            }
        }
      else
        {
          R2 = R*R;
          Q3 = Q*Q*Q;
          K=sqrt(R2 - Q3);
          Ap = -exp(log(R + K)/cmplx(3.0));
          Am = -exp(log(R - K)/cmplx(3.0));
          if (abs(Ap) > abs(Am))
            A=Ap;
          else 
            A=Am;
          if (A==cmplx(0.0))
            B=cmplx(0.0);
          else
            B = Q/A;
          ApB=A+B;
          AmB=A-B;
          sol[0] = -a/cmplx(3.0) + ApB; 
          sol[1] = -a/cmplx(3.0) - cmplx(0.5)*ApB + I*sqrt3*(AmB);
          sol[2] = -a/cmplx(3.0) - cmplx(0.5)*ApB - I*sqrt3*(AmB);
        }
    }
  //
  // find roots by default uses aberth method which is faster than laguerre implicit method
  inline void find_roots(pvector<cmplx,N>& roots)
    {
      if constexpr (N > 0)
        {  
          if constexpr (N==2)
            {
              solve_quadratic(roots); 
            }
          else if constexpr (N==3)
            {
              solve_cubic(roots);
            }
          else if constexpr (N==4)
            {
              quar.set_coeff(coeff);
              quar.find_roots(roots);
            }
          else
            {
              aberth(roots, gpolish);
            }
        }
      else
        {
          if (n==2)
            {
              solve_quadratic(roots); 
            }
          else if (n==3)
            {
              solve_cubic(roots);
            }
          else if (n==4)
            {
              quar.set_coeff(coeff);
              quar.find_roots(roots);
            }
          else
            {
              aberth(roots, gpolish);
            }
        }
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
  void set_output_prec(ntype e)
    {
      goaleps=e;
    } 
  ntype get_output_prec(void)
    {
      return goaleps;
    } 
  void iniguess_slow(void)
    {
      use_dbl_iniguess=false;
    }
  void iniguess_fast(void)
    {
      use_dbl_iniguess=true;
    }
  void set_show_digits(int p)
    {
      maxdigits=p;
    }

  void init_const(void)
    {
      meps = epsilon();
      eps05 = pow(numeric_limits<ntype>::epsilon(),0.5);
      maxf= getmax();
      minf = numeric_limits<ntype>::min();
      maxdigits =(numeric_limits<ntype>::digits10)-1;
      //cout << "numeric digits=" << maxdigits << "\n";

      maxf2 = pow(maxf,0.5)/10.0;
      maxf3 = pow(maxf,1.0/3.0)/10.0;
      scalfact = pow(maxf,1.0/4.0)/1.618034;
      cubic_rescal_fact = pow(maxf, 1.0/3.0)/1.618034;
      goaleps=numeric_limits<ntype>::epsilon(); 
      //Kconv=sqrt(ntype(2.0))*ntype(2.0)+ntype(1.0);  
      Kconv=3.8;
      gpolish=false;
      guess_provided=false;
      // calculate initial guess using long double precision (much faster!)
      // MPSolve calculate the initial guess using double.
      // It is not recommended to set this to false if using multiple precision
      // since you will experience a slowing down by a factor 3.
      use_dbl_iniguess=true; 
    }
  cpoly()
    {
      init_const();  
    }

  cpoly(int nc): cpolybase<cmplx,ntype,dcmplx,N>(nc)
  {
    init_const();
  }
};
#endif

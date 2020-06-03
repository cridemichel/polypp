#ifndef _RPOLY_
#define _RPOLY_
#include "./pvector.hpp"
#include "./cpoly.hpp"
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
template<class T>
inline const T MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

template <class ntype, int N, class cmplx> 
class rpoly_base_static 
{
public:
  int n;
  constexpr static int dynamic = false;
  pvector<ntype, N+1> coeff;
  pvector<ntype, N+1> cmon;
  cpoly<cmplx,N,ntype> cpol;
  quartic<ntype,cmplx,false> quar;

  void set_coeff(pvector<ntype,N+1> v)
    {
      coeff = v;
      cmon[n]=1.0;
      for (int i=n-1; i >=0; i--)
        {
          cmon[i]=coeff[i]/coeff[n];
        }

    }
  rpoly_base_static()
    {
      n=N;
    }
  rpoly_base_static(int nc): rpoly_base_static()
  {
    n=nc;
  }
};
template <class ntype, int N, class cmplx> 
class rpoly_base_dynamic 
{
public:
  int n;
  constexpr static int dynamic = true;
  pvector<ntype> coeff;
  pvector<ntype> cmon;
  cpoly<cmplx,-1,ntype> cpol;
  quartic<ntype,cmplx,true> quar;

  rpoly_base_dynamic() = default;
  void set_coeff(pvector<ntype,-1> v)
    {
      coeff = v;
      cmon[n]=1.0;
      for (int i=n-1; i >=0; i--)
        {
          cmon[i]=coeff[i]/coeff[n];
        }
    }

  void use_vec(int nc, ntype* coeffv, ntype* cmonv)
    {
      coeff.use_vec(nc+1,coeffv);
      cmon.use_vec(nc+1,cmonv);
      n=nc;
    }
 
  rpoly_base_dynamic(int nc): coeff(nc+1), cmon(nc+1), cpol(nc)
    {
      n=nc;
    }
  ~rpoly_base_dynamic() = default;
  void deallocate(void)
    {
      coeff.deallocate();
      cmon.deallocate();
      cpol.deallocate();
    }
  void allocate(int nc)
    {
      n=nc;
      coeff.allocate(n+1);
      cmon.allocate(n+1);
      cpol.allocate(n);
    }
};

template <class ntype, int N, class cmplx> using rpolybase = 
typename std::conditional<(N>0), rpoly_base_static <ntype, N, cmplx>,
	 rpoly_base_dynamic <ntype, N, cmplx>>::type;
 
template <class ntype, int N=-1, class cmplx=complex<ntype>> 
class rpoly: public numeric_limits<ntype>, public rpolybase<ntype,N, cmplx> {
  using rpolybase<ntype,N,cmplx>::n;
  using rpolybase<ntype,N,cmplx>::coeff;
  using rpolybase<ntype,N,cmplx>::cmon;
  using rpolybase<ntype,N,cmplx>::cpol;
  using rpolybase<ntype,N,cmplx>::quar;

  const ntype pigr=acos(ntype(-1.0));
  const cmplx I = cmplx(0.0,1.0);

  template <class vtype>
  using pvecNm1 = typename std::conditional<(N>0), pvector<vtype, N-1>,
	 pvector<vtype, -1>>::type;
  template <class vtype>
  using pvecNp1 = typename std::conditional<(N>0), pvector<vtype, N+1>,
	 pvector<vtype, -1>>::type;
  template <class vtype>
  using pvecNm2 = typename std::conditional<(N>0), pvector<vtype, N-2>,
	 pvector<vtype, -1>>::type;
  using rpolyNm1 = typename std::conditional<(N>0), rpoly<ntype, N-1,cmplx>,
	 rpoly<ntype, -1,cmplx>>::type;
  using rpolyNm2 = typename std::conditional<(N>0), rpoly<ntype, N-2,cmplx>,
	 rpoly<ntype, -1,cmplx>>::type;

 template <class vtype, int NT>
  using stlarr = typename std::conditional<(NT>0), std::array<vtype,NT-1>, std::vector<vtype>>::type;

  ntype px, dpx; 
  const int maxiter_polish=8;
  int imaxarg1,imaxarg2;
  ntype eps05, meps, maxf, maxf2, maxf3, scalfact, cubic_rescal_fact;
  int maxdigits;
  ntype goaleps;
  bool deflated;
  ntype oqs_max2(ntype a, ntype b)
    {
      if (a >= b)
	return a;
      else
	return b;
    }
  ntype oqs_max3(ntype a, ntype b, ntype c)
    {
      ntype t;
      t = oqs_max2(a,b);
      return oqs_max2(t,c);
    }
  //void oqs_quartic_solver(pvector<cmplx,N>& roots);
  inline void solve_quadratic(pvector<cmplx,N>& sol);
  inline void solve_cubic_analytic(pvector<cmplx,N>& sol);
  inline void oqs_solve_quadratic(ntype a, ntype b, cmplx roots[2]);
public:
  void show(void)
    {
      show(NULL);
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

  pvecNp1<ntype> get_coeff()
    {
      return coeff;
    }
  int degree()
    {
      return n; 
    }

  void aberth(pvector<cmplx,N>& roots)
    {
      /* use aberth method */
      cpol.set_coeff(coeff);
      cpol.find_roots(roots);
    }

  void find_roots(pvector<cmplx,N>& roots)
    {
      if constexpr (N < 0)
        {
          if (n==1)
            {
              cout << "What?!? You are not able to solve a linear equation, come on!";
              exit(-1);
            }
          else if (n==2)
            {
              solve_quadratic(roots);
            }
          else if (n==3) 
            {
              solve_cubic_analytic(roots);
            }
          else if (n==4)
            {
              quar.set_coeff(coeff);
              quar.find_roots(roots);
              //oqs_quartic_solver(roots);
            }
          else 
            {
              /* use aberth method */
              cpol.set_coeff(coeff);
              cpol.find_roots(roots);
            }
        }
      else
        {
          if constexpr (N==1)
            {
              cout << "What?!? You are not able to solve a linear equation, come on!";
              exit(-1);
            }
          else if constexpr (N==2)
            {
              solve_quadratic(roots);
            }
          else if constexpr (N==3) 
            {
              solve_cubic_analytic(roots);
            }
          else if constexpr (N==4)
            {
              quar.set_coeff(coeff);
              quar.find_roots(roots);
              //oqs_quartic_solver(roots);
            }
          else 
            {
              /* use aberth method */
              cpol.set_coeff(coeff);
              cpol.find_roots(roots);
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
   void init_const(void)
    {
      meps = epsilon();
      eps05 = pow(numeric_limits<ntype>::epsilon(),0.5);
      maxf= getmax();
      maxdigits = (numeric_limits<ntype>::digits10)-1;
      maxf2 = pow(maxf,0.5)/10.0;
      maxf3 = pow(maxf,1.0/3.0)/10.0;
      scalfact = pow(maxf,1.0/4.0)/1.618034;
      cubic_rescal_fact = pow(maxf, 1.0/3.0)/1.618034;
      goaleps=numeric_limits<ntype>::epsilon();   
      deflated=false;
   }
  void set_output_prec(ntype e)
    {
      goaleps=e;
    } 
  ntype get_output_prec(void)
    {
      return goaleps;
    } 
  rpoly(): rpolybase<ntype,N,cmplx>()
    {
      init_const();
    }

  rpoly(int nc): rpolybase<ntype,N,cmplx>(nc)
    {
      init_const(); 
    }
};
// quadratic equation
template<class ntype, int N, class cmplx> void rpoly<ntype,N,cmplx>::solve_quadratic(pvector<cmplx,N>&sol)
{
  cmplx r[2];
  ntype a,b;
  a = coeff[1]/coeff[2];
  b = coeff[0]/coeff[2];
  oqs_solve_quadratic(a, b, r);
  sol[0] = r[0];
  sol[1] = r[1];
}
// cubic polynomial
template <class ntype, int N, class cmplx> void rpoly<ntype,N, cmplx>::solve_cubic_analytic(pvector<cmplx,N>& sol)
{
  /* solve the cubic coeff[3]*x^3 + coeff[2]*x^2 +  coeff[1]*x + coeff[0] = 0
   * according to the method described in Numerical Recipe book */  
  ntype a, b, c, Q, R, theta, Q3, R2, A, B;
  const ntype sqrt32=sqrt((ntype)3.0)/2.0;
  a = coeff[2]/coeff[3];
  b = coeff[1]/coeff[3];
  c = coeff[0]/coeff[3];
  Q = (Sqr(a) - 3.0*b)/9.0;
  R = (2.0*Sqr(a)*a - 9.0*a*b + 27.0*c)/54.0;
  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
      theta = acos(R/sqrt(Q3));
      sol[0] = -2.0*sqrt(Q)*cos(theta/3.0)- a/3.0;
      sol[1] = -2.0*sqrt(Q)*cos((theta+2.0*M_PI)/3.0) - a/3.0;
      sol[2] = -2.0*sqrt(Q)*cos((theta-2.0*M_PI)/3.0) - a/3.0;
    }
  else
    {
      A = -copysign((ntype)1.0,R)*pow(abs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol[0] = (A+B) - a/3.0;
      sol[1] = cmplx(-0.5*(A+B)-a/3.0,sqrt32*(A-B));
      sol[2] = conj(sol[1]);
      //sol[1] = -0.5*(A+B)-a/3.0+cmplx(0,1)*sqrt32*(A-B);
      //sol[2] = -0.5*(A+B)-a/3.0-cmplx(0,1)*sqrt32*(A-B);
    }
}

template <class ntype, int N, class cmplx> void  rpoly<ntype,N, cmplx>::oqs_solve_quadratic(ntype a, ntype b, cmplx roots[2])
{ 
  ntype div,sqrtd,diskr,zmax,zmin;
  diskr=a*a-4*b;   
  if(diskr>=0.0)
    {
      if(a>=0.0)
	div=-a-sqrt(diskr);
      else
	div=-a+sqrt(diskr);

      zmax=div/2;

      if(zmax==0.0)
	zmin=0.0;
      else
	zmin=b/zmax;
      roots[0]=cmplx(zmax,0.0);
      roots[1]=cmplx(zmin,0.0);
    } 
  else
    {   
      sqrtd = sqrt(-diskr);
      roots[0]=cmplx(-a/2,sqrtd/2);
      roots[1]=cmplx(-a/2,-sqrtd/2);      
    }   
}
#endif

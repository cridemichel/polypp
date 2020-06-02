#ifndef _RPOLY_
#define _RPOLY_
#include "./pvector.hpp"
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
#define Sqr(x) ((x)*(x))
//#define FAST_MATH
using namespace std;
template<class T>
inline const T MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

template <class ntype, int N> 
class rpoly_base_static 
{
public:
  int n;
  constexpr static int dynamic = false;
  pvector<ntype, N+1> coeff;
  pvector<ntype, N+1> cmon;
  
  void set_coeff(pvector<ntype,N+1> v)
    {
      coeff = v;
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
template <class ntype, int N> 
class rpoly_base_dynamic 
{
public:
  int n;
  constexpr static int dynamic = true;
  pvector<ntype> coeff;
  pvector<ntype> cmon;
  
  rpoly_base_dynamic() = default;
  void set_coeff(pvector<ntype,-1> v)
    {
      coeff = v;
    }

  void use_vec(int nc, ntype* coeffv, ntype* cmonv)
    {
      coeff.use_vec(nc+1,coeffv);
      cmon.use_vec(nc+1,cmonv);
      n=nc;
    }
 
  rpoly_base_dynamic(int nc): coeff(nc+1), cmon(nc+1)
    {
      n=nc;
    }
  ~rpoly_base_dynamic() = default;
  void deallocate(void)
    {
      coeff.deallocate();
      cmon.deallocate();
    }
  void allocate(int nc)
    {
      n=nc;
      coeff.allocate(n+1);
      cmon.allocate(n+1);
    }
};

template <class ntype, int N> using rpolybase = 
typename std::conditional<(N>0), rpoly_base_static <ntype, N>,
	 rpoly_base_dynamic <ntype, N>>::type;
 
template <class ntype, int N=-1, class cmplx=complex<ntype>> 
class rpoly: public numeric_limits<ntype>, public rpolybase<ntype,N> {
  using rpolybase<ntype,N>::n;
  using rpolybase<ntype,N>::coeff;
  using rpolybase<ntype,N>::cmon;
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
  inline void oqs_quartic_solver(pvector<cmplx,N>& roots);
  inline void oqs_solve_cubic_analytic_depressed_handle_inf(ntype b, ntype c, ntype& sol);
  inline void oqs_solve_cubic_analytic_depressed(ntype b, ntype c, ntype& sol);
  inline void oqs_calc_phi0(ntype a, ntype b, ntype c, ntype d, ntype& phi0, int scaled);
  inline ntype oqs_calc_err_ldlt(ntype b, ntype c, ntype d, ntype d2, ntype l1, ntype l2, ntype l3);
  inline ntype oqs_calc_err_abcd_cmplx(ntype a, ntype b, ntype c, ntype d, 
	    			       cmplx aq, cmplx bq, cmplx cq, cmplx dq);
  inline ntype oqs_calc_err_abcd(ntype a, ntype b, ntype c, ntype d, ntype aq, ntype bq, ntype cq, ntype dq);
  inline ntype oqs_calc_err_abc(ntype a, ntype b, ntype c, ntype aq, ntype bq, ntype cq, ntype dq); 
  inline void oqs_NRabcd(ntype a, ntype b, ntype c, ntype d, ntype& AQ, ntype& BQ, ntype& CQ, ntype& DQ);
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
  inline void find_roots(pvector<cmplx,N>& roots)
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
              oqs_quartic_solver(roots);
            }
          else 
            {
              /* use aberth method */
              pvector<cmplx, -1> co(n+1);
              pvector<cmplx, -1> ro(n);
              cpoly<cmplx,-1,ntype> pol(n);
              for (int i=0; i <= n; i++)
                co[i] = cmplx(coeff[i],0.0);
              pol.set_coeff(co);
              pol.find_roots(ro);
              for (int i=0; i < n; i++)
                roots[i] = ro[i].real();

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
              oqs_quartic_solver(roots);
            }
          else 
            {
              /* use aberth method */
              pvector<cmplx, N+1> co;
              pvector<cmplx, N> ro;
              cpoly<cmplx,N,ntype> pol;
              for (int i=0; i <= N; i++)
                co[i] = cmplx(coeff[i],0.0);
              pol.set_coeff(co);
              pol.find_roots(ro);
              for (int i=0; i < N; i++)
                roots[i] = ro[i].real();
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
  rpoly(): rpolybase<ntype,N>()
    {
      init_const();
    }

  rpoly(int nc): rpolybase<ntype,N>(nc)
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

// quartics with OQS
template <class ntype, int N, class cmplx> void rpoly<ntype,N,cmplx>::oqs_solve_cubic_analytic_depressed_handle_inf(ntype b, ntype c, ntype& sol)
{
 /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
  * where coefficients b and c are large (see sec. 2.2 in the manuscript) */ 
  ntype Q, R, theta, A, B, QR, QRSQ, KK, sqrtQ, RQ;;
  const ntype PI2=M_PI/2.0, TWOPI=2.0*M_PI;
#ifdef FAST_MATH
  ntype rq3;
#endif
 
  Q = -b/3.0;
  R = 0.5*c;
  if (R==0)
    {
      if (b <= 0)
	{
	  sol=sqrt(-b);
	}
      else
	{
	  sol=0;
	}
      return;
    }
  
  if (abs(Q) < abs(R))
    {
      QR=Q/R;
      QRSQ=QR*QR; 
      KK=1.0 - Q*QRSQ;
    }
  else
    {
      RQ = R/Q;
      KK = copysign((ntype)1.0,Q)*(RQ*RQ/Q-1.0);
    }

  if (KK < 0.0)
    {
      sqrtQ=sqrt(Q);
#ifdef FAST_MATH
      // se si use -Ofast rq3 può essere >1 o < -1 di quantità minori di machine epsilon 
      // causando dei NaN
      rq3 = (R/abs(Q))/sqrtQ;
      if (rq3 > 1.0)
        theta = 1.0;
      else if (rq3 < -1.0)
        theta  = M_PI;
      else
        theta = acos(rq3);
#else
      theta = acos((R/abs(Q))/sqrtQ);
#endif
      if (theta < PI2) 
	sol = -2.0*sqrtQ*cos(theta/3.0);
      else 
	sol = -2.0*sqrtQ*cos((theta+TWOPI)/3.0);
    }
  else
    {
      if (abs(Q) < abs(R))
	A = -copysign((ntype)1.0,R)*cbrt(abs(R)*(1.0+sqrt(KK)));
      else
	{
	  A = -copysign((ntype)1.0,R)*cbrt(abs(R)+sqrt(abs(Q))*abs(Q)*sqrt(KK));
	}
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol = A+B;
    }
}
template <class ntype, int N, class cmplx> void rpoly<ntype,N,cmplx>::oqs_solve_cubic_analytic_depressed(ntype b, ntype c, ntype& sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * (see sec. 2.2 in the manuscript) */ 
  ntype Q, R, theta, Q3, R2, A, B, sqrtQ;
#ifdef FAST_MATH
  ntype rq3;
#endif
  Q = -b/3.0;
  R = 0.5*c;
  // these number could be made larger for long double */
  //if (abs(Q) > 1E102 || abs(R) > 1E154)
  if (abs(Q) > maxf3 || abs(R) > maxf2)
    {
      oqs_solve_cubic_analytic_depressed_handle_inf(b, c, sol);
      return;
    }
  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
#ifdef FAST_MATH
      // se si use -Ofast rq3 può essere >1 o < -1 di quantità minori di machine epsilon 
      // causando dei NaN
      rq3 = R/sqrt(Q3);
      if (rq3 > 1.0)
        theta = 1.0;
      else if (rq3 < -1.0)
        theta  = M_PI;
      else
        theta = acos(rq3);
#else
      theta = acos(R/sqrt(Q3));
#endif
      sqrtQ=-2.0*sqrt(Q);
      if (theta < M_PI/2) 
	sol = sqrtQ*cos(theta/3.0);
      else 
	sol = sqrtQ*cos((theta+2.0*M_PI)/3.0);
    }
  else
    {
      A = -copysign((ntype)1.0,R)*cbrt(abs(R) + sqrt(R2 - Q3));
       if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol = A+B; /* this is always largest root even if A=B */
    }
}
template <class ntype, int N, class cmplx> void  rpoly<ntype,N, cmplx>::oqs_calc_phi0(ntype a, ntype b, ntype c, ntype d, ntype& phi0, int scaled)
{
  /* find phi0 as the dominant root of the depressed and shifted cubic 
   * in eq. (64) (see also the discussion in sec. 2.2 of the manuscript) */
  ntype rmax, g,h,gg,hh,aq,bq,cq,dq,s,diskr;
  ntype maxtt, xxx, gx, x, xold, f, fold, df, xsq;
  ntype ggss, hhss, dqss, aqs, bqs, cqs, rfact, rfactsq; 
  int iter;

  diskr=9*a*a-24*b;                    
  /* eq. (67) */
  if(diskr > 0.0)
    { 
      diskr=sqrt(diskr);
      if(a > 0.0)
	s=-2*b/(3*a+diskr);                     
      else
	s=-2*b/(3*a-diskr);                      
    }
  else
    {      
      s=-a/4;                                    
    }
  /* eqs. (63) */
  aq=a+4*s;                                      
  bq=b+3*s*(a+2*s);                              
  cq=c+s*(2*b+s*(3*a+4*s));                      
  dq=d+s*(c+s*(b+s*(a+s)));                      
  gg=bq*bq/9;
  hh=aq*cq;     
  
  g=hh-4*dq-3*gg;                       /* eq. (60) */  
  h=(8*dq+hh-2*gg)*bq/3-cq*cq-dq*aq*aq; /* eq. (61) */          
  oqs_solve_cubic_analytic_depressed(g, h, rmax);
  if (isnan(rmax) || isinf(rmax))
    {
      oqs_solve_cubic_analytic_depressed_handle_inf(g, h, rmax);
      if ((isnan(rmax) || isinf(rmax)) && scaled)
	{
	  // try harder: rescale also the depressed cubic if quartic has been already rescaled
	  rfact = cubic_rescal_fact; 
	  rfactsq = rfact*rfact;
	  ggss = gg/rfactsq;
	  hhss = hh/rfactsq;
	  dqss = dq/rfactsq;
	  aqs = aq/rfact;
	  bqs = bq/rfact;
	  cqs = cq/rfact;
	  ggss=bqs*bqs/9.0;
	  hhss=aqs*cqs;   
	  g=hhss-4.0*dqss-3.0*ggss;                       
	  h=(8.0*dqss+hhss-2.0*ggss)*bqs/3-cqs*(cqs/rfact)-(dq/rfact)*aqs*aqs; 
	  oqs_solve_cubic_analytic_depressed(g, h, rmax);
	  if (isnan(rmax) || isinf(rmax))
	    {
	      oqs_solve_cubic_analytic_depressed_handle_inf(g, h, rmax);
	    }
	  rmax *= rfact;
	}
    }
  /* Newton-Raphson used to refine phi0 (see end of sec. 2.2 in the manuscript) */
  x = rmax;
  xsq=x*x;
  xxx=x*xsq;
  gx=g*x;
  f = x*(xsq + g) + h;
  if (abs(xxx) > abs(gx))
    maxtt = abs(xxx);
  else
    maxtt = abs(gx);
  if (abs(h) > maxtt)
    maxtt = abs(h);

  if (abs(f) > maxtt)
    {
      for (iter=0; iter < 8; iter++)
	{   
	  df =  3.0*xsq + g;
	  if (df==0)
	    {
	      break;
	    }
	  xold = x;
	  x += -f/df;
	  fold = f;
	  xsq = x*x;
	  f = x*(xsq + g) + h;
	  if (f==0)
	    {
	      break;
	    } 

	  if (abs(f) >= abs(fold))
	    {
	      x = xold;
	      break;
	    }
    	}
    }
  phi0 = x;
}
template <class ntype, int N, class cmplx> ntype  rpoly<ntype,N,cmplx>::oqs_calc_err_ldlt(ntype b, ntype c, ntype d, ntype d2, ntype l1, ntype l2, ntype l3)
{
  /* Eq. (21) in the manuscript */
  ntype sum;
  if (b==0.0)
    sum =  abs(d2 + l1*l1 + 2.0*l3);
  else 
    sum =  abs(((d2 + l1*l1 + 2.0*l3)-b)/b);
  if (c==0.0)
    sum += abs(2.0*d2*l2 + 2.0*l1*l3);
  else 
    sum += abs(((2.0*d2*l2 + 2.0*l1*l3)-c)/c);
  if (d==0.0)
    sum += abs(d2*l2*l2 + l3*l3);
  else 
    sum += abs(((d2*l2*l2 + l3*l3)-d)/d);
  return sum;
}
template <class ntype, int N, class cmplx> 
ntype rpoly<ntype,N, cmplx>::oqs_calc_err_abcd_cmplx(ntype a, ntype b,  ntype c, ntype d, cmplx aq, 
                                                         cmplx bq, cmplx cq, cmplx dq)
{
  /* Eq. (53) in the manuscript for complex alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq) */
  ntype sum;
  sum = (d==0)?abs(bq*dq):abs((bq*dq-d)/d);
  sum += (c==0)?abs(bq*cq + aq*dq):abs(((bq*cq + aq*dq) - c)/c);
  sum +=(b==0)?abs(bq + aq*cq + dq):abs(((bq + aq*cq + dq) - b)/b);
  sum +=(a==0)?abs(aq + cq):abs(((aq + cq) - a)/a);
  return sum;
}
template <class ntype, int N, class cmplx> ntype rpoly<ntype,N,cmplx>::oqs_calc_err_abcd(ntype a, ntype b, ntype c, ntype d, ntype aq, ntype bq, ntype cq, ntype dq)
{
  /* Eq. (53) in the manuscript for real alpha1 (aq), beta1 (bq), alpha2 (cq) and beta2 (dq)*/
  ntype sum;

  if (d==0.0)
    sum = abs(bq*dq);
  else
    sum = abs((bq*dq-d)/d);

  if (c==0.0)
    sum += abs(bq*cq + aq*dq);
  else 
    sum += abs(((bq*cq + aq*dq) - c)/c);

  if (b==0.0)
    sum +=abs(bq + aq*cq + dq);
  else 
    sum +=abs(((bq + aq*cq + dq) - b)/b);
  if (a==0.0)
    sum +=abs(aq + cq);
  else 
    sum +=abs(((aq + cq) - a)/a);
  return sum;
}
template <class ntype, int N, class cmplx> ntype  rpoly<ntype,N,cmplx>::oqs_calc_err_abc(ntype a, ntype b, ntype c, ntype aq, ntype bq, ntype cq, ntype dq)
{
  /* Eq. (40) in the manuscript */
  ntype sum;
  if (c==0.0)
    sum = abs(bq*cq + aq*dq);
  else
    sum = abs(((bq*cq + aq*dq) - c)/c);
  if (b==0.0)
  sum +=abs(bq + aq*cq + dq);
  else 
  sum +=abs(((bq + aq*cq + dq) - b)/b);
  if (a==0.0)
    sum +=abs(aq + cq);
  else 
    sum +=abs(((aq + cq) - a)/a);
  return sum;
}
template <class ntype, int N, class cmplx> void rpoly<ntype,N, cmplx>::oqs_NRabcd(ntype a, ntype b, ntype c, ntype d, ntype& AQ, ntype& BQ, ntype& CQ, ntype& DQ)
{
  /* Newton-Raphson described in sec. 2.3 of the manuscript for complex
   * coefficients a,b,c,d */
  int iter, k1, k2;
  ntype x02, errf, errfold, xold[4], x[4], dx[4], det, Jinv[4][4], fvec[4], vr[4];
  x[0] = AQ;
  x[1] = BQ;
  x[2] = CQ;
  x[3] = DQ;
  vr[0] = d;
  vr[1] = c;
  vr[2] = b;
  vr[3] = a;
  fvec[0] = x[1]*x[3] - d;
  fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
  fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
  fvec[3] = x[0] + x[2] - a; 
  errf=0;
  for (k1=0; k1 < 4; k1++)
    {
      if (vr[k1]==0)
        errf += abs(fvec[k1]);
      else
        errf +=abs(fvec[k1]/vr[k1]);
    }
  for (iter = 0; iter < 8; iter++)
    {
      x02 = x[0]-x[2];
      det = x[1]*x[1] + x[1]*(-x[2]*x02 - 2.0*x[3]) + x[3]*(x[0]*x02 + x[3]);
      if (det==0.0)
	break;
      Jinv[0][0] = x02;
      Jinv[0][1] = x[3] - x[1];
      Jinv[0][2] = x[1]*x[2] - x[0]*x[3];
      Jinv[0][3] = -x[1]*Jinv[0][1] - x[0]*Jinv[0][2]; 
      Jinv[1][0] = x[0]*Jinv[0][0] + Jinv[0][1];
      Jinv[1][1] = -x[1]*Jinv[0][0];
      Jinv[1][2] = -x[1]*Jinv[0][1];   
      Jinv[1][3] = -x[1]*Jinv[0][2];
      Jinv[2][0] = -Jinv[0][0];
      Jinv[2][1] = -Jinv[0][1];
      Jinv[2][2] = -Jinv[0][2];
      Jinv[2][3] = Jinv[0][2]*x[2] + Jinv[0][1]*x[3];
      Jinv[3][0] = -x[2]*Jinv[0][0] - Jinv[0][1];
      Jinv[3][1] = Jinv[0][0]*x[3];
      Jinv[3][2] = x[3]*Jinv[0][1];
      Jinv[3][3] = x[3]*Jinv[0][2];
      for (k1=0; k1 < 4; k1++)
	{
	  dx[k1] = 0;
	  for (k2=0; k2 < 4; k2++)
	    dx[k1] += Jinv[k1][k2]*fvec[k2];
	}
      for (k1=0; k1 < 4; k1++)
      	xold[k1] = x[k1];

      for (k1=0; k1 < 4; k1++)
	{
	  x[k1] += -dx[k1]/det;
	}
      fvec[0] = x[1]*x[3] - d;
      fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
      fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
      fvec[3] = x[0] + x[2] - a; 
      errfold = errf;
      errf=0;
      for (k1=0; k1 < 4; k1++)
	{
        if (vr[k1]==0)
	  errf += abs(fvec[k1]);
        else
	  errf += abs(fvec[k1]/vr[k1]);
	}
      if (errf==0)
	break;
      if (errf >= errfold)
	{
	  for (k1=0; k1 < 4; k1++)
	    x[k1] = xold[k1];
	  break;
	}
    }
  AQ=x[0];
  BQ=x[1];
  CQ=x[2];
  DQ=x[3];
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

template <class ntype, int N, class cmplx> void rpoly<ntype,N,cmplx>::oqs_quartic_solver(pvector<cmplx,N>& roots)
{
  /* USAGE:
   *
   * This routine calculates the roots of the quartic equation
   *
   * coeff[4]*x^4 + coeff[3]*x^3 + coeff[2]*x^2 + coeff[1]*x + coeff[0] = 0
   * 
   * if coeff[4] != 0 
   *
   * the four roots will be stored in the complex array roots roots[] 
   *
   * */
  cmplx acx1, bcx1, ccx1, dcx1,acx,bcx,ccx,dcx,cdiskr,zx1,zx2,zxmax,zxmin, qroots[2];
  ntype l2m[12], d2m[12], res[12], resmin, bl311, dml3l3, err0=0, err1=0, aq1, bq1, cq1, dq1; 
  ntype a,b,c,d,phi0,aq,bq,cq,dq,d2,d3,l1,l2,l3, errmin, errv[3], aqv[3], cqv[3],gamma,del2;
  int realcase[2], whichcase, k1, k, kmin, nsol;
  ntype rfactsq, rfact=1.0;

  if (coeff[4]==0.0)
    {
      cout << "That's not a quartic!\n";
      return;
    }

  a=coeff[3]/coeff[4];
  b=coeff[2]/coeff[4];
  c=coeff[1]/coeff[4];
  d=coeff[0]/coeff[4];
  oqs_calc_phi0(a,b,c,d,phi0,0);
  //cout << "phi0=" << phi0 << "\n";
  // simple polynomial rescaling
  if (isnan(phi0)||isinf(phi0))
    {
      rfact = scalfact;
      a /= rfact;
      rfactsq = rfact*rfact;
      b /= rfactsq;
      c /= rfactsq*rfact;
      d /= rfactsq*rfactsq;
      oqs_calc_phi0(a,b,c,d,phi0,1);
    }

  l1=a/2;          /* eq. (4) */                                        
  l3=b/6+phi0/2;   /* eq. (6) */                                
  del2=c-a*l3;     /* defined just after eq. (20) */                             
  nsol=0;
  bl311 =2.*b/3.-phi0-l1*l1;   /* This is d2 as defined in eq. (18)*/ 
  dml3l3 = d-l3*l3;            /* dml3l3 is d3 as defined in eq. (9) with d2=0 */ 
  
  /* Three possible solutions for d2 and l2 (see eqs. (18)-(20) and discussion which follows) */
  if (bl311!=0.0)
    {
      d2m[nsol] = bl311;  
      l2m[nsol] = del2/(2.0*d2m[nsol]);   
      res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }
  if (del2!=0)
    {
      l2m[nsol]=2*dml3l3/del2;
      if (l2m[nsol]!=0)
	{
  	  d2m[nsol]=del2/(2*l2m[nsol]);
	  res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
	  nsol++;
	}

      d2m[nsol] = bl311;
      l2m[nsol] = 2.0*dml3l3/del2;
      res[nsol] = oqs_calc_err_ldlt(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }

  if (nsol==0)
    {
      l2=d2=0.0;
    }
  else
    {
      /* we select the (d2,l2) pair which minimizes errors */
      for (k1=0; k1 < nsol; k1++)
	{
	  if (k1==0 || res[k1] < resmin)
	    {
	      resmin = res[k1];
	      kmin = k1;	
	    }
	}
      d2 = d2m[kmin];
      l2 = l2m[kmin];
    }
  whichcase = 0; 
  if (d2 < 0.0) 
    {
      /* Case I eqs. (27)-(30) */
      gamma=sqrt(-d2);                               
      aq=l1+gamma;                                  
      bq=l3+gamma*l2;                              

      cq=l1-gamma;                                
      dq=l3-gamma*l2;                            
      if(abs(dq) < abs(bq))
	dq=d/bq;                                
      else if(abs(dq) > abs(bq))
	bq=d/dq;                               
      if (abs(aq) < abs(cq))
	{
	  nsol=0;
	  if (dq !=0)
	    {
	      aqv[nsol] = (c - bq*cq)/dq;    /* eq. (37) */
	      errv[nsol]=oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
	      nsol++;
	    }
	  if (cq != 0) 
	    {
	      aqv[nsol] = (b - dq - bq)/cq;  /* eq. (38) */
	      errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
	      nsol++;
	    }
	  aqv[nsol] = a - cq;                /* eq. (39) */
	  errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
	  nsol++;
	  /* we select the value of aq (i.e. alpha1 in the manuscript) which minimizes errors */
	  for (k=0; k < nsol; k++)
	    {
	      if (k==0 || errv[k] < errmin)
		{
		  kmin = k;
		  errmin = errv[k];
		}
	    }
	  aq = aqv[kmin];
	}
      else 
	{
	  nsol = 0;
	  if (bq != 0)
	    { 
	      cqv[nsol] = (c - aq*dq)/bq;              /* eq. (44) */
	      errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
	      nsol++;
	    }
	  if (aq != 0)
	    {
	      cqv[nsol] = (b - bq - dq)/aq;            /* eq. (45) */
	      errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
	      nsol++;
	    }
	  cqv[nsol] = a - aq;                          /*  eq. (46) */
	  errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
	  nsol++;	  
	  /* we select the value of cq (i.e. alpha2 in the manuscript) which minimizes errors */
	  for (k=0; k < nsol; k++)
	    {
	      if (k==0 || errv[k] < errmin)
		{
		  kmin = k;
		  errmin = errv[k];
		}
	    }
	  cq = cqv[kmin];
	}
      realcase[0]=1;
    }
   else if (d2 > 0)   
    {
      /* Case II eqs. (47)-(50) */
      gamma=sqrt(d2); 
      acx=cmplx(l1,gamma);  
      bcx=cmplx(l3,gamma*l2);
      ccx = conj(acx);
      dcx = conj(bcx);
      realcase[0] = 0; 
    }
  else 
    realcase[0] = -1; // d2=0
  /* Case III: d2 is 0 or approximately 0 (in this case check which solution is better) */
  if (realcase[0]==-1 || (abs(d2) <= meps*oqs_max3(abs(2.*b/3.), abs(phi0), l1*l1))) 
    {
      d3 = d - l3*l3;
      if (realcase[0]==1)
	err0 = oqs_calc_err_abcd(a, b, c, d, aq, bq, cq, dq);
      else if (realcase[0]==0)
	err0 = oqs_calc_err_abcd_cmplx(a, b, c, d, acx, bcx, ccx, dcx);
      if (d3 <= 0)
	{
	  realcase[1] = 1;
	  aq1 = l1;   
	  bq1 = l3 + sqrt(-d3);
	  cq1 = l1;
	  dq1 = l3 - sqrt(-d3);
	  if(abs(dq1) < abs(bq1))  
	    dq1=d/bq1;                                        
	  else if(abs(dq1) > abs(bq1))
	    bq1=d/dq1;                                       
	  err1 = oqs_calc_err_abcd(a, b, c, d, aq1, bq1, cq1, dq1); /* eq. (53) */
	}
      else /* complex */
	{
	  realcase[1] = 0;
	  acx1 = l1;
	  bcx1 = cmplx(l3,sqrt(d3));
	  ccx1 = l1;
	  dcx1 = conj(bcx1);
	  err1 = oqs_calc_err_abcd_cmplx(a, b, c, d, acx1, bcx1, ccx1, dcx1); 
	}
      if (realcase[0]==-1 || err1 < err0)
	{
          whichcase=1; // d2 = 0
	  if (realcase[1]==1)
	    {
	      aq = aq1;
	      bq = bq1;
	      cq = cq1;
	      dq = dq1;
	    }
	  else
	    {
	      acx = acx1;
	      bcx = bcx1;
	      ccx = ccx1;
	      dcx = dcx1;
	    }
	}
    }
  if (realcase[whichcase]==1)
    {
      /* if alpha1, beta1, alpha2 and beta2 are real first refine 
       * the coefficient through a Newton-Raphson */
      oqs_NRabcd(a,b,c,d,aq,bq,cq,dq);      
      /* finally calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
      oqs_solve_quadratic(aq,bq,qroots);
      roots[0]=qroots[0];
      roots[1]=qroots[1];        
      oqs_solve_quadratic(cq,dq,qroots);
      roots[2]=qroots[0];
      roots[3]=qroots[1];
    }
  else
    {
      /* complex coefficients of p1 and p2 */
      if (whichcase==0) // d2!=0
	{
	  cdiskr=acx*acx/((cmplx)4.0)-bcx;               
	  /* calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
	  zx1=-acx/((cmplx)2.0)+sqrt(cdiskr);
	  zx2=-acx/((cmplx)2.0)-sqrt(cdiskr);
	  if(abs(zx1) > abs(zx2))
	    zxmax=zx1;
	  else
	    zxmax=zx2;
	  zxmin=bcx/zxmax;        
	  roots[0]=zxmin;
	  roots[1]=conj(zxmin);
	  roots[2]=zxmax;
	  roots[3]=conj(zxmax);
	}
      else // d2 ~ 0
	{
	  /* never gets here! */
	  cdiskr=sqrt(acx*acx-cmplx(4.0)*bcx);
	  zx1 = -cmplx(0.5)*(acx+cdiskr);
	  zx2 = -cmplx(0.5)*(acx-cdiskr);
	  if (abs(zx1) > abs(zx2))
	    zxmax = zx1;
	  else
	    zxmax = zx2;
	  zxmin = bcx/zxmax;
	  roots[0] = zxmax;
	  roots[1] = zxmin;
	  cdiskr=sqrt(ccx*ccx-cmplx(4.0)*dcx);
	  zx1 = -cmplx(0.5)*(ccx+cdiskr);
	  zx2 = -cmplx(0.5)*(ccx-cdiskr);
	  if (abs(zx1) > abs(zx2))
	    zxmax = zx1;
	  else
	    zxmax = zx2;
	  zxmin = dcx/zxmax;
	  roots[2]= zxmax;
	  roots[3]= zxmin;
	}
    }
  if (rfact!=1.0)
    {
      for (k=0; k < 4; k++)
        {
          roots[k] *= rfact;
        }
    }
}
#endif

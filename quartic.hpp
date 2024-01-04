#ifndef _QUARTIC_
#define _QUARTIC_
#include "./pvector.hpp"
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
using namespace std;

template <class ntype, class cmplx> 
class quartic_base_static 
{
  const int n=4;
public:
  constexpr static int dynamic = false;
  pvector<ntype, 5> coeff;
  pvector<ntype, 5> cmon;
  pvector<cmplx, 5> coeffc;
  pvector<cmplx, 5> cmonc;
  int is_cmplx;

  void set_coeff(pvector<ntype,5> v)
    {
      is_cmplx = 0;
      coeff = v;
      cmon[n] = 1.0;
      cmonc[n] = 1.0;
      for (int i=n-1; i >=0; i--)
        {
          cmon[i] = coeff[i]/coeff[n];
          cmonc[i] = cmon[i];
        }
    }
  void set_coeff(pvector<cmplx,5> v)
    {
      is_cmplx = 1;
      coeffc = v;
      cmonc[n] = 1.0;
      for (int i=n-1; i >=0; i--)
        cmonc[i] = coeffc[i]/coeffc[n];
    }

  quartic_base_static()
    {
       // empty
    }
};
template <class ntype, class cmplx> 
class quartic_base_dynamic 
{
  const int n=4;
public:
  constexpr static int dynamic = true;
  pvector<ntype> coeff;
  pvector<ntype> cmon;
  pvector<cmplx> cmonc;
  pvector<cmplx> coeffc;
  int is_cmplx;
  quartic_base_dynamic()
    {
      coeff.allocate(5);
      cmon.allocate(5);
      cmonc.allocate(5);
      coeffc.allocate(5);
    }
  void set_coeff(pvector<ntype,-1> v)
    {
      is_cmplx = 0;
      coeff = v;
      cmon[n] = 1.0;
      cmonc[n] = 1.0;
      for (int i=n-1; i >=0; i--)
        {
          cmon[i] = coeff[i]/coeff[n];
          cmonc[i] = cmon[i];
        }
    }
  void set_coeff(pvector<cmplx,-1> v)
    {
      is_cmplx = 1;
      coeffc = v;
      cmonc[n] = 1.0;
      for (int i=n-1; i >=0; i--)
        cmonc[i] = coeffc[i]/coeffc[n];
    }


  ~quartic_base_dynamic() = default;
};

template <class ntype, class cmplx, bool dynamic> using quarticbase = 
typename std::conditional<(dynamic==false), quartic_base_static <ntype, cmplx>,
	 quartic_base_dynamic <ntype, cmplx>>::type;

template <class ntype, class cmplx=complex<ntype>, bool dynamic=false> 
class quartic: public numeric_limits<ntype>, public quarticbase<ntype,cmplx, dynamic> {
  const int n=4, N=4;
  using quarticbase<ntype,cmplx,dynamic>::coeff;
  using quarticbase<ntype,cmplx,dynamic>::cmon;
  using quarticbase<ntype,cmplx,dynamic>::coeffc;
  using quarticbase<ntype,cmplx,dynamic>::cmonc;
  using quarticbase<ntype,cmplx,dynamic>::is_cmplx;
  using roots_vtype = typename std::conditional<(dynamic==false), pvector<cmplx, 4>,
	 pvector<cmplx, -1>>::type;

  const ntype pigr=acos(ntype(-1.0));
  ntype eps05, meps, maxf, maxf2, maxf3, scalfact, cubic_rescal_fact;
  int maxdigits;
  ntype goaleps;
  const cmplx I = cmplx(0.0,1.0);
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

   ntype oqs_min2(ntype a, ntype b)
    {
      if (a <= b)
	return a;
      else
	return b;
    }

   ntype oqs_min3(ntype a, ntype b, ntype c)
     {
       ntype t;
       t = oqs_min2(a,b);
       return oqs_min2(t,c);
    }

  inline void oqs_quartic_solver(roots_vtype& roots);
  inline void oqs_quartic_solver_cmplx(roots_vtype& roots);      
  inline void oqs_solve_cubic_analytic_depressed_handle_inf(ntype b, ntype c, ntype& sol);
  inline void oqs_solve_cubic_analytic_depressed_handle_inf_cmplx(cmplx b, cmplx c, cmplx& sol);
  inline void oqs_solve_cubic_analytic_depressed(ntype b, ntype c, ntype& sol);
  inline void oqs_solve_cubic_analytic_depressed_cmplx(cmplx b, cmplx c, cmplx& sol);
  inline ntype oqs_calc_phi0(ntype a, ntype b, ntype c, ntype d, ntype& phi0, int scaled);
  inline cmplx oqs_calc_phi0_cmplx(cmplx a, cmplx b, cmplx c, cmplx d, cmplx& phi0, int scaled);
  inline ntype oqs_calc_err_ldlt(ntype b, ntype c, ntype d, ntype d2, ntype l1, ntype l2, ntype l3);
  inline ntype oqs_calc_err_abcd_cmplx(ntype a, ntype b, ntype c, ntype d, 
	    			       cmplx aq, cmplx bq, cmplx cq, cmplx dq);
  inline ntype oqs_calc_err_abcd(ntype a, ntype b, ntype c, ntype d, ntype aq, ntype bq, ntype cq, ntype dq);
  inline ntype oqs_calc_err_abc(ntype a, ntype b, ntype c, ntype aq, ntype bq, ntype cq, ntype dq); 
  inline void NRabcdCCmplx(cmplx a, cmplx b, cmplx c, cmplx d, cmplx& AQ, cmplx& BQ, cmplx& CQ, cmplx& DQ);
  inline void oqs_NRabcd(ntype a, ntype b, ntype c, ntype d, ntype& AQ, ntype& BQ, ntype& CQ, ntype& DQ);
  inline void oqs_solve_quadratic(ntype a, ntype b, cmplx roots[2]);
  ntype oqs_calc_err_ldlt_cmplx(cmplx b, cmplx c, cmplx d, cmplx d2, 
                                 cmplx l1, cmplx l2, cmplx l3)
    {
      /* Eqs. (29) and (30) in the manuscript */
      ntype sum;
      sum =  (b==cmplx(0))?abs(d2 + l1*l1 + cmplx(2.0)*l3):abs(((d2 + l1*l1 + cmplx(2.0)*l3)-b)/b);
      sum += (c==cmplx(0))?abs(cmplx(2.0)*d2*l2 + cmplx(2.0)*l1*l3):abs(((cmplx(2.0)*d2*l2 + cmplx(2.0)*l1*l3)-c)/c);
      sum += (d==cmplx(0))?abs(d2*l2*l2 + l3*l3):abs(((d2*l2*l2 + l3*l3)-d)/d);
      return sum;
    }
  ntype oqs_calc_err_abc_cmplx(cmplx a, cmplx b, cmplx c, cmplx aq, 
                              cmplx bq, cmplx cq, cmplx dq)
    {
      /* Eqs. (48)-(51) in the manuscript */
      ntype sum;
      sum = (c==cmplx(0))?abs(bq*cq + aq*dq):abs(((bq*cq + aq*dq) - c)/c);
      sum +=(b==cmplx(0))?abs(bq + aq*cq + dq):abs(((bq + aq*cq + dq) - b)/b);
      sum +=(a==cmplx(0))?abs(aq + cq):abs(((aq + cq) - a)/a);
      return sum;
    }
  ntype oqs_calc_err_abcd_ccmplx(cmplx a, cmplx b, cmplx c, cmplx d, 
                                 cmplx aq, cmplx bq, cmplx cq, cmplx dq)
    {
      /* Eqs. (68) and (69) in the manuscript */
      ntype sum;
      sum = (d==cmplx(0))?abs(bq*dq):abs((bq*dq-d)/d);
      sum += (c==cmplx(0))?abs(bq*cq + aq*dq):abs(((bq*cq + aq*dq) - c)/c);
      sum +=(b==cmplx(0))?abs(bq + aq*cq + dq):abs(((bq + aq*cq + dq) - b)/b);
      sum +=(a==cmplx(0))?abs(aq + cq):abs(((aq + cq) - a)/a);
      return sum;
    }
  void solve_quadratic(ntype a, ntype b, cmplx roots[2])
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
          bn = cmonc[i] + bn*x;
        }
      
      return bn;
    }
   cmplx evaldpoly(cmplx x)
    {
      // evaluate polynomail via Horner's formula 
      cmplx bn=0.0;
      for (int i=n-1; i >= 0; i--)
        {
          bn = (i+1)*cmonc[i+1] + bn*x;
        }
      return bn;
    }
   ntype evalddpoly(cmplx x)
    {
      // evaluate second derivative of polynomail via Horner's formula 
      ntype bn=0.0;
      if (n == 1)
        return 0;
      for (int i=n-2; i >= 0; i--)
        {
          bn = (i+2)*(i+1)*cmonc[i+2] + bn*x;
        }
      return bn;
    }
  
   ntype evalpoly(ntype x)
    {
      // evaluate polynomail via Horner's formula 
      if (is_cmplx)
        {
          cout << "[evalpoly] for a complex polynomial argument must be a complex number";
          return 0.0;
        }
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
      if (is_cmplx)
        {
          cout << "[evaldpoly] for a complex polynomial argument must be a complex number";
          return 0.0;
        }
   
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
      if (is_cmplx)
        {
          cout << "[evalddpoly] for a complex polynomial argument must be a complex number";
          return 0.0;
        }
      
      if (n == 1)
        return 0;
      for (int i=n-2; i >= 0; i--)
        {
          bn = (i+2)*(i+1)*cmon[i+2] + bn*x;
        }
      return bn;
    }
  inline void find_roots(roots_vtype& roots)
    {
      if (is_cmplx == -1)
        {
          cout << "[ERROR] You have to set the quartic coefficients with set_coeff method first!\n";
          exit(-1);
        }
      else if (is_cmplx == 0)
        {
          oqs_quartic_solver(roots);
        }
      else
        {
          oqs_quartic_solver_cmplx(roots);
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
      // cout << setprecision(50) << "meps=" << meps << "\n";
      eps05 = pow(numeric_limits<ntype>::epsilon(),0.5);
      maxf= getmax();
      //cout << setprecision(50) << "maxf=" << maxf << "\n";
      maxdigits = (numeric_limits<ntype>::digits10)-1;
      //cout << "numeric digits=" << maxdigits << "\n";
      maxf2 = pow(maxf,0.5)/10.0;
      maxf3 = pow(maxf,1.0/3.0)/10.0;
      scalfact = pow(maxf,1.0/4.0)/1.618034;
      cubic_rescal_fact = pow(maxf, 1.0/3.0)/1.618034;
      goaleps=numeric_limits<ntype>::epsilon();   
      is_cmplx=-1;
   }

   quartic() 
    {
      init_const();
    }
};
// quartics with OQS
template <class ntype, class cmplx, bool dynamic> void quartic<ntype,cmplx,dynamic>::oqs_solve_cubic_analytic_depressed_handle_inf(ntype b, ntype c, ntype& sol)
{
 /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
  * where coefficients b and c are large (see sec. 2.2 in the manuscript) */ 
  ntype Q, R, theta, A, B, QR, QRSQ, KK, sqrtQ, RQ;;
  const ntype PI2=pigr/2.0, TWOPI=2.0*pigr;
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
        theta  = pigr;
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

template <class ntype, class cmplx, bool dynamic> void quartic<ntype,cmplx, dynamic>::oqs_solve_cubic_analytic_depressed_handle_inf_cmplx(cmplx b, cmplx c, cmplx& sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * where coefficients b and c are large (see sec. 2.2 in the manuscript) */ 
  cmplx asol[3], Am, Ap, ApB, AmB, Q, R, A, B, QR, QRSQ, KK, RQ;
  const ntype PI2=pigr/2.0, TWOPI=2.0*pigr;
  ntype theta, sqrtQr, Ar, Br, QRr, QRSQr, KKr, RQr;
  const ntype sqrt3=sqrt(3.0)/2.0;
  int arereal=0;
  Q = -b/cmplx(3.0);
  R = cmplx(0.5)*c;
  if (R==cmplx(0))
    {
      sol=sqrt(-b);
      return;
    }
  if (Q.imag()==0 && R.imag()==0)
    {
      arereal=1;
    }
  else
    {
      arereal=0;
    }
  if (arereal)
    {
      if (abs(Q.real()) < abs(R.real()))
        {
          QRr=cmplx(Q/R).real();
          QRSQr=QRr*QRr; 
          KKr=1.0 - Q.real()*QRSQr;
        }
      else
        {
          RQr = cmplx(R/Q).real();
          KKr = copysign(ntype(1.0),ntype(Q.real()))*(RQr*RQr/Q.real()-1.0);
        }

      if (KKr < 0.0)
        {
          sqrtQr=sqrt(Q.real());
          theta = acos((R.real()/abs(Q.real()))/sqrtQr);
          if (theta < PI2) 
            sol = -2.0*sqrtQr*cos(theta/3.0);
          else 
            sol = -2.0*sqrtQr*cos((theta+TWOPI)/3.0);
        }
      else
        {
          if (abs(Q.real()) < abs(R.real()))
            Ar = -copysign(ntype(1.0),ntype(R.real()))*cbrt(abs(R.real())*(1.0+sqrt(KKr)));
          else
            {
              Ar = -copysign(ntype(1.0),ntype(R.real()))*cbrt(abs(R.real())+sqrt(abs(Q.real()))*abs(Q.real())*sqrt(KKr));
            }
          if (Ar==0.0)
            Br=0.0;
          else
            Br = Q.real()/Ar;
          sol = Ar+Br;
        }
    }
  else
    {
      if (abs(Q) < abs(R))
        {
          QR=Q/R;
          QRSQ=QR*QR; 
          KK=cmplx(1.0) - Q*QRSQ;
          Ap = -exp(log(R*(cmplx(1.0)+sqrt(KK)))/cmplx(3.0));
          Am = -exp(log(R*(cmplx(1.0)-sqrt(KK)))/cmplx(3.0));
          if (abs(Ap) > abs(Am))
            A = Ap;
          else
            {
              A = Am;
            }
        }
      else
        {
          RQ = R/Q;
          KK = RQ*RQ/Q-cmplx(1.0);
          KK *= Q*Q*Q;
          Ap = -exp(log(R+sqrt(KK))/cmplx(3.0));
          Am = -exp(log(R-sqrt(KK))/cmplx(3.0));
          if (abs(Ap) > abs(Am))
            A = Ap;
          else
            {
              A = Am;
            }
        }
      if (A==cmplx(0.0))
        B=cmplx(0.0);
      else
        B = Q/A;
      sol = A+B;
      ApB=A+B;
      AmB=A-B;
      asol[0] = ApB; /* this is always largest root even if A=B */
      asol[1] = -cmplx(0.5)*ApB + I*sqrt3*(AmB);
      asol[2] = -cmplx(0.5)*ApB - I*sqrt3*(AmB);
      sol=asol[0];
      if (abs(sol) < abs(asol[1]))
        sol=asol[1];
      if (abs(sol) < abs(asol[2]))
        sol=asol[2];
    }
}
template <class ntype, class cmplx, bool dynamic> void quartic<ntype,cmplx, dynamic>::oqs_solve_cubic_analytic_depressed_cmplx(cmplx b, cmplx c, cmplx& sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c 
   * (see sec. 2.2 in the manuscript) */ 
  cmplx K, Q, R, Q3, R2, A, B, Ap, Am, asol[3], ApB, AmB;
  int arereal=0;
  ntype theta, Q3r, R2r, sqrtQr, Ar, Br;
  const ntype sqrt3=sqrt(3.0)/2.0;
  Q = -b/cmplx(3.0);
  R = cmplx(0.5)*c;
  if (abs(Q) > maxf3 || abs(R) > maxf2)
    {
      oqs_solve_cubic_analytic_depressed_handle_inf_cmplx(b, c, sol);
      return;
    }
  if (Q.imag()==0 && R.imag()==0)
    {
      arereal=1;
      Q3r = cmplx(Sqr(Q)*Q).real();
      R2r = cmplx(Sqr(R)).real();
    }
  else
    {
      arereal=0;
      Q3 = Sqr(Q)*Q;
      R2 = Sqr(R);
    }
  if (arereal)
    {
      if (R2r < Q3r)
        {
          theta = acos(R.real()/sqrt(Q3r));
          sqrtQr=-2.0*sqrt(Q.real());
          if (theta < pigr/2) 
            sol = sqrtQr*cos(theta/3.0);
          else 
            sol = sqrtQr*cos((theta+2.0*pigr)/3.0);
        }
      else
        {
          Ar = -copysign(ntype(1.0),ntype(R.real()))*cbrt(abs(R.real()) + sqrt(R2r - Q3r));
          if (Ar==0.0)
            Br=0.0;
          else
            Br = Q.real()/Ar;
          sol = Ar+Br; /* this is always largest root even if A=B */
        }
    }
  else
    {
      K=sqrt(R2 - Q3);
      Ap = -exp(log(R + K)/cmplx(3.0));
      Am = -exp(log(R - K)/cmplx(3.0));
      //Ap = -pow(R + K,1.0/3.0);
      //Am = -pow(R - K,1.0/3.0);
      if (abs(Ap) > abs(Am))
        A=Ap;
      else 
        A=Am;
      if (A==cmplx(0.0))
        B=0.0;
      else
        B = Q/A;
      ApB=A+B;
      AmB=A-B;
      asol[0] = ApB; /* this is always largest root even if A=B */
      asol[1] = -cmplx(0.5)*ApB + I*sqrt3*(AmB);
      asol[2] = -cmplx(0.5)*ApB - I*sqrt3*(AmB);
      sol=asol[0];
      if (abs(sol) < abs(asol[1]))
        sol=asol[1];
      if (abs(sol) < abs(asol[2]))
        sol=asol[2];
    }
}


template <class ntype, class cmplx, bool dynamic> void quartic<ntype,cmplx,dynamic>::oqs_solve_cubic_analytic_depressed(ntype b, ntype c, ntype& sol)
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
      // if one uses -Ofast rq3 can be > 1 o < -1 
      // causing NaN
      rq3 = R/sqrt(Q3);
      if (rq3 > 1.0)
        theta = 1.0;
      else if (rq3 < -1.0)
        theta  = pigr;
      else
        theta = acos(rq3);
#else
      theta = acos(R/sqrt(Q3));
#endif
      sqrtQ=-2.0*sqrt(Q);
      if (theta < pigr/2.0) 
	sol = sqrtQ*cos(theta/3.0);
      else 
	sol = sqrtQ*cos((theta+2.0*pigr)/3.0);

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
template <class ntype, class cmplx, bool dynamic> cmplx quartic<ntype,cmplx, dynamic>::oqs_calc_phi0_cmplx(cmplx a, cmplx b, cmplx c, cmplx d, cmplx& phi0, int scaled)
{
  /* find phi0 as the dominant root of the depressed and shifted cubic 
   * in eq. (79) (see also the discussion in sec. 2.2 of the manuscript) */
  cmplx rmax, g,h,gg,hh,aq,bq,cq,dq,s,diskr, sp, sm;
  cmplx xxx, gx, x, xold, f, fold, df, xsq;
  cmplx ggss, hhss, dqss, aqs, bqs, cqs;
  ntype rfact, rfactsq;
  ntype maxtt, diskrr;
  int iter;

  /* eq. (87) */ 
  if (a.imag()==0 && b.imag()==0)
    {
      diskrr=cmplx(cmplx(9.0)*a*a-cmplx(24.0)*b).real();                    
      if(diskrr > 0.0)
        { 
          diskrr=sqrt(diskrr);
          if(a.real() > 0.0)
            s=-2.0*b.real()/(3.0*a.real()+diskrr);                     
          else
            s=-2.0*b.real()/(3.0*a.real()-diskrr);                      
        }
      else
        {      
          s=-a.real()/4.0;                                    
        }
    }
  else
    {
      diskr=sqrt(cmplx(9.0)*a*a-cmplx(24.0)*b);                    
      sp = -cmplx(3.0)*a+diskr;
      sm = -cmplx(3.0)*a-diskr;
      if (abs(sp) > abs(sm))
        s = cmplx(2.0)*b/sp;
      else
        s = cmplx(2.0)*b/sm;
    }
  /* eqs. (83) */
  aq=a+cmplx(4.0)*s;                                      
  bq=b+cmplx(3.0)*s*(a+cmplx(2.0)*s);                              
  cq=c+s*(cmplx(2.0)*b+s*(cmplx(3.0)*a+cmplx(4.0)*s));                      
  dq=d+s*(c+s*(b+s*(a+s)));                      
  gg=bq*bq/cmplx(9.0);
  hh=aq*cq;      
  g=hh-cmplx(4.0)*dq-cmplx(3.0)*gg;                       /* eq. (85) */                             
  h=(cmplx(8.0)*dq+hh-cmplx(2.0)*gg)*bq/cmplx(3.0)-cq*cq-dq*aq*aq; /* eq. (86) */         
  oqs_solve_cubic_analytic_depressed_cmplx(g, h, rmax);
  if (isnan(rmax.real()) || isinf(rmax.real())||
      isnan(rmax.real()) || isinf(rmax.real()))

    {
      oqs_solve_cubic_analytic_depressed_handle_inf_cmplx(g, h, rmax);
      if ((isnan(rmax.real()) || isinf(rmax.real())||
           isnan(rmax.real()) || isinf(rmax.real())) && scaled)
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
          ggss=bqs*bqs/cmplx(9.0);
          hhss=aqs*cqs;  
          g=hhss-cmplx(4.0)*dqss-cmplx(3.0)*ggss;                       
          h=(cmplx(8.0)*dqss+hhss-cmplx(2.0)*ggss)*bqs/cmplx(3.0)-cqs*(cqs/rfact)-(dq/rfact)*aqs*aqs; 
          oqs_solve_cubic_analytic_depressed_cmplx(g, h, rmax);
          if (isnan(rmax.real()) || isinf(rmax.real())||
              isnan(rmax.imag()) || isinf(rmax.imag()))
            {
              oqs_solve_cubic_analytic_depressed_handle_inf_cmplx(g, h, rmax);
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
  if (abs(f) > meps*maxtt)
    {
      for (iter=0; iter < 8; iter++)
        {   
          df =  cmplx(3.0)*xsq + g;
          if (df==cmplx(0))
            {
              break;
            }
          xold = x;
          x += -f/df;
          fold = f;
          xsq = x*x;
          f = x*(xsq + g) + h;
          if (f==cmplx(0))
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
  return f;
}
template <class ntype, class cmplx, bool dynamic> ntype quartic<ntype, cmplx, dynamic>::oqs_calc_phi0(ntype a, ntype b, ntype c, ntype d, ntype& phi0, int scaled)
{
  /* find phi0 as the dominant root of the depressed and shifted cubic 
   * in eq. (64) (see also the discussion in sec. 2.2 of the manuscript) */
  ntype rmax, g,h,gg,hh,aq,bq,cq,dq,s,diskr;
  ntype maxtt, xxx, gx, x, xold, f, fold, df, xsq;
  ntype ggss, hhss, dqss, aqs, bqs, cqs, rfact, rfactsq; 
  int iter;

  diskr=9.0*a*a-24.0*b;                    
  /* eq. (67) */
  if(diskr > 0.0)
    { 
      diskr=sqrt(diskr);
      if(a > 0.0)
	s=-2.0*b/(3.0*a+diskr);                     
      else
	s=-2.0*b/(3.0*a-diskr);                      
    }
  else
    {      
      s=-a/4.0;                                    
    }
  /* eqs. (63) */
  aq=a+4.0*s;                                      
  bq=b+3.0*s*(a+2.0*s);                              
  cq=c+s*(2.0*b+s*(3.0*a+4.0*s));                      
  dq=d+s*(c+s*(b+s*(a+s)));                      
  gg=bq*bq/9.0;
  hh=aq*cq;     
  g=hh-4.0*dq-3.0*gg;                       /* eq. (60) */  
  h=(8.0*dq+hh-2.0*gg)*bq/3.0-cq*cq-dq*aq*aq; /* eq. (61) */          
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
	  h=(8.0*dqss+hhss-2.0*ggss)*bqs/3.0-cqs*(cqs/rfact)-(dq/rfact)*aqs*aqs; 
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
  // Note that f is det(M) as in Eq. (22) of my ACM 2020
  return f;
}
template <class ntype, class cmplx, bool dynamic> ntype  quartic<ntype,cmplx, dynamic>::oqs_calc_err_ldlt(ntype b, ntype c, ntype d, ntype d2, ntype l1, ntype l2, ntype l3)
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
template <class ntype, class cmplx, bool dynamic> 
ntype quartic<ntype, cmplx, dynamic>::oqs_calc_err_abcd_cmplx(ntype a, ntype b,  ntype c, ntype d, cmplx aq, 
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
template <class ntype, class cmplx, bool dynamic> ntype quartic<ntype, cmplx, dynamic>::oqs_calc_err_abcd(ntype a, ntype b, ntype c, ntype d, ntype aq, ntype bq, ntype cq, ntype dq)
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
template <class ntype, class cmplx, bool dynamic> ntype  quartic<ntype,cmplx, dynamic>::oqs_calc_err_abc(ntype a, ntype b, ntype c, ntype aq, ntype bq, ntype cq, ntype dq)
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
template <class ntype, class cmplx, bool dynamic> void quartic<ntype,cmplx, dynamic>::oqs_NRabcd(ntype a, ntype b, ntype c, ntype d, ntype& AQ, ntype& BQ, ntype& CQ, ntype& DQ)
{
  /* Newton-Raphson described in sec. 2.3 of the manuscript for complex
   * coefficients a,b,c,d */
  int iter, k1, k2;
  ntype x02, errf, errfold, x[4], dx[4], det, Jinv[4][4], fvec[4], vr[4];
  ntype  errfa, errfmin;
  ntype xmin[4];

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
  errfa=0;
  for (k1=0; k1 < 4; k1++)
    {
      if (vr[k1]==0)
        errf += abs(fvec[k1]);
      else
        errf += abs(fvec[k1]/vr[k1]);
      errfa += abs(fvec[k1]);
      xmin[k1]=x[k1];
    }
  errfmin=errfa;
  if (errfa==0)
    return;
  for (iter = 0; iter < 20; iter++)
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
	{
	  x[k1] += -dx[k1]/det;
	}
      fvec[0] = x[1]*x[3] - d;
      fvec[1] = x[1]*x[2] + x[0]*x[3] - c;
      fvec[2] = x[1] + x[0]*x[2] + x[3] - b;
      fvec[3] = x[0] + x[2] - a; 
      errfold = errf;
      errf=0;
      errfa = 0.0;
      for (k1=0; k1 < 4; k1++)
	{
          if (vr[k1]==0)
            errf += abs(fvec[k1]);
          else
            errf += abs(fvec[k1]/vr[k1]);
          errfa += abs(fvec[k1]);
        }
      if (errfa < errfmin)
        {
          errfmin = errfa;
          for (k1=0; k1 < 4; k1++)
            xmin[k1]=x[k1];
        }
      if (errfa==0)
	break;
      if (errf >= errfold)
        break;
    }
  AQ=xmin[0];
  BQ=xmin[1];
  CQ=xmin[2];
  DQ=xmin[3];
}

template <class ntype, class cmplx, bool dynamic> void quartic<ntype,cmplx, dynamic>::NRabcdCCmplx(cmplx a, cmplx b, cmplx c, cmplx d, 
                  cmplx& AQ, cmplx& BQ, cmplx& CQ, cmplx& DQ)
{
  /* Newton-Raphson described in sec. 2.3 of the manuscript for complex
   * coefficients a,b,c,d */
  int iter, k1, k2;
  cmplx x02, xold[4], dx[4], x[4], det, Jinv[4][4], fvec[4], vr[4];
  ntype errf, errfold;
  ntype errfa, errfmin;
  cmplx xmin[4];
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
  errfa=0;
  for (k1=0; k1 < 4; k1++)
    {
      errf += (vr[k1]==cmplx(0))?abs(fvec[k1]):abs(fvec[k1]/vr[k1]);
      errfa += abs(fvec[k1]);
      xmin[k1]=x[k1];
    }

  errfmin = errfa;
  if (errfa==0)
    return;

  for (iter = 0; iter < 20; iter++)
    {
      x02 = x[0]-x[2];
      det = x[1]*x[1] + x[1]*(-x[2]*x02 - cmplx(2.0)*x[3]) + x[3]*(x[0]*x02 + x[3]);
      if (det==cmplx(0.0))
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
      errfa = 0.0;
      for (k1=0; k1 < 4; k1++)
        {
          errf += (vr[k1]==cmplx(0))?abs(fvec[k1]):abs(fvec[k1]/vr[k1]);
          errfa += abs(fvec[k1]);
        }
      if (errfa < errfmin)
        {
          errfmin = errfa;
          for (k1=0; k1 < 4; k1++)
            xmin[k1]=x[k1];
        }
      if (errfa==0)
        break;
      if (errf >= errfold)
        {
          break;
        }
    }
  AQ=xmin[0];
  BQ=xmin[1];
  CQ=xmin[2];
  DQ=xmin[3];
}

template <class ntype, class cmplx, bool dynamic> void  quartic<ntype,cmplx, dynamic>::oqs_solve_quadratic(ntype a, ntype b, cmplx roots[2])
{ 
  ntype div,sqrtd,diskr,zmax,zmin;
  diskr=a*a-4.0*b;   
  if(diskr>=0.0)
    {
      if(a>=0.0)
	div=-a-sqrt(diskr);
      else
	div=-a+sqrt(diskr);

      zmax=div/2.0;

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
      roots[0]=cmplx(-a/2.0,sqrtd/2.0);
      roots[1]=cmplx(-a/2.0,-sqrtd/2.0);      
    }   
}
template <class ntype, class cmplx, bool dynamic> void quartic<ntype,cmplx, dynamic>::oqs_quartic_solver(roots_vtype& roots)
{
  /* USAGE:
   *
   * This method calculates the roots of the quartic equation
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
  ntype a,b,c,d,phi0,aq,bq,cq,dq,d2,d3,l1,l2,l3, errmin, errv[3], aqv[3], cqv[3],gamma,del2, detM;
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
  detM = oqs_calc_phi0(a,b,c,d,phi0,0);
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
      detM = oqs_calc_phi0(a,b,c,d,phi0,1);
    }
  //cout << setprecision(50) << "phi0=" << phi0 << "\n";
  l1=a/ntype(2.0);          /* eq. (4) */                                        
  l3=b/ntype(6.0)+phi0/ntype(2.0);   /* eq. (6) */                                
  del2=c-a*l3;     /* defined just after eq. (20) */                             
  nsol=0;
  bl311 =ntype(2.)*b/ntype(3.)-phi0-l1*l1;   /* This is d2 as defined in eq. (18)*/ 
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
  //cout << setprecision(50) << "F=" << abs(d2)/(abs(2.*b/3.) + abs(phi0) + l1*l1)/meps << "\n";
  //cout << setprecision(50) << "d2=" << d2 << " detM=" << detM << " detM/d2=" << detM/d2 <<  "\n";
  //ntype minv = oqs_min3(abs(d2*d),abs(d2*d2*l2*l2),abs(l3*l3*d2));

  // CRITICAL FIX: 
  // if d2 == 0 then d3 != 0 and one has to consider 
  // the case 3 (d2 = 0) to find quartic roots (see pags. 9-10 of my ACM 2020).
  // To verify that d2==0 one can use Eq. (2) in my ACM remark, 
  // but in addition one can also check if d3 != 0 by using Eq. (15). 
  // To do this, we note that according to Eq.(15) 
  // d3 = det(M)/d2 = d - d2*l2^2 + l3^2, i.e.
  // det(M) = d2*d - d2^2*l2^2 + d2*l3^2
  // hence to consider d3 != 0 we require that
  // d3 > meps*min{abs(d2*d),abs(d2*d2*l2*l2),abs(l3*l3*d2)     
  
  if (realcase[0]==-1 || abs(d2) <= meps*(abs(2.*b/3.) + abs(phi0) + l1*l1)
    || abs(detM) > meps*oqs_min3(abs(d2*d),abs(d2*d2*l2*l2),abs(l3*l3*d2) )) 
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
template <class ntype, class cmplx, bool dynamic> void quartic<ntype,cmplx, dynamic>::oqs_quartic_solver_cmplx(roots_vtype& roots)      
{
  /* USAGE:
   *
   * This routine calculates the roots of the quartic equation (coeff[] may be complex here)
   *
   * coeff[4]*x^4 + coeff[3]*x^3 + coeff[2]*x^2 + coeff[1]*x + coeff[0] = 0
   * 
   * if coeff[4] != 0 
   *
   * the four roots will be stored in the complex array roots[] 
   *
   * */
  cmplx acx1, bcx1, ccx1, dcx1,acx,bcx,cdiskr,zx1,zx2,zxmax,zxmin, ccx, dcx;
  cmplx l2m[12], d2m[12], bl311, dml3l3; 
  cmplx a,b,c,d,phi0,d2,d3,l1,l2,l3,acxv[3],ccxv[3],gamma,del2,qroots[2], detM;
  ntype res[12], resmin, err0, err1;
  ntype errmin, errv[3];
  int k1, k, kmin, nsol;
  ntype aq, bq, cq, dq;
  ntype rfactsq, rfact=1.0;
  if (coeffc[4]==cmplx(0.0,0.0))
    {
      printf("That's not a quartic!\n");
      return;
    }
  a=coeffc[3]/coeffc[4];
  b=coeffc[2]/coeffc[4];
  c=coeffc[1]/coeffc[4];
  d=coeffc[0]/coeffc[4];
  detM=oqs_calc_phi0_cmplx(a,b,c,d,phi0,0);
  // simple polynomial rescaling
  if (isnan(phi0.real())||isinf(phi0.real())||
      isnan(phi0.imag())||isinf(phi0.imag()))
    {
      rfact = scalfact;
      a /= rfact;
      rfactsq = rfact*rfact;
      b /= rfactsq;
      c /= rfactsq*rfact;
      d /= rfactsq*rfactsq;
      detM=oqs_calc_phi0_cmplx(a,b,c,d,phi0, 1);
    }

  l1=a/cmplx(2);        /* eq. (16) */                                        
  l3=b/cmplx(6)+phi0/cmplx(2); /* eq. (18) */                                
  del2=c-a*l3;   /* defined just after eq. (27) */                               
  nsol=0;
  bl311 =cmplx(2)*b/cmplx(3)-phi0-l1*l1; /* This is d2 as defined in eq. (20)*/ 
  dml3l3 = d-l3*l3;          /* This is d3 as defined in eq. (15) with d2=0 */
  /* Three possible solutions for d2 and l2 (see eqs. (28)) and discussion which follows) */
  if (bl311!=cmplx(0.0))
    {
      d2m[nsol] = bl311;  
      l2m[nsol] = del2/(cmplx(2.0)*d2m[nsol]);   
      res[nsol] = oqs_calc_err_ldlt_cmplx(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }
  if (del2!=cmplx(0))
    {
      l2m[nsol]=cmplx(2)*dml3l3/del2;
      if (l2m[nsol]!=cmplx(0))
        {
          d2m[nsol]=del2/(cmplx(2)*l2m[nsol]);
          res[nsol] = oqs_calc_err_ldlt_cmplx(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
          nsol++;
        }

      d2m[nsol] = bl311;
      l2m[nsol] = cmplx(2.0)*dml3l3/del2;
      res[nsol] = oqs_calc_err_ldlt_cmplx(b,c,d,d2m[nsol], l1, l2m[nsol], l3);
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
  /* Case I eqs. (37)-(40) */
  gamma=sqrt(-d2);                               
  acx=l1+gamma;                                  
  bcx=l3+gamma*l2;                              
  ccx=l1-gamma;                                
  dcx=l3-gamma*l2;                            
  if(abs(dcx) < abs(bcx))
    dcx=d/bcx;     
  else if(abs(dcx) > abs(bcx))
    bcx=d/dcx;    
  if (abs(acx) < abs(ccx))
    {
      nsol=0;
      if (dcx != cmplx(0))
        {
          acxv[nsol] = (c - bcx*ccx)/dcx;   /* see eqs. (47) */
          errv[nsol]=oqs_calc_err_abc_cmplx(a, b, c, acxv[nsol], bcx, ccx, dcx);
          nsol++;
        }
      if (ccx != cmplx(0)) 
        {
          acxv[nsol] = (b - dcx - bcx)/ccx;  /* see eqs. (47) */
          errv[nsol] = oqs_calc_err_abc_cmplx(a, b, c, acxv[nsol], bcx, ccx, dcx);
          nsol++;
        }
      acxv[nsol] = a - ccx;                  /* see eqs. (47) */ 
      errv[nsol] = oqs_calc_err_abc_cmplx(a, b, c, acxv[nsol], bcx, ccx, dcx);
      nsol++;
      /* we select the value of acx (i.e. alpha1 in the manuscript) which minimizes errors */
      for (k=0; k < nsol; k++)
        {
          if (k==0 || errv[k] < errmin)
            {
              kmin = k;
              errmin = errv[k];
            }
        }
      acx = acxv[kmin];
    }
  else 
    {
      nsol = 0;
      if (bcx != cmplx(0))
        { 
          ccxv[nsol] = (c - acx*dcx)/bcx;      /* see eqs. (53) */
          errv[nsol] = oqs_calc_err_abc_cmplx(a, b, c, acx, bcx, ccxv[nsol], dcx);
          nsol++;
        }
      if (acx != cmplx(0))
        {
          ccxv[nsol] = (b - bcx - dcx)/acx;    /* see eqs. (53) */
          errv[nsol] = oqs_calc_err_abc_cmplx(a, b, c, acx, bcx, ccxv[nsol], dcx);
          nsol++;
        }
      ccxv[nsol] = a - acx;                    /* see eqs. (53) */
      errv[nsol] = oqs_calc_err_abc_cmplx(a, b, c, acx, bcx, ccxv[nsol], dcx);
      nsol++;     
      /* we select the value of ccx (i.e. alpha2 in the manuscript) which minimizes errors */
      for (k=0; k < nsol; k++)
        {
          if (k==0 || errv[k] < errmin)
            {
              kmin = k;
              errmin = errv[k];
            }
        }
      ccx = ccxv[kmin];
    }
  /* Case III: d2 is 0 or approximately 0 (in this case check which solution is better) */

  if (abs(d2) <= meps*(abs(cmplx(2.)*b/cmplx(3.)) + abs(phi0) + abs(l1*l1)) 
      || abs(detM) > meps*oqs_min3(abs(d2*d),abs(d2*d2*l2*l2),abs(l3*l3*d2))) 
    {
      d3 = d - l3*l3;
      err0 = oqs_calc_err_abcd_ccmplx(a, b, c, d, acx, bcx, ccx, dcx);
      acx1 = l1;  
      bcx1 = l3 + sqrt(-d3);
      ccx1 = l1;
      dcx1 = l3 - sqrt(-d3);

      if(abs(dcx1) < abs(bcx1)) 
        dcx1=d/bcx1;                                        
      else if(abs(dcx1) > abs(bcx1))
        bcx1=d/dcx1;                                       
      err1 = oqs_calc_err_abcd_ccmplx(a, b, c, d, acx1, bcx1, ccx1, dcx1);
      if (d2==cmplx(0) || err1 < err0)
        {
          acx = acx1;
          bcx = bcx1;
          ccx = ccx1;
          dcx = dcx1;
        }
    }
  if (acx.imag()==0 && bcx.imag()==0 && ccx.imag()==0 && dcx.imag()==0)
    {
      /* if acx, bcx, ccx and dxc are all real do calculations with real numbers... */
      aq=acx.real();
      bq=bcx.real();
      cq=ccx.real();
      dq=dcx.real();
      oqs_NRabcd(a.real(),b.real(),c.real(),d.real(),aq,bq,cq,dq);      
      solve_quadratic(aq,bq,qroots);
      roots[0]=qroots[0];
      roots[1]=qroots[1];        
      solve_quadratic(cq,dq,qroots);
      roots[2]=qroots[0];
      roots[3]=qroots[1];
    }
  else
    {
      /* first refine the coefficient through a Newton-Raphson */
      NRabcdCCmplx(a,b,c,d,acx,bcx,ccx,dcx);
      /* finally calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1) */
      cdiskr=sqrt(acx*acx-cmplx(4.0)*bcx);
      zx1 = -cmplx(0.5)*(acx+cdiskr);
      zx2 = -cmplx(0.5)*(acx-cdiskr);
      if (abs(zx1) > abs(zx2))
        zxmax = zx1;
      else
        zxmax = zx2;
      if (zxmax==cmplx(0))
        zxmin=0;
      else
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
      if (zxmax==cmplx(0))
        zxmin=0;
      else
        zxmin = dcx/zxmax;
      roots[2]= zxmax;
      roots[3]= zxmin;
    }
  if (rfact!=1.0)
    {
      for (k=0; k < 4; k++)
        roots[k] *= rfact;
    }
}
#endif

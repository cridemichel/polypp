#include <stdlib.h>
#include <stdio.h>
#define CPOLY
#ifdef CPOLY
#include "./cpolyvp.hpp"
#else
#include "./rpolyvp.hpp"
#endif
#include<iostream>
#include<fstream>
#include<complex>
#include<list>
#include<string>
#include <iomanip>
using namespace std;
#define WP 1000
#define WPO 50
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl=mpfr_float;
using cmplx=mpc_complex;
using pdbl=mpfr_float;
using pcmplx=mpc_complex;
//template <int N, int digits=200>
//using rpolymp = rpoly<number<mpfr_float_backend<digits>>,N,false,number<mpc_complex_backend<digits>>>;
#define Complex(x,y) cmplx("x","y"))
bool allreal=false, doswap=false;
#undef M_PI
#define M_PI 3.1415926535897932384626433832795029L
//#define PRINTOUT_COEFF
using numty = vldbl;
vldbl *c;
cmplx *er;
void calc_coeff(vldbl c[], cmplx er[]);
numty gauss(void)
{
  numty  a1=3.949846138, a3 = 0.252408784, a5 = 0.076542912, 
    a7 = 0.008355968, a9 = 0.029899776;
  numty sum, r, r2;
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
int NDEG=0;
void calc_coeff_dep_on_case(int CASO)
{
  int i;
  if (CASO==1)
    {
      // wilkinson
      NDEG=10;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      allreal=true;
      for (i=0; i < NDEG; i++)
        { 
          er[i] = i+1;
        }
      calc_coeff(c, er);
    }
  else if (CASO==2)
    {
      // wilkinson
      //
      NDEG = 15;
      allreal=true;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      for (i=0; i < NDEG; i++)
        { 
          er[i] = i+1;
        }
      calc_coeff(c, er);
    }
  else if (CASO==3)
    {
      NDEG = 20;
      allreal=true;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      for (i=0; i < NDEG; i++)
        { 
          er[i] = i+1;
        }
      calc_coeff(c, er);
    }
  else if (CASO==4)
    {
      NDEG = 20;
      allreal=true;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      er[0] = vldbl(-2.1L);
      for (i=1; i < NDEG; i++)
        { 
          er[i] = er[i-1]+vldbl(0.2L);
        }
      calc_coeff(c, er);
    }
  else if (CASO==5)
    {
      NDEG = 10;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      allreal=true;
      for (i=1; i < NDEG+1; i++)
        { 
          er[i-1] = vldbl(1.0L)/vldbl(i);
        }
      calc_coeff(c, er);
    }
  else if (CASO==6)
    {
      NDEG = 15;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      allreal=true;
      for (i=1; i < NDEG+1; i++)
        { 
          er[i-1] = vldbl(1.0L)/vldbl(i);
        }
      calc_coeff(c, er);
    }
  else if (CASO==7)
    {
      NDEG = 20;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      allreal=true;
      for (i=1; i < NDEG+1; i++)
        { 
          er[i-1] = cmplx(1.0)/cmplx(i);
        }
      calc_coeff(c, er);
    }
  else if (CASO==8)
    {
      NDEG = 20;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];
 
      allreal=true;
      for (i=0; i < NDEG; i++)
        { 
          er[i] = vldbl(1.0L)/pow(vldbl(2),NDEG/2-i);
        }
      calc_coeff(c, er);
    }
  else if (CASO==9)
    {
      NDEG = 20;
      allreal=true;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      for (i=0; i < NDEG; i++)
        { 
          er[i] = 1.0L/pow(((vldbl)2),NDEG/2-i)-3.0L;
        }
      calc_coeff(c, er);
    }
  else if (CASO==10)
    {
      NDEG = 20;
      allreal=true;
      doswap=true;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      static complex<long double> erl[20]=

        {-0.98883082622512854506974288293400861L - 
          0.14904226617617444692935471527721756L*1il, \
            -0.98883082622512854506974288293400861L + 
            0.14904226617617444692935471527721756L*1il, \
            -0.90096886790241912623610231950744505L - 
            0.43388373911755812047576833284835875L*1il, \
            -0.90096886790241912623610231950744505L + 
            0.43388373911755812047576833284835875L*1il, \
            -0.73305187182982632852243148927067191L - 
            0.68017273777091939018735870103374024L*1il, \
            -0.73305187182982632852243148927067191L + 
            0.68017273777091939018735870103374024L*1il, \
            -0.50000000000000000000000000000000000L - 
            0.86602540378443864676372317075293618L*1il, \
            -0.50000000000000000000000000000000000L + 
            0.86602540378443864676372317075293618L*1il, \
            -0.22252093395631440428890256449679476L - 
            0.97492791218182360701813168299393122L*1il, \
            -0.22252093395631440428890256449679476L + 
            0.97492791218182360701813168299393122L*1il, 
          0.07473009358642425429093974573476665L - 
            0.99720379718118014822502987087811927L*1il, 
          0.07473009358642425429093974573476665L + 
            0.99720379718118014822502987087811927L*1il, 
          0.36534102436639501454473799892976880L - 
            0.93087374864420425563779924195127531L*1il, 
          0.36534102436639501454473799892976880L + 
            0.93087374864420425563779924195127531L*1il, 
          0.62348980185873353052500488400423981L - 
            0.78183148246802980870844452667405775L*1il, 
          0.62348980185873353052500488400423981L + 
            0.78183148246802980870844452667405775L*1il, 
          0.82623877431599487194516257377267840L - 
            0.56332005806362202774926153802976051L*1il, 
          0.82623877431599487194516257377267840L + 
            0.56332005806362202774926153802976051L*1il, 
          0.95557280578614073281133405376746667L - 
            0.29475517441090421683077298196019097L*1il, 
          0.95557280578614073281133405376746667L + 
            0.29475517441090421683077298196019097L*1il};  
      for (i=0; i < NDEG; i++)
        er[i]=cmplx(erl[i]);

      calc_coeff(c, er);
    }
  else if (CASO==11)
    {
      NDEG = 20;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      static complex<long double> erl[20]=
        {-0.98883082622512854506974288293400861L- 
          0.14904226617617444692935471527721756L*1il, 
          -0.98883082622512854506974288293400861L+ 
            0.14904226617617444692935471527721756L*1il, 
          -0.90096886790241912623610231950744505L- 
            0.43388373911755812047576833284835875L*1il, 
          -0.90096886790241912623610231950744505L+ 
            0.43388373911755812047576833284835875L*1il, 
          -0.73305187182982632852243148927067191L- 
            0.68017273777091939018735870103374024L*1il, 
          -0.73305187182982632852243148927067191L+ 
            0.68017273777091939018735870103374024L*1il, 
          -0.50000000000000000000000000000000000L- 
            0.86602540378443864676372317075293618L*1il, 
          -0.50000000000000000000000000000000000L+ 
            0.86602540378443864676372317075293618L*1il, 
          -0.22252093395631440428890256449679476L- 
            0.97492791218182360701813168299393122L*1il, 
          -0.22252093395631440428890256449679476L+ 
            0.97492791218182360701813168299393122L*1il, 
          0.07473009358642425429093974573476665L- 
            0.99720379718118014822502987087811927L*1il, 
          0.07473009358642425429093974573476665L+ 
            0.99720379718118014822502987087811927L*1il, 
          0.36534102436639501454473799892976880L- 
            0.93087374864420425563779924195127531L*1il, 
          0.36534102436639501454473799892976880L+ 
            0.93087374864420425563779924195127531L*1il, 
          0.62348980185873353052500488400423981L- 
            0.78183148246802980870844452667405775L*1il, 
          0.62348980185873353052500488400423981L+ 
            0.78183148246802980870844452667405775L*1il, 
          0.82623877431599487194516257377267840L- 
            0.56332005806362202774926153802976051L*1il, 
          0.82623877431599487194516257377267840L+ 
            0.56332005806362202774926153802976051L*1il, 
          0.95557280578614073281133405376746667L- 
            0.29475517441090421683077298196019097L*1il, 
          0.95557280578614073281133405376746667L+ 
            0.29475517441090421683077298196019097L*1il};
      doswap=true;
      allreal=true;
      for (i=0; i < NDEG; i++)
        er[i]=cmplx(erl[i]);

      calc_coeff(c, er);
    }
  else if (CASO==12)
    {
      NDEG = 24;
      doswap=true;
      allreal=true;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      static cmplx erl[24];

      erl[0]=cmplx("-3.52E2","0");
      erl[1]=cmplx("-3.52E2", "0");
      erl[2]=cmplx("-2.8371450777E2", "-2.9920517772E2");
      erl[3]=cmplx("-2.8371450777E2", "2.9920517772E2");
      erl[4]=cmplx("-2.7867414048E2", "6.1005469197E2");
      erl[5]=cmplx("-2.7867414048E2", "-6.1005469197E2");
      erl[6]=cmplx("-2.74892372E2", "0");
      erl[7]=cmplx("-2.014171531E2", "0");
      erl[8]=cmplx("-1.255366582E2", "0");
      erl[9]=cmplx("-9.599999999E1", "0");
      erl[10]=cmplx("-8.8692435121E1", "5.5009607430E2");
      erl[11]=cmplx("-8.8692435121E1", "-5.5009607430E2");
      erl[12]=cmplx("-1.6000000000E1", "0");
      erl[13]=cmplx("8.23178509855E1", "0");
      erl[14]=cmplx("8.8692435121E1", "-5.50096074303E2");
      erl[15]=cmplx("8.8692435121E1", "5.50096074303E2");
      erl[16]=cmplx("1.9293739373E2", "1.60865921259E3");
      erl[17]=cmplx("1.9293739373E2", "-1.60865921259E3");
      erl[18]=cmplx("2.0141715312E2", "0");
      erl[19]=cmplx("2.7489237213E2", "0");
      erl[20]=cmplx("7.52E2", "0");
      erl[21]=cmplx("7.52E2", "0");
      erl[22]=cmplx("9.1106065E2", "1.5722");
      erl[23]=cmplx("9.1106065E2", "-1.5722");
      for (i=0; i < NDEG; i++)
        er[i]=erl[i];
      calc_coeff(c, er);
    }
  else if (CASO==13)
    {
      NDEG = 12;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      //roots and coefficients were calculated by Wolfram Mathematica with a precision of 1000 digits
      //print*, "Vanni Noferini's example degree 12 or 35"
      for ( i=0; i < NDEG; i++)
        {
          er[i]=gauss();
        }
      er[0] *= 1E9;
      er[1] *=1E12;

      calc_coeff(c, er);
    }
  else if (CASO==14)
    {
      // Noferini
      NDEG=35;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      //roots and coefficients were calculated by Wolfram Mathematica with a precision of 1000 digits
      //print*, "Vanni Noferini's example degree 12 or 35"
      for ( i=0; i < NDEG; i++)
        {
          er[i]=gauss();
        }
      er[0] *= 1E9;
      er[1] *=1E12;

      calc_coeff(c, er);

    }
  else if (CASO==15)
    {
      NDEG = 10;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      allreal=true;
      er[0]=vldbl(0.1);
      for (i=1; i < NDEG; i++)
        { 
          er[i] = er[i-1]/vldbl(10.0L);
        }
      calc_coeff(c, er);
    }
  else if (CASO==16)
    {
      NDEG = 20;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      allreal=true;
      er[0]=vldbl(0.1);
      for (i=1; i < NDEG; i++)
        { 
          er[i] = er[i-1]/vldbl(10.0L);
        }
      calc_coeff(c, er);
    }
  else if (CASO==17)
    {
      NDEG = 60;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];

      int ii;
      int m=15;
      vldbl pi = 2.0*acos(vldbl(0.0));
      for  (ii=-m+1; ii <= 0; ii++)
        er[m-1+ii]= vldbl("0.9")*exp(cmplx("0",vldbl(ii)*pi/vldbl("2.0")/vldbl(m)));
      for  (ii=1; ii <= m; ii++)
        er[m-1+ii] = conj(er[m-ii]);
      for  (ii=m+1; ii <= 2*m; ii++)
        er[m-1+ii] = exp(cmplx(0,vldbl(ii)*pi/vldbl("2.0")/vldbl(m)));
      for  (ii=2*m+1; ii <= 3*m; ii++)
        er[m-1+ii] = conj(er[5*m-ii]);
      calc_coeff(c, er);
    }
  else if (CASO==18)
    {
      NDEG = 20;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];
      er[0]=1;
      for (i=1; i < NDEG; i++)
        { 
          er[i] = er[i-1]/vldbl("10.0");
        }
      calc_coeff(c, er);
    }
}

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
void sort_sol_opt(cmplx csol[], cmplx exsol[], vldbl allrelerr[])
{
  int k1, k2, k2min=0;
  int perm[NDEG];
  vldbl relerr, relerrmin, relerrmax;
  cmplx diff, solt[NDEG];
  bool used_exsol[NDEG];
  for (k1=0; k1 < NDEG; k1++)
    used_exsol[k1]=false;
  for (k1=0; k1 < NDEG; k1++)
    {
      bool ini = true;
      for (k2=0; k2 < NDEG; k2++)
        {
          if (used_exsol[k2]==true)
            continue;
          diff = csol[k1] - exsol[k2];
          relerr = (exsol[k2]==cmplx(0.0,0.0))?abs(diff):abs(diff/exsol[k2]);
          if (ini==true || relerr <= relerrmin)
           {
             ini=false;
             k2min=k2;
             relerrmin = relerr;
           } 
        }
      perm[k1] = k2min;
      //cout << "perm[" << k1 << "]=" << k2min << "\n";
      allrelerr[k2min] = relerrmin;
      used_exsol[k2min]=true;
    }

  for (k1=0; k1 < NDEG; k1++)
    solt[k1] = csol[k1];

  for (k1=0; k1 < NDEG; k1++)
    csol[perm[k1]] = solt[k1];
}
numty print_accuracy_at(char *str, cmplx csol[], cmplx exsol[], vldbl allrelerr[])
{
  /* we follow FLocke here */
  int k1;
  vldbl relerrmax;
  for (k1=0; k1 < NDEG; k1++)
    {
      if (k1==0 || allrelerr[k1] > relerrmax)
        {
          relerrmax=allrelerr[k1];
        }
    }

  //printf("[%s] relative accuracy=%.16LG\n", str, relerrmax);
  cout << setprecision(WPO) << "[" << str << "relative accuracy=" <<  relerrmax << "]\n";
  return relerrmax;
}

void print_roots(char *str, cmplx er[], cmplx cr[], vldbl allrelerr[])
{
  printf("CASE %s\n", str);
  for (auto i=0; i < NDEG; i++)
    {
      cout << setprecision(WPO) << "root #" << i << " EX: "<< er[i] << " C:" << cr[i];
      cout << setprecision(WPO) << " [ eps: " << allrelerr[i] << " ]\n"; 
    }
}

void print( list<int> l){
    for(list<int>::iterator it=l.begin(); it!=l.end() ; ++it)
            cout << " " << *it;
    cout<<endl;
}


void calc_coeff(vldbl co[], cmplx er[])
{
  vldbl alpha, beta, zero;
  int ii, jj;
  vldbl *rr, *ir, *c;
  rr = new vldbl[NDEG];
  ir = new vldbl[NDEG];
  c = new vldbl[NDEG+1];
  zero = 0.0;
  for (ii=0; ii < NDEG; ii++)
    {
      rr[ii] = er[ii].real();
      ir[ii] = er[ii].imag();
      c[ii]  = 0.0;
    }
  c[NDEG]=1.0;
  ii=0;
  
  while (ii < NDEG)
    { 
      if (ir[ii] == zero) 
        {
          alpha = -rr[ii];
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
      else
        {
          alpha = -rr[ii]*2.0;
          beta = rr[ii]*rr[ii] + ir[ii]*ir[ii];
          for (jj=ii+1; jj >= 0; jj--)
            { 
              //cout << "jj=" << jj << "\n";
              //do jj=ii+1,1,-1
              if (jj == 1)
                {
                  c[jj] = c[jj] + alpha*c[jj-1] + beta;
                }
              else if (jj == 0) 
                {
                  c[jj] = c[jj] + alpha;
                }
              else 
                c[jj] = c[jj] + alpha*c[jj-1] + beta*c[jj-2];
            }
          ii=ii+2;
        }
    }
  for (ii=0; ii < NDEG; ii++)
     co[ii] = c[NDEG-ii-1];
  co[NDEG]=1.0;
#ifdef PRINTOUT_COEFF
  fstream f;
  f.open("coeff.dat", ios::out|ios::trunc);
  for (ii=0; ii < NDEG+1; ii++)
    f << setprecision(WPO) << co[ii] << "\n";
  f.close();
#endif
  delete [] ir;
  delete [] rr;
  delete [] c;
}

int main(int argc, char *argv[])
{
  numty::default_precision(WP);
  cmplx::default_precision(WP);

  char testo2[256];
  int i, CASO;

  if (argc == 2)
    {
      CASO = atoi(argv[1]);
    }
  else
    {
      CASO = 1;
    }
  if (CASO < 1 || CASO > 18)
    {
      printf("Case must be between 1 and 18\n");
      exit(-1);
    }
  calc_coeff_dep_on_case(CASO);

  cout << "NDEG=" << NDEG << "\n";
  pvector<pcmplx> roots(NDEG);
  cmplx* cr = new cmplx[NDEG];
  numty *allrelerr= new numty[NDEG];
  srand48(0);

#ifdef CPOLY
  pvector<pcmplx> ca(NDEG+1);
  for (i=0; i < NDEG+1; i++)
    ca[i]=pcmplx(vldbl(c[i]),0.0);
  cpolyvp<pcmplx,pdbl> rp(NDEG);
#else
  pvector<pdbl> ca(NDEG+1);
  for (i=0; i < NDEG+1; i++)
    ca[i]=pdbl(c[i]);
  rpolyvp<pdbl,pcmplx> rp(NDEG);
#endif
  rp.set_initial_precision(WPO+10);
  rp.set_output_precision(WPO);
  rp.set_coeff(ca);
  rp.find_roots(roots);
  sprintf(testo2, "OPS");
  for (i=0; i < NDEG; i++)
    cr[i] = cmplx(roots[i]);
  // sort roots and calculate relative error
  sort_sol_opt(cr, er, allrelerr);
  print_roots(testo2, er, cr, allrelerr);
  cout << "Forward relative error:\n";
  print_accuracy_at(testo2, cr, er, allrelerr);
  return 0;
}

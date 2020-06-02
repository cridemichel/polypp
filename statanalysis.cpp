#include<stdlib.h>
#include<stdio.h>
#include "./rpoly.hpp"
#include <complex>
#define MPC_MP
//#define GMP_MP
//#define CPP_MP
#ifndef NDEG
#define NDEG 6
#endif
#ifdef CPP_MP
#include <boost/multiprecision/cpp_bin_float.hpp> 
#include <boost/multiprecision/cpp_complex.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl = number<cpp_bin_float<100>>;
using cmplx = cpp_complex<100>;
using pdbl=vldbl;
using pcmplx=cmplx;
#elif defined(GMP_MP)
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl=number<gmp_float<100>>;
using cmplx=number<complex_adaptor<gmp_float<100>>>;
using pdbl=vldbl;
using pcmplx=cmplx;
#elif defined(MPC_MP)
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl=number<mpfr_float_backend<100>>;
using cmplx=number<mpc_complex_backend<100>>;
using pdbl=vldbl;
using pcmplx=cmplx;
#else
using vldbl=long double;
using cmplx=complex<vldbl>;
using pdbl=double;
using pcmplx=complex<pdbl>;
#endif
void calc_coeff(pvector<pdbl,NDEG+1>& co, pvector<cmplx,NDEG> er)
{
  vldbl rr[NDEG], ir[NDEG], c[NDEG+1], alpha, beta, zero;
  int ii, jj;
  /* 
   * il seguente algoritmo per generare i coefficienti vale solo se sono reali */
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
          alpha = -rr[ii]*vldbl(2.0L);
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
     co[ii] = pdbl(c[NDEG-ii-1]);
}
double ranf(void)
{
  return drand48();
}
void sort_sol_opt(pvector<pcmplx,NDEG>& csol, pvector<cmplx,NDEG>& exsol, vldbl allrelerr[])
{
  int k1, k2, k2min;
  int perm[NDEG];
  vldbl relerr, relerrmin;
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
          diff = cmplx(csol[k1]) - exsol[k2];
          relerr = (exsol[k2]==cmplx(0.0,0.0))?abs(diff):abs(diff/exsol[k2]);
          if (ini==true || relerr <= relerrmin)
           {
             ini=false;
             k2min=k2;
             relerrmin = relerr;
           } 
        }

      perm[k1] = k2min;
      allrelerr[k2min] = relerrmin;
      used_exsol[k2min]=true;
    }

  for (k1=0; k1 < NDEG; k1++)
    solt[k1] = cmplx(csol[k1]);

  for (k1=0; k1 < NDEG; k1++)
    csol[perm[k1]] = pcmplx(solt[k1]);
}
#define MAXSOLV 10
#define PEPTS 500
int cmplxreal=0, restart, dojust=-1;
double cumPEall[PEPTS], PEall[PEPTS];
pvector<pcmplx,NDEG> csolall[MAXSOLV];
void print_legend(FILE *f)
{
  int ic;
  if (dojust < 0)
    {
      for (ic=0; ic < 2; ic++)
	fprintf(f,"@    s%d legend \"%s\"\n", ic, ic2algo(ic));
    }
  else
    fprintf(f,"@    s%d legend \"%s\"\n", dojust, ic2algo(dojust));
}
int maxic=3, icref;
char fname[256];
void save_PE(long long int numtrials, int numpts, vldbl dlogdE, vldbl logdEmin)
{
  FILE *f;
  int k, kk, ic;
  sprintf(fname,"PE.dat");
  f = fopen(fname, "w+");
  for (k=0; k < numpts; k++)
    {
      if (PEall[ic][k]==0) 
        continue;
      fprintf(f, "%.32G %.32G\n", double(vldbl(k)*dlogdE+logdEmin), double(PEall[ic][k]/((double)numtrials)/double(NDEG)));
    }
  fclose(f);
  for (k=0; k < numpts; k++)
    {
      cumPEall[k] = 0.0;
      for (kk=k; kk < numpts; kk++) 
        {
          cumPEall[k] += PEall[kk]/((double)numtrials)/((double)NDEG);
        }
    }
  for (ic=0; ic < maxic; ic++)
    {
      if (dojust >= 0 && ic != dojust)
	continue;
      
      sprintf(fname,"cumPE-%s.dat", ic2algo(ic));
      f = fopen(fname, "w+");
      for (k=0; k < numpts; k++)
	{
	  if (cumPEall[ic][k]==0)
	    continue;
	  fprintf(f, "%.32G %.32G\n", double(k*dlogdE+logdEmin), double(cumPEall[ic][k]));
	}
      fclose(f);
    }
}

vldbl allrelerr[NDEG];
int main(int argc, char **argv)
{
  pdbl  sig, sig2, x1, y1;
  vldbl logdE, dlogdE, logdEmax, logdEmin;
  vldbl dE;
  pvector<cmplx,NDEG> exsol;
  pvector<pcmplx,NDEG> csol;
  pvector<pdbl,NDEG+1> co;
  long long int numtrials, its=0, numout, itsI;
  int numpts, ilogdE;
  int k, i;
  pvector<pcmplx,-1> csold(NDEG);
  pvector<pdbl,-1> cod(NDEG+1);
  rpoly<pdbl,-1,pcmplx> oqs(NDEG);
  srand48(4242);
  
  for (k=0; k < PEPTS; k++)
    PEall[k] = 0;

  sig = 1.0;
  sig2= 1.0;
  logdEmax=10.0;
  logdEmin=((int)log10(numeric_limits<vldbl>::epsilon()))-6;
  numpts = PEPTS; 
  dlogdE = (logdEmax -logdEmin)/numpts;

  if (argc>=2)
    numtrials=atoll(argv[1]);
  else 
    numtrials=1000000000;

  restart = 0;
  itsI = 0;
  if (argc>=3)
    numout=atoll(argv[2]);
  else
    numout=100;
  if (argc>=4)
    cmplxreal = atoi(argv[3]);
  if (cmplxreal < 0 || cmplxreal > 5)
    {
      printf("cmplxreal must be between 0 and 5!\n");
      exit(-1);
    }
  if (cmplxreal==3)
    {
      sig = 1.0;
      sig2= 1E6;
      cmplxreal=1;
    }
  else if (cmplxreal==4)
    {
      sig = 1E6;
      sig2 = 1E6;
      cmplxreal = 2;
    } 
  if (argc  >= 5)
    dojust=atoi(argv[4]);
  if (dojust > 2)
    {
      printf("which test should I have to perform?!?\n Last arg is too big!\n");
      exit(-1);
    }
  if (numtrials < 0)
    {
      printf("number of trials must be a positive integer!\n");
      exit(-1);
    }
  else
    {
      printf("numtrials=%lld numout=%lld cmplxreal=%d dojust=%d\n", 
	     numtrials, numout, cmplxreal, dojust);
    }
  //oqs.set_output_prec(1E-40);
  for (its=itsI; its < numtrials; its++)
    {
      if (its > 0 && (its % (numtrials/numout) == 0))
	{
          if (cmplxreal == 0)
	    printf("[SAMPLE A sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
          else if (cmplxreal==1)
	    printf("[SAMPLE B sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
	  else if (cmplxreal==2)
	    printf("[SAMPLE C sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
          else if (cmplxreal==3)
	    printf("[SAMPLE D sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
	  else if (cmplxreal==4)
            printf("[SAMPLE E sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
          else if (cmplxreal==5)
            printf("[SAMPLE F sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
	  save_PE(its, numpts, dlogdE, logdEmin);
	}
      /* generate 4 random roots */
      if (cmplxreal==2) /* all complex */
	{
          for (i=0; i < 2*(NDEG/2); i+=2)
            {
              x1 = sig*(ranf()-0.5);
              y1 = sig*(ranf()-0.5);
              exsol[i] = cmplx(x1, y1);
              exsol[i+1]=cmplx(x1,-y1);
            }
          if (NDEG%2==1)
            {
              exsol[NDEG-1] = sig2*(ranf()-0.5);
            }
        }
      else if (cmplxreal==1) /* half complex half real */
	{
          

          for (i=0; i < NDEG/2; i=i+2)
            { 
              x1 = sig2*(ranf()-0.5);
              y1 = sig2*(ranf()-0.5);
              exsol[i] = cmplx(x1,2.0*y1);
              exsol[i+1] = cmplx(x1,-0.2*y1);
            }

          for (; i < NDEG; i++)
            {
              exsol[i] = cmplx(sig*(ranf()-0.5),0);
            }
	}
      else if (cmplxreal==0)/* all real */
	{
          for (i=0; i < NDEG; i++)
            exsol[i] = cmplx(sig*(ranf()-0.5),0.0);
	}
      if (cmplxreal == 5)
	{
	  co[5]=1.0;
          for (i=0; i < NDEG; i++)
            co[i]= ranf()-0.5;
	}
      else
	{
          calc_coeff(co, exsol);
        }

      
      for (i=0; i <= NDEG; i++)
        cod[i] = co[i];
      oqs.set_coeff(cod);
      oqs.set_polish(false);
      oqs.find_roots(csold);
      for (i=0; i < NDEG; i++)
        csolall[i] = csold[i];
      sort_sol_opt(csolall, exsol, allrelerr);
      for (k=0; k < NDEG; k++)
        {	
          dE = allrelerr[ic][k]; 
          if (dE > 0.0)
            {
              logdE=log10(dE)-logdEmin;
              if (log10(dE) >=10)
                {
                  cout << "exsol[" << k << "]=" << exsol[k] << " csol=" << csolall[k] << "\n";
                  cout << "allrelel=" << allrelerr[k] << "\n";
                  cout << "dE=" << dE << "\n";
                  exit(1);
                }
              ilogdE=(int)(logdE/dlogdE);

              if (ilogdE >= 0 && ilogdE < numpts)
                {
                  (PEall[ilogdE])++;
                }
            }
        }
    }
  save_PE(numtrials, numpts, dlogdE, logdEmin);
  printf("Finished\n");
  exit(0);
}

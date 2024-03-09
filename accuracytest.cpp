#include <stdlib.h>
#include <stdio.h>
#define CPOLY
#ifdef CPOLY
#include "./cpoly.hpp"
#else
#include "./rpoly.hpp"
#endif
#include<complex>
#include<iostream>
#include<fstream>
#include<list>
#include<string>
#include <iomanip>
using namespace std;
#define MPC_MP
//#define GMP_MP
#define WP 200
#define WPP 200
#ifdef CPP_MP
#include <boost/multiprecision/cpp_bin_float.hpp> 
#include <boost/multiprecision/cpp_complex.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl = number<cpp_bin_float<WP>>;
using cmplx = cpp_complex<WP>;
using pdbl=number<cpp_bin_float<WPP>>;
using pcmplx=cpp_complex<WPP>;
#elif defined(GMP_MP)
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl=number<gmp_float<WP>>;
using cmplx=number<complex_adaptor<gmp_float<WP>>>;
using pdbl=number<gmp_float<WPP>>;
using pcmplx=number<complex_adaptor<gmp_float<WPP>>>;
#elif defined(MPC_MP)
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using vldbl=number<mpfr_float_backend<WP>>;
using cmplx=number<mpc_complex_backend<WP>>;
using pdbl=number<mpfr_float_backend<WPP>>;
using pcmplx=number<mpc_complex_backend<WPP>>;
#else
using vldbl=long double;
using cmplx=complex<vldbl>;
using pdbl=double;
using pcmplx=complex<pdbl>;
#define WP 16
#endif
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
  else if (CASO==19)
    {
      NDEG=20;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];
      er[0]=vldbl("1E-15");
      for (i=1; i < NDEG; i++)
        { 
          er[i] = vldbl("1000.0");
        }
      calc_coeff(c, er);
    }
  else if (CASO==20)
    {
      NDEG=252;
      c = new vldbl[NDEG+1];
      er = new cmplx[NDEG];
      c[0]=vldbl("1");
      c[1]=vldbl("0");
      c[2]=vldbl("0");
      c[3]=vldbl("4");
      c[4]=vldbl("10");
      c[5]=vldbl("24");
      c[6]=vldbl("84");
      c[7]=vldbl("287");
      c[8]=vldbl("951");
      c[9]=vldbl("2997");
      c[10]=vldbl("9178");
      c[11]=vldbl("27687");
      c[12]=vldbl("82556");
      c[13]=vldbl("243856");
      c[14]=vldbl("714585");
      c[15]=vldbl("2078280");
      c[16]=vldbl("5998812");
      c[17]=vldbl("17186929");
      c[18]=vldbl("48889915");
      c[19]=vldbl("138114909");
      c[20]=vldbl("387573696");
      c[21]=vldbl("1080521975");
      c[22]=vldbl("2993201892");
      c[23]=vldbl("8239628896");
      c[24]=vldbl("22541699645");
      c[25]=vldbl("61292110635");
      c[26]=vldbl("165648278291");
      c[27]=vldbl("444993617790");
      c[28]=vldbl("1188285124297");
      c[29]=vldbl("3154274100432");
      c[30]=vldbl("8323365625136");
      c[31]=vldbl("21833611738404");
      c[32]=vldbl("56935515935122");
      c[33]=vldbl("147595804188348");
      c[34]=vldbl("380362965100104");
      c[35]=vldbl("974440000395468");
      c[36]=vldbl("2481664860255580");
      c[37]=vldbl("6282905072767168");
      c[38]=vldbl("15812679064764308");
      c[39]=vldbl("39561848766737000");
      c[40]=vldbl("98394972185596224");
      c[41]=vldbl("243273059572206656");
      c[42]=vldbl("597915751348686336");
      c[43]=vldbl("1460869301670406400");
      c[44]=vldbl("3548207176726510080");
      c[45]=vldbl("8567100294413053952");
      c[46]=vldbl("20563028996235599872");
      c[47]=vldbl("49064758617351913472");
      c[48]=vldbl("116381070793653321728");
      c[49]=vldbl("274426190054250119168");
      c[50]=vldbl("643278853540344102912");
      c[51]=vldbl("1499004645324515442688");
      c[52]=vldbl("3472446557055584567296");
      c[53]=vldbl("7996415371966953816064");
      c[54]=vldbl("18305431438069837332480");
      c[55]=vldbl("41656854387498307026944");
      c[56]=vldbl("94234787788137850470400");
      c[57]=vldbl("211909348279228041265152");
      c[58]=vldbl("473694408447087991586816");
      c[59]=vldbl("1052568436777891030630400");
      c[60]=vldbl("2324879849307308480790528");
      c[61]=vldbl("5104373412871347225231360");
      c[62]=vldbl("11139593594723103393447936");
      c[63]=vldbl("24164291696373902792982528");
      c[64]=vldbl("52101361251502338774925312");
      c[65]=vldbl("111657134497566932808499200");
      c[66]=vldbl("237836287050758294069051392");
      c[67]=vldbl("503518607597445749023965184");
      c[68]=vldbl("1059473209511917061794693120");
      c[69]=vldbl("2215605877396960945316560896");
      c[70]=vldbl("4604833195011426891960680448");
      c[71]=vldbl("9511406079677090563654418432");
      c[72]=vldbl("19524292560534994205148184576");
      c[73]=vldbl("39828666195714777692629893120");
      c[74]=vldbl("80741344554785612573215555584");
      c[75]=vldbl("162654321296352351206435192832");
      c[76]=vldbl("325607251868991338477894238208");
      c[77]=vldbl("647695000342969251316757430272");
      c[78]=vldbl("1280218095487171208167548780544");
      c[79]=vldbl("2514334119407674428612822433792");
      c[80]=vldbl("4906547729325344180676568547328");
      c[81]=vldbl("9513300235029901647092291469312");
      c[82]=vldbl("18326365856000555274851628089344");
      c[83]=vldbl("35075108112717762154127788343296");
      c[84]=vldbl("66693950162339302062002502893568");
      c[85]=vldbl("125986945792112789790299961425920");
      c[86]=vldbl("236430200489990805496042138632192");
      c[87]=vldbl("440762802678721830515691455774720");
      c[88]=vldbl("816239227408104183250564289658880");
      c[89]=vldbl("1501502833203250379860951907172352");
      c[90]=vldbl("2743571645769631815459708484976640");
      c[91]=vldbl("4979352399564455799232439573807104");
      c[92]=vldbl("8975950305969163218962282820141056");
      c[93]=vldbl("16070274849360973522304285450174464");
      c[94]=vldbl("28574943948058538856813091060973568");
      c[95]=vldbl("50460351882391898970962784208027648");
      c[96]=vldbl("88491430062273439471679863051517952");
      c[97]=vldbl("154106554106558019948014499436429312");
      c[98]=vldbl("266497125222998597768899131112161280");
      c[99]=vldbl("457612308230538312643860357065998336");
      c[100]=vldbl("780222557140420494375925900657033216");
      c[101]=vldbl("1320798271021936915091765993945432064");
      c[102]=vldbl("2219897253784442962955699750323617792");
      c[103]=vldbl("3704149579364669041625579761466081280");
      c[104]=vldbl("6135980096438938379287269652142489600");
      c[105]=vldbl("10090189738586806503777788386181382144");
      c[106]=vldbl("16470793417469466685589527965842014208");
      c[107]=vldbl("26687567124836071614160683109846089728");
      c[108]=vldbl("42920211312746756367063073247987761152");
      c[109]=vldbl("68509671385517767040534587471330017280");
      c[110]=vldbl("108531918568469332879801996482912976896");
      c[111]=vldbl("170630499708668241530265247767282581504");
      c[112]=vldbl("266211660290453694151517509710686191616");
      c[113]=vldbl("412141140152976250539039674805900017664");
      c[114]=vldbl("633126077448923733342974299403412570112");
      c[115]=vldbl("965019801895808146374559595055626059776");
      c[116]=vldbl("1459352049289286810995420551824451567616");
      c[117]=vldbl("2189461698425262217944576072592950034432");
      c[118]=vldbl("3258691420167757100164258410823560462336");
      c[119]=vldbl("4811189386004281275988190015434767990784");
      c[120]=vldbl("7045945286624261765055759355116886425600");
      c[121]=vldbl("10234755621571758859992587645078728605696");
      c[122]=vldbl("14744851525146654482665753510523336392704");
      c[123]=vldbl("21066911571041400178689631622315315822592");
      c[124]=vldbl("29849097615744058262036918446561781874688");
      c[125]=vldbl("41937565309697396252573434798373806800896");
      c[126]=vldbl("58423581459647120679505253302554293436416");
      c[127]=vldbl("80696897271874605460890884150107972829184");
      c[128]=vldbl("110504353330505868305411554307482685800448");
      c[129]=vldbl("150011812562140627270303805258897352032256");
      c[130]=vldbl("201866431716627522494646539026214617612288");
      c[131]=vldbl("269255014741534253041731641335670951116800");
      c[132]=vldbl("355952799665216325013644146528870396329984");
      c[133]=vldbl("466355609884295839132228271941245091184640");
      c[134]=vldbl("605486989712979058452765860912115675561984");
      c[135]=vldbl("778970924164994874198768646082318909308928");
      c[136]=vldbl("992960232256249032485244768203090012667904");
      c[137]=vldbl("1254010962349185448786573279743186403590144");
      c[138]=vldbl("1568894348038134518783397007325216116834304");
      c[139]=vldbl("1944340313672929438982927208545983286738944");
      c[140]=vldbl("2386710290672834973745175360491547176992768");
      c[141]=vldbl("2901602249256742790644486097724846471380992");
      c[142]=vldbl("3493397243661753161273988713215146648928256");
      c[143]=vldbl("4164764106885675813665124537142252346015744");
      c[144]=vldbl("4916146706938452238935086555467217845092352");
      c[145]=vldbl("5745265688926011629652129432133116255797248");
      c[146]=vldbl("6646673014599326964285276948327598848475136");
      c[147]=vldbl("7611401919759035200650815460341001104654336");
      c[148]=vldbl("8626756193531783862226927284642308356571136");
      c[149]=vldbl("9676280126159358933845113582818615237804032");
      c[150]=vldbl("10739943524440086529045311636179566228144128");
      c[151]=vldbl("11794564704560525876199464137987995656519680");
      c[152]=vldbl("12814478690342377208897566319437749759770624");
      c[153]=vldbl("13772438877482643268574553907663154475696128");
      c[154]=vldbl("14640719623275077521154445624778614321971200");
      c[155]=vldbl("15392366493912262657218958449523639768317952");
      c[156]=vldbl("16002522440229883626433602004348944972251136");
      c[157]=vldbl("16449744216963850380312202468669591058907136");
      c[158]=vldbl("16717215915016695763893227970111814961725440");
      c[159]=vldbl("16793767022807396408179700007422512697704448");
      c[160]=vldbl("16674611678186333431892591228080091524235264");
      c[161]=vldbl("16361743474512921424042848940476933127274496");
      c[162]=vldbl("15863945093374810969001816781866002443927552");
      c[163]=vldbl("15196401970917802932755434351574230792929280");
      c[164]=vldbl("14379941265884801841729224034291086568456192");
      c[165]=vldbl("13439948297574846856541452442242502252036096");
      c[166]=vldbl("12405039087829421916587128860441661927325696");
      c[167]=vldbl("11305586834509163800196259405738392611192832");
      c[168]=vldbl("10172210036085060851789246830147781036867584");
      c[169]=vldbl("9034329629688998439522441905668205270007808");
      c[170]=vldbl("7918892153562551338857855103338447542157312");
      c[171]=vldbl("6849337013239235415200503544369888232472576");
      c[172]=vldbl("5844860785148600951544295841517998299938816");
      c[173]=vldbl("4920003115154670463753477050350448800694272");
      c[174]=vldbl("4084550353907705282245218116270109315039232");
      c[175]=vldbl("3343727594225425263185622084572499968786432");
      c[176]=vldbl("2698629637235855202056382198717471908691968");
      c[177]=vldbl("2146828173097010294873456596265900655509504");
      c[178]=vldbl("1683086715798839331300492859057215394283520");
      c[179]=vldbl("1300116237494349314023738047511302684803072");
      c[180]=vldbl("989311875821176700393870984103251140411392");
      c[181]=vldbl("741422863731658968783929080263566249951232");
      c[182]=vldbl("547122026447575869733711486401822028988416");
      c[183]=vldbl("397455904810095173143348679589044718403584");
      c[184]=vldbl("284170166302694553075186557246369718861824");
      c[185]=vldbl("199916251641761889101855522988826531201024");
      c[186]=vldbl("138353477114239222813333574351865560170496");
      c[187]=vldbl("94165866155856300446720956778414952415232");
      c[188]=vldbl("63015031255763107492509971451159785766912");
      c[189]=vldbl("41449981596817408305079312754403650306048");
      c[190]=vldbl("26792470274130959031999195208700533932032");
      c[191]=vldbl("17013133556865719075593870650932205715456");
      c[192]=vldbl("10609867226949952735093248294399120506880");
      c[193]=vldbl("6496156256555611593367103533932351586304");
      c[194]=vldbl("3903788548934274822671070465559873191936");
      c[195]=vldbl("2301745351369922606679541083287159046144");
      c[196]=vldbl("1331136136949366146707489541650311544832");
      c[197]=vldbl("754796472288282473585429586304579928064");
      c[198]=vldbl("419490453862031049598547679680128876544");
      c[199]=vldbl("228421063935254058058402680632281923584");
      c[200]=vldbl("121815664156547198985350858414009352192");
      c[201]=vldbl("63598640934015290489279417378972631040");
      c[202]=vldbl("32492893522631782775131930881463156736");
      c[203]=vldbl("16238125637316784496455155365131059200");
      c[204]=vldbl("7934032971656751534988919055488909312");
      c[205]=vldbl("3788422497234760307876358665988472832");
      c[206]=vldbl("1766918012141675464728346105533693952");
      c[207]=vldbl("804537498221816915238613554768642048");
      c[208]=vldbl("357451274396562029151565329387749376");
      c[209]=vldbl("154876901626909503677695413748498432");
      c[210]=vldbl("65403861944406301212223947937939456");
      c[211]=vldbl("26903150172348839661387682894839808");
      c[212]=vldbl("10772325932769008965854733824688128");
      c[213]=vldbl("4195973050453187506646315780014080");
      c[214]=vldbl("1588797906223844658361896573861888");
      c[215]=vldbl("584383407679091485472892173418496");
      c[216]=vldbl("208633583083104706882371353313280");
      c[217]=vldbl("72239316450388826421336203067392");
      c[218]=vldbl("24237771243762963775168259817472");
      c[219]=vldbl("7873121699317726657162449518592");
      c[220]=vldbl("2473541944423162821141126447104");
      c[221]=vldbl("750875002985972505460042891264");
      c[222]=vldbl("220000191260424911291603222528");
      c[223]=vldbl("62142246279022240376387796992");
      c[224]=vldbl("16901619952372627042170044416");
      c[225]=vldbl("4420593332141187577132613632");
      c[226]=vldbl("1110293065896595462415187968");
      c[227]=vldbl("267392905433875403384029184");
      c[228]=vldbl("61648022056008524167643136");
      c[229]=vldbl("13582916117819982048395264");
      c[230]=vldbl("2854684982536065129644032");
      c[231]=vldbl("571126909604570217840640");
      c[232]=vldbl("108531611296180597686272");
      c[233]=vldbl("19542488870688730906624");
      c[234]=vldbl("3325460310307654598656");
      c[235]=vldbl("533213823047865270272");
      c[236]=vldbl("80300708407938711552");
      c[237]=vldbl("11317031962306707456");
      c[238]=vldbl("1486526260201046528");
      c[239]=vldbl("181149002760529248");
      c[240]=vldbl("20371911509824784");
      c[241]=vldbl("2101421454112735");
      c[242]=vldbl("197416860698809");
      c[243]=vldbl("16748342652361");
      c[244]=vldbl("1270113848134");
      c[245]=vldbl("85022516111");
      c[246]=vldbl("4944696540");
      c[247]=vldbl("244703088");
      c[248]=vldbl("10016937");
      c[249]=vldbl("325563");
      c[250]=vldbl("7875");
      c[251]=vldbl("126");
      c[252]=vldbl("1");
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
  cout << setprecision(WP) << "[" << str << "relative accuracy=" <<  relerrmax << "\n";
  return relerrmax;
}

void print_roots(char *str, cmplx er[], cmplx cr[], vldbl allrelerr[])
{
  printf("CASE %s\n", str);
  for (auto i=0; i < NDEG; i++)
    {
      cout << setprecision(WP) << "root #" << i << " EX: "<< er[i] << " C:" << cr[i];
      cout << setprecision(WP) << " [ eps: " << allrelerr[i] << " ]\n"; 
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
    f << setprecision(WP) << co[ii] << "\n";
  f.close();
#endif
  delete [] ir;
  delete [] rr;
  delete [] c;
}

int main(int argc, char *argv[])
{
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
  if (CASO < 1 || CASO > 20)
    {
      printf("Case must be between 1 and 19\n");
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
#else
  pvector<pdbl> ca(NDEG+1);
  for (i=0; i < NDEG+1; i++)
    ca[i]=pdbl(c[i]);
#endif
#ifdef CPOLY
  cpoly<pcmplx,-1,pdbl> rp(NDEG);
#else
  rpoly<pdbl,-1,false,pcmplx> rp(NDEG);
#endif
  rp.set_coeff(ca);
  rp.set_calc_errb(true);
  rp.show("poly");
  auto t1=std::chrono::high_resolution_clock::now();
  rp.find_roots(roots);
  auto t2=std::chrono::high_resolution_clock::now();
  std::cout << "finding roots took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
  snprintf(testo2, 256, "OPS");
  for (i=0; i < NDEG; i++)
    cr[i] = cmplx(roots[i]);
  // sort roots and calculate relative error
  pdbl maxrelerr=0.0, relerr;
  for (i=0; i < NDEG; i++)
    {
      if (roots[i]==cmplx(0,0))
        relerr = abs(rp.get_error_bound(i));
      else
        relerr = rp.get_error_bound(i)/abs(roots[i]);
      if (i==0 || relerr > maxrelerr)
        maxrelerr = relerr;
      //cout << "root["<< i << "]=" << roots[i] << " relerr=" << rp.get_error_bound(i)/abs(roots[i]) << "\n";
    }
  cout << "MAXIMUM ESTIMATED RELATIVE ERROR=" << maxrelerr << "\n";
  if (CASO!=20)
    {
      sort_sol_opt(cr, er, allrelerr);
      print_roots(testo2, er, cr, allrelerr);
    }
  else
    {
      int cc=1;
      cout << "roots:\n";
      for (auto v: roots)
        {
          cout << "#" << cc << ": " << setprecision(WP) << real(v); 
          if (imag(v) < 0) 
            cout << " - ";
          else
            cout << " + ";
          cout << setprecision(WP) << abs(imag(v)) << "*I\n";
          cc++;   
        }
    }
  if (CASO!=20)
    {
      cout << "Forward relative error:\n";
      print_accuracy_at(testo2, cr, er, allrelerr);
    }  
  return 0;
}

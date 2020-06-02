#include"./rpoly.hpp"
#define NDEG 20
int main(void)
{
  rpoly<double,NDEG> P;
  pvector<double,NDEG+1> c;
  pvector<complex<double>,NDEG> r;
  // This is the Case 2 among accuracy tests in ACM Trans. Math. Softw. 46, 2, Article 20 (May 2020),
  // https://doi.org/10.1145/3386241
  for (int i=0; i <= NDEG; i++)
    c[i] = 1.0;
  P.set_coeff(c);
  P.find_roots(r);
  r.show("roots");
  int cc=0;
  for (auto& r0: r)
    {
      cout << setprecision(16) << "root #" << cc <<  "=" << r0 << "\n";
      cout << setprecision(16) << "p(#" << cc << ")=" << P.evalpoly(r0) << "\n\n";
      cc++;
    }

  return 0;
}

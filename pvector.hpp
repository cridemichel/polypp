#ifndef _PVECTOR_
#define _PVECTOR_
#include <limits>
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<algorithm>
#include<iterator>
using namespace std;
#define VEC_MOVE_SEMANTIC
#define VEC_LAZY_EVAL
#ifdef USE_LAPACK
#define USE_VEC_LAPACK //N.B. le lapack usano il multithreading che non è compatibile con MOSIX 
#endif
#define VEC_COMMA_INIT
#ifdef USE_VEC_LAPACK
void wrap_dgemv(char ta, double *Y, const double *A, const double *X, int n, double alpha, double beta);
void wrap_sgemv(char ta, float *Y, const float *A, const float *X, int n, float alpha, float beta);
#endif
#include<complex>
#include<cmath>

#ifdef VEC_LAZY_EVAL
// possible types of operation that can be made lazy
struct VOpTypes
{  
  static const int VecPlusVec=0, VecMinusVec=1, VecTimesScal=2, ScalTimesVec=3, VecDivScal=4;
};
//REMARK on Lazyness Logic:
//each class it is an operation of Lhs and Rhs which store the values of the operands
//only when the method get_v is called the operation is performed 
//and this happens only when assign operator of pvector class is called (see below)
// or constructor with AnOpV argument is called
template<typename ntype, int NT, int tipo, typename Lhs, typename Rhs>
struct AnOpV {
    Lhs const& lhs;
    Rhs const& rhs;
    AnOpV(Lhs const& lhs, Rhs const& rhs):lhs(lhs), rhs(rhs) {
      // empty body
      //cout << "qui\n";
    }
    // lazy operations
    inline ntype get_v(int i) const {
      if constexpr (tipo == VOpTypes::VecPlusVec)// vec + vec
	return lhs.get_v(i) + rhs.get_v(i);
      else if constexpr (tipo == VOpTypes::VecMinusVec) // vec - vec
	return lhs.get_v(i) - rhs.get_v(i);
      else if constexpr (tipo == VOpTypes::VecTimesScal) // vec * scalar
	{
	  return lhs.get_v(i)*rhs;
	}
      else if constexpr (tipo == VOpTypes::ScalTimesVec) // scalar * vec 
	{
	  return rhs.get_v(i)*lhs;
	}
      else if constexpr (tipo == VOpTypes::VecDivScal) // scalar div vec 
	{
	  return lhs.get_v(i)/rhs;
	}
      return 0;
#ifdef _I_DO_NOT_WANT_TO_IMPLEMENT_THIS_
      else if constexpr (tipo == 5) // mat * vec 
	{
	  // get two matrices and the multiply them
	  Lhs m1((double (*)[N])lhs.get_mptr()),
	  Rhs m2((double*)lhs.get_vptr());
	  return matvec(Lhs,Rhs)[i];
	}
      else if constexpr (tipo == 6) // vec * mac 
	return 0;
#endif
    }
    int get_N() const {return lhs.get_N();}
    Lhs const& get_lhs() const { return lhs; }
    Rhs const& get_rhs() const { return rhs; }
};
#endif
//NMAX_=45, int MAXSTRA_=40, int NMAXINV_=40, int NMAXMUL_=10, int NSTA_=8, int NLAZY_=5
struct vecpars
{
  // params for my class
  static const int NMAX = 1000; // if N>NMAX use dynamic allocation of vectors
  static const int NSTA = 50;// value below which declares static object in operator funcs (only if m is not dynamically allocated)
  static const int NLAZY = 20;// value below which does not use lazy evaluation for addition
  static const int NMAXMUL = 100; // value above which it uses LAPACK routines for muladd method
};

// template for a static vector
template <class ntype,int NT> class pvecsta {
public:
  ntype v[NT];
  static const int N=NT;
  constexpr static int dynamic = false;
  pvecsta() = default;
  pvecsta(int NN): pvecsta()
    {
      NN=NT;// just to avoid warnings
      //empty body
    }
};
// template for a dynamic vector (i.e. allocated with new)
template <class ntype, int NT> class pvecdyns {
  //int nr, nc;
public:
  constexpr static int dynamic = true;
  static const int N=NT;
  ntype *v; 
  // copy assignment
  pvecdyns<ntype,NT>& operator=(const pvecdyns<ntype,NT>& v1)
    {
      int i;
      for (i=0; i < NT; i++)
	  v[i] = v1.v[i];
      return (*this);
    }
#ifdef VEC_MOVE_SEMANTIC
  // move assignment
  pvecdyns<ntype,NT>& operator=(pvecdyns<ntype,NT>&& v1)
    {
      swap(v, v1.v);
      return (*this);
    }
#endif
  // copy constructor 
  pvecdyns(const pvecdyns<ntype,NT>& v1)
    {
      v = new ntype[NT];
      (*this) = v1;
    }
#ifdef VEC_MOVE_SEMANTIC
  // move constructor
  pvecdyns(pvecdyns<ntype,NT>&& v1)
    {
      (*this).v = v1.v;
      v1.v=nullptr;
    }
#endif
#if 0
  // constructor
  pvecdyns(int NN)
    {
      N=NT;
      v = new ntype[NT];
    }
#endif
  // default constructor
  pvecdyns()
    {
      v = new ntype[NT];
    }
  // destructor
  ~pvecdyns()
    {
      if (v!=nullptr)
        delete[] v;
    }
};	
// template for a dynamic vector using dynamic allocation ptr (i.e. allocated with new)
template <class ntype, int NT> class pvecdynp {
  //int nr, nc;
public:
  constexpr static int dynamic = true;
  ntype *v=nullptr; 
  int N=0;
  bool dealloc=true;
  // copy assignment
  pvecdynp<ntype,NT>& operator=(const pvecdynp<ntype,NT>& v1)
    {
      int i;
      for (i=0; i < N; i++)
	  v[i] = v1.v[i];
      return (*this);
    }
#ifdef VEC_MOVE_SEMANTIC
  // move assignment
  pvecdynp<ntype,NT>& operator=(pvecdynp<ntype,NT>&& v1)
    {
      swap(v, v1.v);
      return (*this);
    }
#endif
  // copy constructor 
  pvecdynp(const pvecdynp<ntype,NT>& v1)
    {
      N=v1.N;
      v = new ntype[N];
      (*this) = v1;
    }
#ifdef VEC_MOVE_SEMANTIC
  // move constructor
  pvecdynp(pvecdynp<ntype,NT>&& v1)
    {
      N=v1.N;
      (*this).v = v1.v;
      v1.v=nullptr;
    }
#endif
  // default constructor
  pvecdynp() = default;
#if 0
    {
      v = new ntype[N];
    }
#endif
  void allocate(int NN)
    {
      N=NN;
      v = new ntype[N];
    }
  void use_vec(int nc, ntype *vv)
    {
      dealloc=false;
      N=nc;
      v=vv;
    }
  void deallocate()
    {
      if (v!=nullptr)
        delete [] v;
    }
  void resize(int NN)
    {
      //pvecdynp<ntype,NT> vt;
      ntype *vt;
      vt = new ntype[NN];
      // initialize with all elements and 0 if larger
      for (auto i=0; i < NN; i++)
	vt[i] = (i < N)?v[i]:0.0;
      if (v!=nullptr)
        delete [] v;
      N=NN;
      v=vt;
      //v = vt;
    }
  pvecdynp(int NN)
    {
      N=NN;
      v = new ntype[N];
    }
  // destructor
  ~pvecdynp()
    {
      if (dealloc==true)
        {
          if (v!=nullptr)
            delete[] v;
        }
    }
};

template <class ntype, int NT> using pvecdyn = 
typename std::conditional<(NT>0), pvecdyns <ntype, NT>,
	 pvecdynp <ntype, NT>>::type;

// dynamic base class
template <class ntype,int NT> class pvecbasedyn: public pvecdyn<ntype,NT> {
  //int nr, nc;
public:
  constexpr static int dynamic = true;
  pvecbasedyn<ntype,NT>(int NN): pvecdyn<ntype,NT>(NN)
  {
    // empty body
  }
  pvecbasedyn<ntype,NT>(): pvecdyn<ntype,NT>()
  {
    // empty body
  }
 void allocate(int NN)
   {
     if constexpr (NT < 0)
       pvecdyn<ntype,NT>::allocate(NN);
   } 
 void deallocate()
   {
     if constexpr (NT < 0)
       pvecdyn<ntype,NT>::deallocate();
   }
 void resize(int NN)
   {
      if constexpr (NT < 0)
       pvecdyn<ntype,NT>::resize(NN);
   }
};
//choose dynamnic or static base class depending on N
template <class ntype,int NT> using pvecbase = typename std::conditional<(NT>vecpars::NMAX||NT<0),pvecbasedyn <ntype, NT>,pvecsta <ntype, NT>>::type;
template <class ntype,int NT=-1> class pvector : public pvecbase<ntype,NT>, vecpars {
  //int vsize;
#ifdef VEC_COMMA_INIT
  unsigned int curidx=0;
#endif
public:
  using pvecbase<ntype,NT>::v;
  using pvecbase<ntype,NT>::N;
  // /numeric_limits<complex<double>>::digits10=0 that's why I set explicitly 16;//fix this!
  int maxdigits=std::numeric_limits<ntype>::digits10-1;
  auto begin()
    {
      // un puntatore è un iteratore!
      return v;
    }

  auto end()
    {
      // un puntatore è un iteratore!
      return v+N;
    }
  void normalize(void)
    {
      ntype invn=1.0/norm();
      *this = invn*(*this);  
    }
  void set_show_digits(int p)
    {
      maxdigits=p;
    }
  pvector<ntype,NT>(int NN): pvecbase<ntype,NT>(NN)
    {
      // empty body
    }
  pvector<ntype,NT>() = default;
  pvector<ntype,NT>(ntype a, ntype b=0, ntype c=0, ntype d=0)
    {
      if (N >= 1)
        v[0]=a;
      if (N>=2)
        v[1]=b;
      if (N>=3)
        v[2]=c;
      if (N>=4)
        v[3]=d;
    }
#ifdef VEC_LAZY_EVAL
  template<typename Lhs, typename Rhs, int tipo>
  // += operator triggers evaluation of lazy expression
  inline pvector<ntype,NT> operator += (AnOpV<ntype,NT, tipo, Lhs, Rhs> const& op)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] += op.get_v[i];
      return *this;
    }

  // -= operator triggers evaluation of lazy expression
  template<typename Lhs, typename Rhs, int tipo>
  inline pvector<ntype,NT> operator -= (AnOpV<ntype,NT, tipo, Lhs, Rhs> const& op)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] -= op.get_v[i];
      return *this;
    }
  // assignment operator triggers evaluation of lazy expression
  template<typename Lhs, typename Rhs, int tipo>
    pvector<ntype,NT>& operator=(AnOpV<ntype,NT, tipo, Lhs, Rhs> const& op) {
      for (int i=0; i < N; i++)
	v[i] = op.get_v(i);
      return (*this);
    }
  // constructor with AnOpV argument triggers evaluation of lazy expression
  template<typename Lhs, typename Rhs, int tipo>
    pvector<ntype,NT>(AnOpV<ntype,NT, tipo, Lhs, Rhs> const& op) {
      for (int i=0; i < N; i++)
	v[i] = op.get_v(i);
    }
  inline ntype get_v(int i) const {return v[i];}
  inline ntype get_vptr() const {return v;}
  inline int get_N() const {return N;}
#endif
#if 0
  int operator!= (const pvector<ntype,NT>& v2)
    {
      for (auto i=0; i < N; i++)
        if (v[i] != v2[i])
          return 1;
      return 0;
    }
  int operator== (const pvector<ntype,NT>& v2)
    {
      for (auto i=0; i < N; i++)
        if (v[i] != v2[i])
          return 0;
      return 1;
    }
#endif

#ifndef VEC_LAZY_EVAL
  // these methods are used if lazyness is not used at all
  inline pvector<ntype,NT> operator +(const pvector<ntype,NT>& param)
    {
      return addition(param);
    }
  inline pvector<ntype,NT> operator -(const pvector<ntype,NT>& param)
    {
      return subtraction(param);
    }
  #endif
  inline pvector<ntype,NT> mulcw(const pvector<ntype,NT>& v1)
    {
      /*  component-wise multiplication, e.g. in 2d 
       *  (ax,ay)*(bx,by) = (ax*bx,ay,by) 
       */
      pvector<ntype,NT> v2;
      if constexpr (NT < 0)
        v2.allocate(N);
      
      for (auto i=0; i < N; i++)
        v2[i] = v[i]*v1.v[i];
      return v2;
    }
  inline pvector<ntype,NT> divcw(const pvector<ntype,NT>& v1)
    {
      /* component-wise division */
      pvector<ntype,NT> v2;
      if constexpr (NT < 0)
        v2.allocate(N);
      
      for (auto i=0; i < N; i++)
        v2[i] = v[i]/v1.v[i];
      return v2;
    }
  // methods used if lazyness is disabled (e.g. for vectors with few elements)
  // addition
  inline pvector<ntype,NT> addition(const pvector<ntype,NT>& v1) const
    {
      if constexpr (NT <= NSTA && pvector<ntype,NT>::dynamic==false)
	{
	  pvector<ntype, NT> v2;
	  int i;
	  if constexpr (NT < 0)
	    v2.allocate(N);
	  for (i=0; i < N; i++)
	    v2.v[i] = v[i] + v1.v[i];
	  return v2;
	}
      else
	{
	  pvector<ntype, NT> v2;
	  int i;
	  if constexpr (NT < 0)
	    v2.allocate(N);

	  for (i=0; i < N; i++)
	    v2.v[i] = v[i] + v1.v[i];

          // N.B. local variable are automatically moved by compiler
	  // if more efficient (c++11) hence std::move is not necessary
	  // and it can prevent compiler optimizations
	  //if (N > NMAX)
	  //  return std::move(v2);
	  //else
          return v2;
	}
    }
  inline pvector<ntype,NT> subtraction(const pvector<ntype,NT>& v1) const
    {
      int i;
      if constexpr (NT <= NSTA && pvector<ntype,NT>::dynamic==false)
	{
	  pvector<ntype, NT> v2;
	  for (i=0; i < N; i++)
	    v2.v[i] = v[i] - v1.v[i];
	  return v2;
	}
      else
	{
	  pvector<ntype, NT> v2;
	  if (NT < 0)
	    v2.allocate(N);
	  for (i=0; i < N; i++)
	    v2.v[i] = v[i] - v1.v[i];
	  // N.B. local variable are automatically moved by compiler
	  // if more efficient (c++11) hence std::move is not necessary
	  // and it can prevent compiler optimizations
	  //if (N > NMAX)
	  //  return std::move(v2);
	  //else
          return v2;
	}
    }
  // *= for pvectors
  inline pvector<ntype,NT> operator *= (const ntype& param)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] *= param;
      return *this;
    }
  // /= for pvectors
  inline pvector<ntype,NT> operator /= (const ntype& param)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] /= param;
      return *this;
    }
  
  // += for pvectors 
  inline pvector<ntype,NT> operator += (const pvector<ntype,NT>& param)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] += param.v[i];
      return *this;
    }
  // -= for pvectors
  inline pvector<ntype,NT> operator -= (const pvector<ntype,NT>& param)
    {
      int i;
      for (i=0; i < N; i++)
	v[i] -= param.v[i];
      return *this;
    }

#ifndef VEC_LAZY_EVAL
  // used if lazyned is completely disabled
  // vec times scalar
  inline pvector<ntype,NT> operator *(const ntype& param) 
    {
      int i;
      pvector<ntype,NT> lv;
      if constexpr (NT < 0)
	lv.allocate(N);
      for (i=0; i < N; i++)
	lv.v[i] = v[i]*param;

      //if (N > NMAX)
      //return std::move(lv);
      //else
      return lv;
    }
#endif
  inline pvector<ntype,NT> vecscal(const ntype& param) const
    {
      int i;
      pvector<ntype,NT> lv;
      if constexpr (NT < 0)
	lv.allocate(N);
      for (i=0; i < N; i++)
	lv.v[i] = v[i]*param;
      //if (N > NMAX)
	//return std::move(lv);
      //else
      return lv;
    }
#ifdef VEC_COMMA_INIT
  /* NOTE: comma operator has the lowest priority among operators 
   * hence it has less operators than <<.
   * In this case << operator initialize first element of matrix
   * and return vector itself, therefore we need to overload comma operator
   * for pvector object returned by << operator et voilà we have the comma initializer!  
   * */
  inline pvector<ntype,NT>& operator,(const ntype& mr)
    {
      curidx++;
      if (curidx >= N)
        {
          cout << "Too many elements in comma initialization!\n";
          exit(-1);
        }
      v[curidx]=mr;
      return (*this);
    }
  inline pvector<ntype,NT>& operator<<(const ntype& mr)
    {
      v[0] = mr;
      curidx=0;
      return (*this);
    } 
#endif

#ifndef VEC_LAZY_EVAL
  friend pvector<ntype,NT> operator *(const ntype& param, const pvector<ntype,NT>& v) 
    {
      int i;
      pvector<ntype,NT> lv;
      if constexpr (NT < 0)
	lv.allocate(N);
   
      for (i=0; i < N; i++)
	lv.v[i] = v.v[i]*param;

      //if (N > NMAX)
      //return std::move(lv);
      //else
      return lv;
    }
#endif
  inline ntype& operator[](const int& i)
    {
      return v[i];
    }
  // scalar product
  inline ntype dot(const pvector<ntype,NT>& v1) const
    {
      int i;
      ntype sum=0;
      for (i=0; i < N; i++)
	{
	  sum+=v[i]*v1.v[i]; 
	}
      return sum;
    }
  // norm 
  inline ntype norm() const
    {
      return sqrt((*this).dot(*this));
    }
  // cross product
  inline pvector<ntype,NT> cross(const pvector<ntype,NT>& v1)
    {
      static pvector<ntype,NT> v2;
      v2.v[0] = v[1]*v1.v[2] - v[2]*v1.v[1];
      v2.v[1] = v[2]*v1.v[0] - v[0]*v1.v[2];
      v2.v[2] = v[0]*v1.v[1] - v[1]*v1.v[0];
      return v2;
    }
  //^ = cross product 
  inline pvector<ntype,NT> operator ^(const pvector<ntype,NT> &vv)
    {
      return (*this).cross(vv);
    }
  inline friend ntype abs(const pvector<ntype,NT> &vv)
    {
      return vv.norm();
    }
  //* dot product
  inline ntype operator *(const pvector<ntype,NT> &vv)
    {
      return (*this).dot(vv);
    }
  bool operator==(const pvector<ntype,NT>& vb)
    {
#if 0
      int i;
      for (i=0; i < N; i++)
        if (v[i]!=vb.v[i])
          {
            return false;
          }
      return true;
#else
      return std::equal(v,v+N, vb.v);
#endif
    }
  // uso un template poiché qua non posso usare pmatrixq al posto di T in quanto ancora non ho dichiarato la classo
  template <typename T>
    void muladd(T& m1, const pvector<ntype,NT>& v1, ntype alpha=1.0, ntype beta=1.0)
    //quando si userà muladd se m1 è pmatrix istanzierà automaticamente questo templated method
    //
      {
	// v=alpha*m1*v1 + beta*v
	//using pvector<ntype,NT>::v;
#ifdef USE_VEC_LAPACK
	if (N > vecpars::NMAXMUL)
	  {
	    if constexpr (std::is_same<ntype, double>::value)
	      {
		wrap_dgemv('t',&(v[0]), &(m1.m[0][0]), &(v1.v[0]), N, alpha, beta);
	      }
	    else if constexpr (std::is_same<ntype, float>::value)
	      {
		wrap_sgemv('t',&(v[0]), &(m1.m[0][0]), &(v1.v[0]), N, alpha, beta);
	      }
	    else
	      {
		if (alpha==1.0 && beta==0)
		  (*this) = m1*v1; 
		else if (alpha==1.0 && beta==1.0)
		  (*this) = m1*v1 + (*this);
		else
		  (*this) = alpha*m1*v1 + beta*(*this);
	      }
	  }
	else
	  {
	    if (alpha==1.0 && beta==0)
		  (*this) = m1*v1; 
		else if (alpha==1.0 && beta==1.0)
		  (*this) = m1*v1 + (*this);
		else
		  (*this) = alpha*m1*v1 + beta*(*this);
	  }
#endif
      }
  ntype ranf(void)
    {
      return drand48();
    }
  void random_box(void)
    {
      for (auto i=0; i < N; i++)
        v[i] = ranf()-0.5;
    }
  void random_orient(void)
    {
      if (N==1)
        {
          v[0] = 1.0;
        }
      else if (N == 2)
        {
          ntype theta=ranf()*2.0*M_PI;
          v[0] = cos(theta);
          v[1] = sin(theta);
        }
      else if (N==3)
        {
          ntype  xisq, xi1, xi2, xi;
          xisq = 1.0;
          while (xisq >= 1.0)
            {
              //cout << "RAND - orient\nRAND - orient\n";
              xi1  = 1.0 - 2.0*drand48();
              xi2  = 1.0 - 2.0*drand48();
              xisq = xi1 * xi1 + xi2 * xi2;
              //cout << setprecision(20) << "xisq=" << xisq << "\n";
            }

          xi = sqrt (abs(1.0 - xisq));
          //cout << setprecision(20) << "xisq=" << xisq << " xi1=" << xi1 << " xi2=" << xi2 << " xi=" << xi << "\n";
          v[0] = 2.0 * xi1 * xi;
          v[1] = 2.0 * xi2 * xi;
          v[2] = 1.0 - 2.0 * xisq;
          //cout << "IN o=" << v[0] << " " << v[1] << " " << v[2] << "\n"; 
        }
      else
        {
          cout << "random orientation for N > 3 not implemented yet\n";
          exit(1);
        }
    }
  pvector<ntype,NT> orto(void)
    {
      pvector<ntype,NT> vt;

      if constexpr (NT < 0)
        {
          vt.allocate(N);
        }
      vt << 0,1;
      if (vt==(*this))
        {
          vt << 1,0;
        }
      vt = vt - (vt*(*this))*(*this);
      vt = (1.0/vt.norm())*vt;
      return vt;
    }         
  void show(void)
    {
      show(NULL);
    }

  void show(const char* str)
    {
      int i;
      if (str!=NULL)
	cout << str;
      if (maxdigits==0)
        maxdigits=32;
      cout << "{";
      for (i=0; i < N; i++)
	{
          cout << setprecision(maxdigits-1) << v[i];
          if (i < N-1)
	    cout << ",";
	}
      cout << "}\n";
    }
  int size()
    {
      return N;
    } 
};
template <class ntype,int NT>
double scalprod(const pvector<ntype,NT>& v1, const pvector<ntype,NT>& v2)
{
  int i;
  double sum=0;
  int N=v1.N;
  for (i=0; i < N; i++)
    {
      sum+=v1.v[i]*v2.v[i]; 
    }
  return sum;
}
#ifdef VEC_LAZY_EVAL
// ADDITION
// AnOpV plus vector 
  template<typename ntype, int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator+(AnOpV<ntype, NT, tipo, Lhs, Rhs> const& lhs, pvector<ntype,NT> const& p)
{
  return AnOpV<ntype, NT, VOpTypes::VecPlusVec, AnOpV<ntype, NT, tipo, Lhs, Rhs>, pvector<ntype,NT>>(lhs, p);
} 
// add expression template with point at the left
template<typename ntype, int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator+(pvector<ntype,NT> const& p, AnOpV<ntype, NT, tipo, Lhs, Rhs> const& rhs) 
{
  return AnOpV< ntype, NT, VOpTypes::VecPlusVec, pvector<ntype,NT>, AnOpV<ntype, NT, tipo, Lhs, Rhs> >(p, rhs);
}

// vector plus vector 
// se N < NSTA restituisce pvector e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce AnOpV e usa la lazy evaluation
template <typename ntype,int NT> 
// conditional restituisce il giusto tipo in base alla condizione (N <= NLAZY)
  typename std::conditional<(NT<=vecpars::NLAZY && NT>=0),pvector<ntype,NT>,AnOpV<ntype, NT, VOpTypes::VecPlusVec, pvector<ntype,NT>, pvector<ntype,NT>>>::type
inline operator+(pvector<ntype,NT> const& lhs, pvector<ntype,NT> const& rhs) 
{
  if constexpr (NT <= vecpars::NLAZY && NT>= 0)
    {
      return lhs.addition(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << vecpars::NSTA << "\n";
      return AnOpV<ntype, NT, VOpTypes::VecPlusVec, pvector<ntype,NT>, pvector<ntype,NT>>(lhs, rhs); // lazy evaluation
    }
}
//AnOpV plus AnOpV
template <typename ntype, int NT, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator+(const AnOpV<ntype, NT, tipoL, LLhs, LRhs> & leftOperandconst, const AnOpV<ntype, NT, tipoR, RLhs, RRhs> & rightOperand)
{
  return  AnOpV<ntype, NT, VOpTypes::VecPlusVec, AnOpV<ntype, NT, tipoL, LLhs, LRhs>, AnOpV<ntype, NT, tipoR, RLhs, RRhs>>(leftOperandconst, rightOperand);
}
// SUBTRACTION
// AnOpV minus vector
  template<typename ntype, int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator-(AnOpV<ntype,NT, tipo, Lhs, Rhs> const& lhs, pvector<ntype,NT> const& p)
{
  return AnOpV<ntype, NT, VOpTypes::VecMinusVec, AnOpV<ntype, NT, tipo, Lhs, Rhs>, pvector<ntype,NT>>(lhs, p);
} 
// vector minus AnOpV 
template<typename ntype, int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator-(pvector<ntype,NT> const& p, AnOpV<ntype, NT, tipo, Lhs, Rhs> const& rhs) 
{
  return AnOpV< ntype, NT, VOpTypes::VecMinusVec, pvector<ntype,NT>, AnOpV<ntype, NT, tipo, Lhs, Rhs> >(p, rhs);
}

// vector minus vector
// se N < NSTA restituisce pvector e usa direttamente il metodo addition di fatto quindi evitando la lazy evaluation 
// altrimenti resituisce AnOpV e usa la lazy evaluation
template <typename ntype,int NT> 
// conditional restituisce il giusto tipo in base alla condizione (N <= NLAZY)
  typename std::conditional<(NT<=vecpars::NLAZY && NT>=0),pvector<ntype,NT>,AnOpV<ntype, NT, VOpTypes::VecMinusVec, pvector<ntype,NT>, pvector<ntype,NT>>>::type
inline operator-(pvector<ntype,NT> const& lhs, pvector<ntype,NT> const& rhs) 
{
  if constexpr (NT <= vecpars::NLAZY && NT >=0)
    {
      return lhs.subtraction(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << NSTA << "\n";
      return AnOpV<ntype, NT, VOpTypes::VecMinusVec, pvector<ntype,NT>, pvector<ntype,NT>>(lhs, rhs); // lazy evaluation
    }
}
// AnOpV time AnOpV
  template <typename ntype, int NT, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline auto operator-(const AnOpV<ntype, NT, tipoL, LLhs, LRhs> & leftOperandconst, const AnOpV<ntype, NT, tipoR, RLhs, RRhs> & rightOperand)
{
  return  AnOpV<ntype, NT, VOpTypes::VecMinusVec, AnOpV<ntype, NT, tipoL, LLhs, LRhs>, AnOpV<ntype, NT, tipoR, RLhs, RRhs>>(leftOperandconst, rightOperand);
}
   
///
// Scalar product of two AnOpV
//
  template <typename ntype, int NT, int tipoL, int tipoR, typename LLhs, typename LRhs, typename RLhs, typename RRhs>
inline ntype operator*(const AnOpV<ntype, NT, tipoL, LLhs, LRhs> & lhs, const AnOpV<ntype, NT, tipoR, RLhs, RRhs> & rhs)
{
  pvector<ntype, NT> vR, vL;
  if constexpr (NT < 0)
    {
      int N=lhs.get_N();
      vR.allocate(N);
      vL.allocate(N);
    }
  for (int i=0; i < vR.N; i++)
    {
      vL[i] = lhs.get_v(i);
      vR[i] = rhs.get_v(i);
    }
  return vL*vR;
}
// AnOpV times pvector
//
   template<typename ntype, int NT, int tipo, typename Lhs, typename Rhs> ntype 
inline operator*(AnOpV<ntype,NT, tipo, Lhs, Rhs> const& lhs, pvector<ntype,NT> const& vR)
{
  pvector<ntype, NT> vL;
  if constexpr (NT < 0)
    vL.allocate(vR.N);
   
  for (int i=0; i < vR.N; i++)
    vL[i] = lhs.get_v(i);
  return vL*vR;
} 
// vector times right AnOpV
template<typename ntype, int NT, int tipo, typename Lhs, typename Rhs>  
inline ntype operator*(pvector<ntype,NT> const& vL, AnOpV<ntype, NT, tipo, Lhs, Rhs> const& rhs) 
{
  pvector<ntype, NT> vR;
  if constexpr (NT < 0)
    vR.allocate(vL.N);
  for (int i=0; i < vL.N; i++)
    {
      vR[i] = rhs.get_v(i);
    }
  return vL*vR;
}
// scalar times AnOvV 
template<typename ntype, int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator*(const ntype& lhs, AnOpV<ntype, NT, tipo, Lhs, Rhs> const& rhs)
{
  return AnOpV<ntype, NT, VOpTypes::ScalTimesVec, ntype, AnOpV<ntype, NT, tipo, Lhs, Rhs>>(lhs, rhs);
}
// scalar times vector
template<typename ntype,int NT>
  typename std::conditional<(NT<=vecpars::NLAZY && NT>=0),pvector<ntype,NT>,AnOpV<ntype,NT, VOpTypes::ScalTimesVec, ntype, pvector<ntype,NT>>>::type
inline operator*(const ntype& lhs, pvector<ntype,NT> const& rhs)
{
  if constexpr (NT <= vecpars::NLAZY && NT >= 0)
    {
      return rhs.vecscal(lhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << vecpars::NSTA << "\n";
      return AnOpV<ntype, NT, VOpTypes::ScalTimesVec, ntype, pvector<ntype,NT>>(lhs, rhs);
    }
}

// AnOpV times scalar
template<typename ntype, int NT, int tipo, typename Lhs, typename Rhs> 
inline auto operator*(AnOpV<ntype, NT, tipo, Lhs, Rhs> const& lhs, const ntype& p)
{
  return AnOpV<ntype, NT, VOpTypes::VecTimesScal, AnOpV<ntype, NT, tipo, Lhs, Rhs>, ntype>(lhs, p);
}
// PmOvV times scalar
template<typename ntype,int NT>
  typename std::conditional<(NT<=vecpars::NLAZY&& NT>=0),pvector<ntype,NT>,AnOpV<ntype,NT, VOpTypes::VecTimesScal, pvector<ntype,NT>, ntype>>::type
inline operator*(pvector<ntype,NT> const& lhs, const ntype& rhs)
{
  if constexpr (NT <= vecpars::NLAZY && NT >=0 )
    {
      return lhs.vecscal(rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << vecpars::NSTA << "\n";
      return AnOpV<ntype, NT, VOpTypes::VecTimesScal, pvector<ntype,NT>, ntype>(lhs, rhs);
    }
}
// vector divided by scalar 
  template <typename ntype,int NT>
typename std::conditional<(NT<=vecpars::NLAZY&&NT>0),pvector<ntype,NT>,AnOpV<ntype,NT, VOpTypes::VecDivScal, pvector<ntype,NT>, ntype>>::type
inline operator /(pvector<ntype,NT> const& lhs, const ntype& rhs) 
{
  if constexpr (NT <= vecpars::NLAZY && NT >= 0)
    {
      return lhs.vecscal(1.0/rhs); //normal evaluation
    }
  else
    {
      //cout << "qui N=" << N << " NSTA=" << vecpars::NSTA << "\n";
      return AnOpV<ntype, NT, VOpTypes::VecDivScal, pvector<ntype,NT>, ntype>(lhs, rhs);
    }
}

// AnOpV divided by scalar
  template <typename ntype, int NT, int tipo, typename Lhs, typename Rhs>
inline auto operator /(AnOpV<ntype, NT, tipo, Lhs, Rhs> const& lhs, const ntype& rhs) 
{
  return AnOpV<ntype, NT, VOpTypes::VecDivScal, AnOpV<ntype, NT, tipo, Lhs, Rhs>, ntype>(lhs, rhs);
}
#endif

/* some predefined vectors useful for simulations */
typedef pvector<double,2> vecd2;
typedef pvector<double,3> vecd3;
typedef pvector<double,4> vecd4;
#endif

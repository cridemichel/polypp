#include "./linked_cell_lists_2d.hpp"
template <class ntype>
class part
{
public:
  ntype x[2];
  ntype r;
  vector<int> bonds;
};
template <class cmplx, class ntype>
class clusters
{
  linked_cell_lists_2d<part<ntype>> ll;
  vector<part<ntype>> parts;
  int Np;
  vector<int>* listneigh;
  ntype maxr; //maximum radius of particles
  pvector<ntype,2> L; // box size 
  pvector<ntype,2> rCM[2];
public:

  void init(pvector<cmplx>& ro, pvector<ntype>&rad, ntype& maxerr)
    {
      Np = (*ro).size();
      parts.resize(Np);
      maxr = maxerr;
      for (int i=0; i < Np; i++)
        {
          parts[i].x[0] = real(ro[i]);
          parts[i].x[1] = imag(ro[i]);
          parts[i].r = rad[i];
        }
    }
  using vvint = vector<vector<int>>;
  void color_algo(vvint& cls)
    {
      /* find clusters by using color algorithm */ 
    
    }
  void find_clusters(vvint& cls)
    {
      create_ll();
      find_bonds();
      color_algo(cls);
    }
  void find_bonds(void)
    {
      int i; 
      for (i=0; i < Np; i++)
        parts[i].bonds.clear();
      /* we need to employ linked cell lists to find bonds
       * otherwise can become very slow for large n */
      for (i=0; i < Np; i++)
        {
          ll.findneigh(i);
          listneigh = &(ll.listneigh);
          for (int j: *listneigh)
            {
              auto sig = parts[i].r + parts[j].r;
              if (abs(parts[i].z - parts[j].z) <= sig)
                {
                  parts[i].bonds.add(j);
                  parts[j].bonds.add(i);
                }
            }
        }
    }
  void calculate_L(void)
    {
      bool first=true;
      ntype absx[2];
      for (auto p: parts)
        {
          absx[0] = abs(p.x[0])+abs(p.r);
          absx[1] = abs(p.x[1])+abs(p.r);
          if (first)
            {
              first = false;
              L[0] = absx[0];
              L[1] = absx[1];
            }
          else
            {
              if (absx[0] > L[0])
                L[0] = absx[0];
              if (absx[1] > L[1])
                L[1] = absx[1];
            }
        }
      L[0] *= 2.0;
    }
  void translate_cm(void)
    {
      rCM << ntype(0.0), ntype(0.0);
      for (auto& p: parts)
        {
          rCM[0] += p.x[0];
          rCM[1] += p.x[1];
        }
      rCM /= ntype(Np);
      for (auto& p: parts)
        {
          p.x[0] -= rCM[0];
          p.x[1] -= rCM[1];
        }
    }
  void create_ll(void)
    {
      translate_cm();
      /* calculate L and cutoff radius */
      calculate_L();
      ll.init(&parts,L,ntype(2.0)*maxr,Np);
      ll.build();
    }
  clusters()
    {}

  ~clusters()
    {}

};


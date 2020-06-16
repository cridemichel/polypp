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

public:

  void init(pvector<cmplx>& ro, pvector<ntype>&rad )
    {
      Np = (*ro).size();
      parts.resize(Np);

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
  void create_ll(void)
    {
      ntype L[2], rc;
      /* calculate L and cutoff radius */
      ll.init(&parts,L,rc,Np);
      ll.build();
    }
  clusters()
    {}

  ~clusters()
    {}

};


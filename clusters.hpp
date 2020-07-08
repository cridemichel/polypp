#include "./linked_cell_lists_2d.hpp"
template <class ntype>
class part
{
public:
  using numtype=ntype;
  pvector<ntype,2> x;
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
  vector<int> colors, clsdim, clscol;
public:

  void init(pvector<cmplx>& ro, pvector<ntype>& rad, ntype& maxerr)
    {
      Np = ro.size();
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

  void change_all_colors(int colorsrc, int colordst)
    {
      int ii;
      for (ii = 0; ii < Np; ii++)
        {
          if (colors[ii] == colorsrc)
            colors[ii] = colordst;
        }
    }

  int findmaxColor(void)
    {
      int i, maxc=-1;
      for (i = 0; i < Np; i++) 
        {
          if (colors[i] > maxc)
            maxc = colors[i];
        }
      return maxc;
    }

  void color_algo(vvint& cls)
    {
      int curcolor=0, i, ncls, nc;
      /* find clusters by using color algorithm */ 
      colors.resize(Np);
      clsdim.resize(Np);
      clscol.resize(Np);
      for (i=0; i < Np; i++)
        colors[i] = -1;
      for (i=0; i < Np; i++)
        {
          if (colors[i] == -1)
            colors[i] = curcolor;
          
          for (int j: parts[i].bonds)
            {
              if (colors[j] == -1)
                colors[j] = colors[i];
              else
                {
                  if (colors[i] < colors[j])
                    change_all_colors(colors[j], colors[i]);
                  else if (colors[i] > colors[j])
                    change_all_colors(colors[i], colors[j]);
                }
            }
          curcolor = findmaxColor()+1;
        }
      /* single particles are clusters of size 1 */
      for (i = 0; i < Np; i++)
	{
	  if (colors[i]==-1)
	    {	    
	      colors[i] = curcolor;
	      curcolor++;
	    } 
	}
      ncls = curcolor;
      for (nc = 0; nc < ncls; nc++)
	{
	  clsdim[nc] = 0; 
	}
      int a;
      for (nc = 0; nc < ncls; nc++)
	{
	  for (a = 0; a < Np; a++)
	    if (colors[a] == nc)
	      {
		clsdim[nc]++;
		clscol[nc] = colors[a];
	      }
	}
      cls.resize(ncls);
      int clsidx=0;
      for (nc = 0; nc < ncls; nc++)
        {
          if (clsdim[nc] == 0)
            continue;
          cls[clsidx].resize(clsdim[nc]);
          for (i=0; i < Np; i++)
            {
              if (colors[i] == clscol[nc])
                cls[clsidx].push_back(i);
            }
          clsidx++;
        }
      cls.resize(clsidx);
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
              if ((parts[i].x - parts[j].x).norm() <= sig)
                {
                  parts[i].bonds.push_back(j);
                  parts[j].bonds.push_back(i);
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


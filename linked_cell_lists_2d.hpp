#ifndef _LINKED_CELL_LISTS_
#define _LINKED_CELL_LISTS_
#include<vector>
#include<list>
#include "./pvector.hpp"
/* 
 * N.B. linked cell list must be used if M >= 3 otherwise particles pairs are conisdered more than once
 * giving an energy greater than the real one. Indeed if M==2 LL are meaningless since we have all vs all 
 * If M==2 use NNL!!!
 * */
using namespace std;
// particles is the particle type, e.g. lj<double> 
template <class particles, class rectangles=int, bool boxes=false>
class linked_cell_lists_2d
{
  using ntype = typename particles::numtype;
  vector<particles>* part;
  pvector<int,-1> cellList, cellListBak;
  int cellsx, cellsy, cellsz;
  pvector<int,-1> inCell[2];
  int parnum;
  pvector<ntype,2> L, L2;// box size
  pvector<ntype,2> Lbak, L2bak;
  int cellsxBak, cellsyBak;
  int parnumBak;
  pvector<int,-1> inCellBak[2];
  int cxini, cyini;
  ntype rcut;
   
  void insert_in_new_cell(int i)
    {
      int n;
      n = inCell[1][i]*cellsx + inCell[0][i] + parnum;
      //if (i==332)
      //  cout << "<<<i=332" << " n=" << n-parnum << "\n";
      cellList[i] = cellList[n];
      cellList[n] = i;

    }
  void remove_from_current_cell(int i)
    {
      int n;
      n = inCell[1][i]*cellsx + inCell[0][i] + parnum;
      while (cellList[n] != i) 
        n = cellList[n];
      cellList[n] = cellList[i];
    }
  // method for GC simulations (see also below)
  void adjLinkedListRemove(int ip)
    {
      int k;
      // remove parnum 
      for (k = parnum; k < cellsx*cellsy + parnum; k++)
        {
          cellList[k] = cellList[k+1];
        }
    }
  // method for GC simulations (see also below)
  void adjLinkedListInsert(void)
    {
      int k;
      for (k = cellsx*cellsy + parnum - 1; k >= parnum; k--)
        {
          cellList[k] = cellList[k-1];
        }
    }
  void add_to_LL(int n, bool apply_pbc=false)
    {
      int cx, cy;
      //cout << "[add_to_LL] part " << n << " =" << (*part)[n].r[0] << " " << (*part)[n].r[1] << "\n";
      pvector<ntype,3> pos;

      pos = (*part)[n].r;
      if (apply_pbc)
        apply_pbc_to_par(pos);
      cx =  (pos[0] + 0.5*L[0]) * cellsx / L[0];
      cy =  (pos[1] + 0.5*L[1]) * cellsy / L[1];
      inCell[0][n] = cx;
      inCell[1][n] = cy;
      insert_in_new_cell(n);
    }
  
public:

 void too_few_cells(void)
    {
      // N.B. if cells are less than 3 the shift must be calculated
      // considering the distance of two particles and *not*
      // by considering the cell the particle belongs to as in 
      // rapaport method which I actually use in find_neigh method.
      // This check is used when calculating energy to use the distance
      // method for calculaing the shift.
      // Consider that energy calculations are wrong with less than 3 cells per size */ 
      if (cellsx < 3 || cellsy < 3)
        {
          cout << "[init LL] Number of linked cell lists along each direction must be >= 3\n";
          cout << "Restart the simulation disabling them\n";
          exit(1);
        }
    } 
  ntype get_rcut(void)
    {
      return rcut;
    }
  int numneigh;
  vector<int> listneigh;
  //pvector<int,-1> listneigh;
  // method for GC simulations
  void insert_par(bool apply_pbc=false)
    {
      resize_vectors(parnum+1);
      adjLinkedListInsert();
      add_to_LL(parnum-1, apply_pbc); //last particle 
    }
  // method for GC simulations
  void remove_par(int ip)
    {
      int np;
      /* last particle parnum-1 becomes ip-th */
      remove_from_current_cell(ip);
      np=parnum-1;
      //cout << "ip=" << ip << " np=" << np << "parnum=" << parnum << "\n";
      if (ip < np)
        remove_from_current_cell(np);
      parnum--; 
      if (ip < np)
        {
          for (auto k=0; k < 2; k++)
            inCell[k][ip] = inCell[k][np];
        }
      adjLinkedListRemove(ip);
      if (ip < np)
        insert_in_new_cell(ip);
      //show_lists(); 
      resize_vectors(parnum);
    } 
  void show_lists(void)
    {
      store();
      build();
      for (auto n=parnum; n < cellsx*cellsy*cellsz+parnum; n++)
        {
          cout << "A cell #" << n-parnum << ":";
          auto k=n;
          while (cellListBak[k]!=-1)
            cout << (k=cellListBak[k]) << " ";
          cout << "\n";
          k=n;
          cout << "B cell #" << n-parnum << ":";
          while (cellList[k]!=-1)
            cout << (k=cellList[k]) << " ";
          cout << "\n\n";
        }
      restore();
    }
  linked_cell_lists_2d(vector<particles>* P, vector<rectangles> *R, pvector<ntype,2>& boxsize)
    {
      init(P,R,boxsize);    
    }
  linked_cell_lists_2d()
    {
      L[0]=L[1]=-1.0;
    }
  void init(vector<particles>* P, pvector<ntype,2>& boxsize, ntype rc,int np=0)
    {
      if (P!=nullptr)
        part = P;
      rcut = rc;
      L = boxsize;
      cellsx = L[0] / rcut;
      cellsy = L[1] / rcut;
      /*
       * cout << "qui np=" << np << "rcut=" << rcut << "\n";
       * cout << "cells=" << cellsx << " " << cellsy << " " << cellsz << "\n";
       * L.show("box size");
       */
      too_few_cells();
      L2 = ntype(0.5)*boxsize;
      parnum=np;//part->size();
      resize_vectors(np);
    }

  void resize_vectors(int N)
    {
      int k;
      parnum=N;
      int oldN=inCell[0].size();
      //cout << "oldN=" << oldN << " new=" << parnum <<  "part=" << (*part).size() << "\n";
      listneigh.resize(parnum);
      if (oldN < parnum)
        {
          for (k=0; k < 2; k++)
            {
              inCell[k].resize(parnum);
              inCellBak[k].resize(parnum);
            }
        }
      int dimN = parnum+cellsx*cellsy;
      if (cellList.size() < dimN)   
        {
          cellList.resize(parnum+cellsx*cellsy);
          cellListBak.resize(parnum+cellsx*cellsy);
        }
    }

  void set_ini_numcells(void)
    {
      cxini = cellsx;
      cyini = cellsy;
    }

  // for NPT simulations
  void update_numcells(pvector<ntype,3> Lb)
    {
      L = Lb;
      L2 = ntype(0.5)*L;
      cellsx = L[0] / rcut;
      cellsy = L[1] / rcut;
      too_few_cells(); 
      if (cellsx*cellsy > cxini*cyini)
        {
          //cout << "qui?!?\n";
          //cout << "cells=" << cellsx << " " << cellsy << " " << cellsz << "\n";
          //cout << "cini="   << cxini << " " << cyini  << " " << czini  << "\n";
          cellList.resize(parnum+cellsx*cellsy);
          cellListBak.resize(parnum+cellsx*cellsy);
        } 
    }
  void store(void)
    {
      parnumBak=parnum;
      Lbak=L;
      L2bak=L2;
      cellsxBak = cellsx;
      cellsyBak = cellsy;
      for (auto k=0; k < parnum + cellsx*cellsy; k++)
        {
          cellListBak[k] = cellList[k];
        }
      for (auto i=0; i < parnum; i++)
        {
          inCellBak[0][i] = inCell[0][i];
          inCellBak[1][i] = inCell[1][i];
        }
    }
  void restore(void)
    {
      L2=L2bak;
      L=Lbak;
      cellsx = cellsxBak;
      cellsy = cellsyBak; 
      for (auto k=0; k < parnum + cellsx*cellsy; k++)
        {
          cellList[k] = cellListBak[k];
        }
      for (auto i=0; i < parnum; i++)
        {
          inCell[0][i] = inCellBak[0][i];
          inCell[1][i] = inCellBak[1][i];
        }
    }
  void findneigh(int ip, int nb=-1, bool debug=false)
    {
      int kk, iX, jX, iY, jY, n, na;
      int cellRange[4];
      pvector<ntype,2> rij;
      pvector<ntype,2> shift;
      na=ip;
      listneigh.clear();
      int cc=0;
      for (kk = 0;  kk < 2; kk++)
        {
          cellRange[2*kk]   = - 1;
          cellRange[2*kk+1] =   1;
        }
      for (iY = cellRange[2]; iY <= cellRange[3]; iY++) 
        {
          jY = inCell[1][na] + iY;    
          shift[1] = 0.0;
          if (jY == -1) 
            {
              jY = cellsy - 1;    
              shift[1] = -L[1];
            } 
          else if (jY == cellsy) 
            {
              jY = 0;    
              shift[1] = L[1];
            }
          for (iX = cellRange[0]; iX <= cellRange[1]; iX++) 
            {
              jX = inCell[0][na] + iX;    
              shift[0] = 0.0;
              if (jX == -1) 
                {
                  jX = cellsx - 1;    
                  shift[0] = - L[0];
                } 
              else if (jX == cellsx) 
                {
                  jX = 0;   
                  shift[0] = L[0];
                }

              n = jY * cellsx + jX + parnum;
              for (n = cellList[n]; n > -1; n = cellList[n]) 
                {
                  if (n != na && n != nb && (nb >= -1 || n < na)) 
                    {
                      // N.B. se sono meno di tre celle
                      // settare lo shift non funziona
                      // poiché la particella n verrebbe
                      // selezionata due volte e in un caso
                      // lo shift è sbagliato!
                      listneigh.push_back(n);
                      //listneigh[cc++] = n;
                      (*part)[n].shift = shift;
                      // minimum image convention
                      cc++;
                    }
                } 
            }
        }
      numneigh=cc;
    }
  void adjust_cell_nums(int& cx, int &cy)
    {
      /* if the particle is on the edge of the box, i.e.
       *
       * (pos[0] + 0.5*L[0]) / L[0] = 1
       * or 
       * (pos[1] + 0.5*L[1]) / L[1] = 1
       * or
       * (pos[2] + 0.5*L[2]) / L[2] = 1
       * it has to belong a valid cell, hence we adjust 
       * cx, cy and cz here.
       * (we are assuming that this can happen only for
       * roundoff errors) */
      if (cx == cellsx )
        cx = cellsx-1;
      if (cy == cellsy )
        cy = cellsy-1;
      if (cx == -1 )
        cx = 0;
      if (cy == -1 )
        cy = 0;
    }
  // after a particle moves...
  void update_LL(int n, bool apply_pbc=false)
    {
      pvector<ntype,3> pos;
      int cox, coy, cx, cy;
      cox=inCell[0][n];
      coy=inCell[1][n];
      pos = (*part)[n].r;
      if (apply_pbc)
        apply_pbc_to_par(pos);
      cx =  (pos[0] + 0.5*L[0]) * cellsx / L[0];
      cy =  (pos[1] + 0.5*L[1]) * cellsy / L[1];
      adjust_cell_nums(cx, cy);
      //L.show("L");
      //pos.show("pos");
      if (cx!=cox || cy!=coy)
        {
          remove_from_current_cell(n);
          inCell[0][n] = cx;
          inCell[1][n] = cy;
          insert_in_new_cell(n);
        }
    }
  // rebuild all linked cell lists
  void apply_pbc_to_par(pvector<ntype,2>& pos)
    {
      pvector<ntype,2> delr;
      delr << L[0]*rint(pos[0]/L[0]), L[1]*rint(pos[1]/L[1]);
      pos -= delr;
    }
  void build(bool apply_pbc=false)
    {
      int j, n;
      pvector<ntype,2> pos;
      for (j = 0; j < cellsx*cellsy + parnum; j++)
        cellList[j] = -1;
      /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */
      
      for (n = 0; n < parnum; n++)
        {
          pos = (*part)[n].r;
          if (apply_pbc)
            apply_pbc_to_par(pos);
          inCell[0][n] =  (pos[0] + L2[0]) * cellsx / L[0];
          inCell[1][n] =  (pos[1] + L2[1]) * cellsy / L[1];
          adjust_cell_nums(inCell[0][n],inCell[1][n],inCell[2][n]);
     
          j = inCell[1][n]*cellsx + 
            inCell[0][n] + parnum;
          cellList[n] = cellList[j];
          cellList[j] = n;
        }
    }
};
#endif

/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National 
   Laboratories ( http://www.mantevo.org ). The primary 
   authors of miniMD are Steve Plimpton (sjplimp@sandia.gov) , Paul Crozier 
   (pscrozi@sandia.gov) and Christian Trott (crtrott@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you 
   can redistribute it and/or modify it under the terms of the GNU Lesser 
   General Public License as published by the Free Software Foundation; 
   either version 3 of the License, or (at your option) any later 
   version.
  
   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
   Lesser General Public License for more details.
    
   You should have received a copy of the GNU Lesser General Public 
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA.  See also: http://www.gnu.org/licenses/lgpl.txt .

   For questions, contact Paul S. Crozier (pscrozi@sandia.gov) or
   Christian Trott (crtrott@sandia.gov). 

   Please read the accompanying README and LICENSE files.
---------------------------------------------------------------------- */

#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include "atom.h"
#include "mcl_wrapper.h"
#include "mcl_data.h"
#include "precision.h"

class Neighbor {
 public:
  int every;                       // re-neighbor every this often
  int nbinx,nbiny,nbinz;           // # of global bins
  MMD_float cutneigh;                 // neighbor cutoff
  MMD_float cutneighsq;               // neighbor cutoff squared
  int ncalls;                      // # of times build has been called
  int max_totalneigh;              // largest # of neighbors ever stored

  int *numneigh;                   // # of neighbors for each atom
  cMCLData<int, xx>* d_numneigh;
  int *neighbors;                  // array of neighbors of each atom
  cMCLData<int, xx>* d_neighbors;
  int maxneighs;				   // max number of neighbors per atom
  int *ilist;                       // ptr to next atom in each bin
  cMCLData<int, xx>* d_ilist;

  MCLWrapper* mcl;
  int ghost_newton;

  Neighbor();
  ~Neighbor();
  int setup(Atom &);                      // setup bins based on box and cutoff
  void resize_buffers(Atom &);
  mcl_handle* binatoms(Atom &);           // bin all atoms
  mcl_handle* resize_and_bin(Atom &);
  mcl_handle* build(Atom &);              // create neighbor list
  mcl_handle* reneigh(Atom &);            // create neighbor list
  cMCLData<int, xx>* d_flag;

  int halfneigh;
  
 private:
  MMD_float xprd,yprd,zprd;           // box size

  int nmax;                        // max size of atom arrays in neighbor
  int *bincount;                    // ptr to 1st atom in each bin
  cMCLData<int, xx>* d_bincount;
  int *bins;                       // ptr to next atom in each bin
  cMCLData<int, xx>* d_bins;
  int *ibins;                       // ptr to next atom in each bin
  cMCLData<int, xx>* d_ibins;
  int atoms_per_bin;

  int nstencil;                    // # of bins in stencil
  int *stencil;                    // stencil list of bin offsets
  cMCLData<int, xx>* d_stencil;


  int mbins;                       // binning parameters
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;
  MMD_float binsizex,binsizey,binsizez;
  MMD_float bininvx,bininvy,bininvz;

  
  MMD_float bindist(int, int, int);   // distance between binx
  int coord2bin(MMD_float, MMD_float, MMD_float);   // mapping atom coord to a bin
};

#endif

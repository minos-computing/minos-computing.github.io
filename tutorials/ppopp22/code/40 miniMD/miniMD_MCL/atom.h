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

#ifndef ATOM_H
#define ATOM_H

#include "mcl_wrapper.h"
#include "mcl_data.h"
#include "precision.h"

struct Box {
  MMD_float xprd,yprd,zprd;
  MMD_float xlo,xhi;
  MMD_float ylo,yhi;
  MMD_float zlo,zhi;
};

class Atom {
 public:
  int natoms;
  int nlocal,nghost;
  int nmax;
  int use_tex;
  MMD_float3 *x;
  MMD_float3 *v;
  MMD_float3 *f;
  MMD_float3 *vold;
  MMD_float mass;
  cMCLData<MMD_float3, xx>* d_x;
  cMCLData<MMD_float3, xx>* d_v;
  cMCLData<MMD_float3, xx>* d_f;
  cMCLData<MMD_float3, xx>* d_vold;
  MCLWrapper* mcl;
  int threads_per_atom;

  int comm_size,reverse_size,border_size;

  struct Box box;

  Atom();
  ~Atom();
  void addatom(MMD_float, MMD_float, MMD_float, MMD_float, MMD_float, MMD_float);
  void pbc();
  void growarray();

  void copy(int, int);

  void pack_comm(int, int *, MMD_float *, int *);
  void unpack_comm(int, int, MMD_float *);
  void pack_reverse(int, int, MMD_float *);
  void unpack_reverse(int, int *, MMD_float *);

  int pack_border(int, MMD_float *, int *);
  int unpack_border(int, MMD_float *);
  int pack_exchange(int, MMD_float *);
  int unpack_exchange(int, MMD_float *);
  int skip_exchange(MMD_float *);

  void cpu_comm_self(int* sendlist, int offset, MMD_float3 pbc, int first, int n);
  
  MMD_float **realloc_2d_MMD_float_array(MMD_float **, int, int, int);
  MMD_float **create_2d_MMD_float_array(int, int);
  void destroy_2d_MMD_float_array(MMD_float **);
};

#endif

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

#include "stdio.h"
#include "math.h"
#include "force.h"

Force::Force() {}
Force::~Force() {}

void Force::setup()
{
  cutforcesq = cutforce*cutforce;
}

mcl_handle* Force::compute(Atom &atom, Neighbor &neighbor, int nwait, mcl_handle** waitlist)
{
  mcl_handle* hdl;
	if(atom.threads_per_atom<0)
	    hdl = mcl->LaunchKernel("force_kernel.h", "force_compute_loop",-(atom.nlocal-atom.threads_per_atom-1)/atom.threads_per_atom, nwait, waitlist, 8,
	    		atom.d_x->devData(),atom.d_x->devSize(), atom.d_x->mclFlags(),
	    		atom.d_f->devData(),atom.d_f->devSize(), atom.d_f->mclFlags(),
	    		neighbor.d_numneigh->devData(),neighbor.d_numneigh->devSize(), neighbor.d_numneigh->mclFlags(),
	    		neighbor.d_neighbors->devData(),neighbor.d_neighbors->devSize(), neighbor.d_neighbors->mclFlags(),
	    		&neighbor.maxneighs,sizeof(neighbor.maxneighs), MCL_ARG_SCALAR,
          &atom.nlocal,sizeof(atom.nlocal), MCL_ARG_SCALAR,
	    		&cutforcesq,sizeof(cutforcesq), MCL_ARG_SCALAR,
          &atom.threads_per_atom,sizeof(atom.threads_per_atom), MCL_ARG_SCALAR);
	else if(atom.threads_per_atom>1)
	    hdl = mcl->LaunchKernel("force_kernel.h", "force_compute_split",atom.nlocal*atom.threads_per_atom, nwait, waitlist, 9,
	    		atom.d_x->devData(),atom.d_x->devSize(), atom.d_x->mclFlags(),
	    		atom.d_f->devData(),atom.d_f->devSize(), atom.d_f->mclFlags(),
	    		neighbor.d_numneigh->devData(),neighbor.d_numneigh->devSize(), neighbor.d_numneigh->mclFlags(),
	    		neighbor.d_neighbors->devData(),neighbor.d_neighbors->devSize(), neighbor.d_neighbors->mclFlags(),
	    		&neighbor.maxneighs,sizeof(neighbor.maxneighs), MCL_ARG_SCALAR,
          &atom.nlocal,sizeof(atom.nlocal), MCL_ARG_SCALAR,
	    		&cutforcesq,sizeof(cutforcesq), MCL_ARG_SCALAR,
          &atom.threads_per_atom,sizeof(atom.threads_per_atom), MCL_ARG_SCALAR,
	    		NULL,sizeof(MMD_float3)*mcl->blockdim, MCL_ARG_LOCAL);
	else if(atom.use_tex)
      //Unsupported image type
      throw "Use TEX unsupported.";
  else {
      // fprintf(stderr, "Launching handle for force compute, nlocal: %d\n", atom.nlocal);
      hdl = mcl->LaunchKernel("force_kernel.h", "force_compute",atom.nlocal, nwait, waitlist, 7,
              atom.d_x->devData(),atom.d_x->devSize(), atom.d_x->mclFlags(),
              atom.d_f->devData(),atom.d_f->devSize(), atom.d_f->mclFlags(),
              neighbor.d_numneigh->devData(),neighbor.d_numneigh->devSize(), neighbor.d_numneigh->mclFlags(),
              neighbor.d_neighbors->devData(),neighbor.d_neighbors->devSize(), neighbor.d_neighbors->mclFlags(),
              &neighbor.maxneighs,sizeof(neighbor.maxneighs), MCL_ARG_SCALAR,
              &atom.nlocal,sizeof(atom.nlocal), MCL_ARG_SCALAR,
              &cutforcesq,sizeof(cutforcesq), MCL_ARG_SCALAR);
  }

  return hdl;
}

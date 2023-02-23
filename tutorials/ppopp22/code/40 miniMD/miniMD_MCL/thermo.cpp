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
#include "stdlib.h"
#include "integrate.h"
#include "thermo.h"

Thermo::Thermo() {}
Thermo::~Thermo() {}

void Thermo::setup(MCLWrapper* w, MMD_float rho_in, Integrate &integrate, Atom &atom,int units,int nparts)
{
  mcl = w;
  rho = rho_in;
  ntimes = integrate.ntimes;
  partitions = nparts;

  int maxstat;
  if (nstat == 0) maxstat = 2;
  else maxstat = ntimes/nstat + 1;
  steparr = (int *) malloc(maxstat*sizeof(int));
  tmparr = (MMD_float *) malloc(maxstat*sizeof(MMD_float));
  engarr = (MMD_float *) malloc(maxstat*sizeof(MMD_float));
  prsarr = (MMD_float *) malloc(maxstat*sizeof(MMD_float));

  if(units == LJ) {
    mvv2e = 1.0;
    dof_boltz = (atom.natoms * 3 - 3);
    t_scale = mvv2e / dof_boltz;
    p_scale = 1.0 / 3 / atom.box.xprd / atom.box.yprd / atom.box.zprd;
    e_scale = 0.5;
  } else if(units == METAL) {
    mvv2e = 1.036427e-04;
    dof_boltz = (atom.natoms * 3 - 3) * 8.617343e-05;
    t_scale = mvv2e / dof_boltz;
    p_scale = 1.602176e+06 / 3 / atom.box.xprd / atom.box.yprd / atom.box.zprd;
    e_scale = 524287.985533;//16.0;
    integrate.dtforce /= mvv2e;

  }

}

void Thermo::compute(int iflag, Atom atom[], Neighbor neighbor[], Force &force, Timer &timer, Comm &comm)
{
  if (iflag > 0 && iflag % nstat) return;
  if (iflag == -1 && nstat > 0 && ntimes % nstat == 0) return;

  t_act=0;
  e_act=0;
  p_act=0;

  MMD_float2** sums = new MMD_float2*[partitions];
  mcl_handle** hdls = new mcl_handle*[partitions];
  MMD_float2 ev = {0.0f, 0.0f};
  for(int i = 0; i < partitions; i++){
    int nblocks = (atom[i].nlocal + mcl->blockdim - 1)/mcl->blockdim;
    sums[i] = new MMD_float2[nblocks];
    hdls[i] = mcl->LaunchKernel("thermo_kernel.h", "energy_virial", atom[i].nlocal, 0, NULL, 8,
      atom[i].d_x->devData(), atom[i].d_x->devSize(), atom[i].d_x->mclFlags(),
      neighbor[i].d_numneigh->devData(), neighbor[i].d_numneigh->devSize(), neighbor[i].d_numneigh->mclFlags(),
      neighbor[i].d_neighbors->devData(), neighbor[i].d_neighbors->devSize(), neighbor[i].d_neighbors->mclFlags(),
      sums[i], nblocks * sizeof(MMD_float2), MCL_ARG_BUFFER | MCL_ARG_OUTPUT,
      NULL, mcl->blockdim * sizeof(MMD_float2), MCL_ARG_BUFFER | MCL_ARG_LOCAL,
      &force.cutforcesq, sizeof(force.cutforcesq), MCL_ARG_SCALAR,
      &neighbor[i].maxneighs, sizeof(neighbor[i].maxneighs), MCL_ARG_SCALAR,
      &atom[i].nlocal, sizeof(atom->nlocal), MCL_ARG_SCALAR
    );
  }
  mcl_wait_all();
  for(int i = 0; i < partitions; i++){
    mcl_hdl_free(hdls[i]);
    int nblocks = (atom[i].nlocal + mcl->blockdim - 1)/mcl->blockdim;
    for(int j = 0; j < nblocks; j++){
      ev.x += sums[i][j].x;
      ev.y += sums[i][j].y;
    }
    delete[] sums[i];
  }
  ev.x  = ev.x/atom[0].natoms;
  delete[] sums;

  MMD_float** temp_sums = new MMD_float*[partitions];
  MMD_float t = 0.0f;
  int nblocks = 64;
  int threads = nblocks * mcl->blockdim;
  for(int i = 0; i < partitions; i++){
    temp_sums[i] = new MMD_float[nblocks];
    hdls[i] = mcl->LaunchKernel("thermo_kernel.h", "temperature", threads, 0, NULL, 4, 
      atom[i].d_v->devData(), atom[i].d_v->devSize(), MCL_ARG_BUFFER | MCL_ARG_INPUT | MCL_ARG_RDONLY,
      temp_sums[i], nblocks * sizeof(MMD_float), MCL_ARG_BUFFER | MCL_ARG_OUTPUT,
      NULL, mcl->blockdim * sizeof(MMD_float3), MCL_ARG_BUFFER | MCL_ARG_LOCAL,
      &atom[i].nlocal, sizeof(atom[i].nlocal), MCL_ARG_SCALAR);
  }
  mcl_wait_all();

  t = 0.0f;
  for(int i = 0; i < partitions; i++){
    mcl_hdl_free(hdls[i]);
    for(int j = 0; j < nblocks; j++){
      //fprintf(stderr, "%f ", temp_sums[i][j]);
      t += temp_sums[i][j];
    }
    //fprintf(stderr, "\nCumulative t: %f\n", t);
    delete[] temp_sums[i];
  }
  delete[] temp_sums;

  t *= t_scale;

  //t = temperature(atom, partitions);

  MMD_float eng = 0.5*ev.x;
  MMD_float p = (t * dof_boltz + 0.5*ev.y) * p_scale;

  int istep = iflag;
  if (iflag == -1) istep = ntimes;

  if (iflag == 0) mstat = 0;

  

  steparr[mstat] = istep;
  tmparr[mstat] = t;
  engarr[mstat] = eng;
  prsarr[mstat] = p;

  mstat++;

  double oldtime=timer.array[TIME_TOTAL];
  timer.stop(TIME_TOTAL);
  fprintf(stdout, "%d ", istep);
  fprintf(stdout, "%e ", t);
  fprintf(stdout, "%e ", eng);
  fprintf(stdout, "%e ", p);
  fprintf(stdout, "%6.3lf\n", istep==0?0.0:timer.array[TIME_TOTAL]);

  // fprintf(stdout,"%i %e %e %e %6.3lf\n",istep,t,eng,p,istep==0?0.0:timer.array[TIME_TOTAL]);

  timer.array[TIME_TOTAL]=oldtime;
  exit(1);
}

MMD_float Thermo::temperature(Atom atom[], int nparts)
{
  int i, j;
  MMD_float vx,vy,vz;

  MMD_float t = 0.0;
  for(j = 0; j < nparts; j++){
    for (i = 0; i < atom[j].nlocal; i++) {
      vx = atom[j].v[i].x;
      vy = atom[j].v[i].y;
      vz = atom[j].v[i].z;
      t += (vx*vx) + (vy*vy) + (vz*vz);
    }
  }

  return t * t_scale;
}

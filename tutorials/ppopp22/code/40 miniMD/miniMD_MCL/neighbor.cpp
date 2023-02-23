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

#include "neighbor.h"

#define FACTOR 0.999
#define SMALL 1.0e-6

Neighbor::Neighbor()
{
  ncalls = 0;
  max_totalneigh = 0;
  numneigh = NULL;
  neighbors = NULL;
  d_numneigh = NULL;
  d_neighbors = NULL;
  maxneighs = 100;
  d_ilist = NULL;
  ilist = NULL;
  nmax = 0;
  bincount = NULL;
  d_bincount = NULL;
  bins = NULL;
  d_bins = NULL;
  ibins = NULL;
  d_ibins = NULL;
  atoms_per_bin = 8;
  stencil = NULL;
  d_stencil = NULL;
  d_flag = NULL;
}

Neighbor::~Neighbor()
{
  delete d_neighbors;
  delete d_numneigh;
  delete d_bincount;
  delete d_bins;
  delete d_ibins;
  delete d_stencil;
  delete d_ilist;
//  if (bincount) free(bincount);
//  if (bins) free(bins);
}

/* binned neighbor list construction with full Newton's 3rd law
   every pair stored exactly once by some processor
   each owned atom i checks its own bin and other bins in Newton stencil */

void Neighbor::resize_buffers(Atom &atom)
{
  int i,j,k,m,n,ibin,jbin,nlocal,nall,npnt;
  MMD_float xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;
  MMD_float3 *x;

  ncalls++;
  nlocal = atom.nlocal;
  nall = atom.nlocal + atom.nghost;

  /* extend atom arrays if necessary */

  if (nall > nmax) {
    if(nmax){
      mcl_unregister_buffer(d_numneigh->devData());
      mcl_unregister_buffer(d_neighbors->devData());
      mcl_unregister_buffer(d_ibins->devData());
      delete d_neighbors;
      delete d_numneigh;
      delete d_ibins;
    }
    
    nmax = nall;
    //printf("Creating buffer for size: %d\n", nmax);
    d_numneigh = new cMCLData<int,xx>(mcl, MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_DYNAMIC, nmax);
    d_neighbors = new cMCLData<int,xx>(mcl, MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_DYNAMIC, nmax*maxneighs);
    d_ibins = new cMCLData<int,xx>(mcl, MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_DYNAMIC, nmax);
    numneigh = d_numneigh->hostData();
    neighbors = d_neighbors->hostData();
    ibins = d_ibins->hostData();
  }
}

mcl_handle* Neighbor::build(Atom &atom) {
  /* loop over each atom, storing neighbors */

  MMD_float3 *x = atom.x;
  d_flag->hostData()[0]=0;
  d_flag->upload();

  mcl_handle* hdl = mcl->LaunchKernel("neighbor_kernel.h", "neighbor_build",atom.nlocal, 0, NULL, 13,
    atom.d_x->devData(),atom.d_x->devSize(), atom.d_x->mclFlags(),
    d_numneigh->devData(),d_numneigh->devSize(), d_numneigh->mclFlags(),
    d_neighbors->devData(),d_neighbors->devSize(), d_neighbors->mclFlags(),
    d_bincount->devData(),d_bincount->devSize(), d_bincount->mclFlags(),
    d_bins->devData(),d_bins->devSize(), d_bins->mclFlags(),
    d_ibins->devData(),d_ibins->devSize(), d_ibins->mclFlags(),
    d_flag->devData(),d_flag->devSize(), d_flag->mclFlags(),
    d_stencil->devData(),d_stencil->devSize(), d_stencil->mclFlags(),
    &nstencil,sizeof(nstencil), MCL_ARG_SCALAR, 
    &cutneighsq,sizeof(cutneighsq), MCL_ARG_SCALAR,
    &atoms_per_bin,sizeof(atoms_per_bin), MCL_ARG_SCALAR, 
    &maxneighs,sizeof(maxneighs), MCL_ARG_SCALAR,
    &atom.nlocal, sizeof(atom.nlocal), MCL_ARG_SCALAR
    );

  return hdl;
}


mcl_handle* Neighbor::reneigh(Atom &atom) {
  d_flag->download();
  if(d_flag->hostData()[0])
  {
    mcl_unregister_buffer(d_neighbors->devData());
    delete d_neighbors;
    maxneighs *= 1.5;
    d_neighbors = new cMCLData<int,xx>(mcl, MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_DYNAMIC, nmax*maxneighs);
    neighbors = d_neighbors->hostData();
    //fprintf(stderr, "Creating new handle, returning...\n");
    return build(atom);
  }

  return NULL;
  
}
      
/* bin owned and ghost atoms */

mcl_handle* Neighbor::binatoms(Atom &atom)
{
  int ibin,nlocal,nall;
  MMD_float3 *x;

  nlocal = atom.nlocal;
  nall = atom.nlocal + atom.nghost;
  x = atom.x;

  xprd = atom.box.xprd;
  yprd = atom.box.yprd;
  zprd = atom.box.zprd;

#if MDPREC == 2
    MMD_floatK3 bininv;
#else
    MMD_float3 bininv;
#endif
  bininv.x=bininvx;bininv.y=bininvy;bininv.z=bininvz;
  MMD_float3 prd;
  prd.x=xprd;prd.y=yprd;prd.z=zprd;
  cl_int3 mbinlo;
  mbinlo.x=mbinxlo;mbinlo.y=mbinylo;mbinlo.z=mbinzlo;
  cl_int3 nbin;
  nbin.x=nbinx;nbin.y=nbiny;nbin.z=nbinz;
  cl_int3 mbin;
  mbin.x=mbinx;mbin.y=mbiny;mbin.z=mbinz;

  for(int i=0; i< mbins; i++) bincount[i]=0;
  d_flag->hostData()[0]=0;
  d_flag->upload();
  d_bincount->upload();

  //fprintf(stderr, "Launching kernel...");
  mcl_handle* hdl = mcl->LaunchKernel("neighbor_kernel.h", "neighbor_bin",nall, 0, NULL, 12,
    atom.d_x->devData(),atom.d_x->devSize(), atom.d_x->mclFlags() | MCL_ARG_REWRITE,
    d_bincount->devData(),d_bincount->devSize(), d_bincount->mclFlags() | MCL_ARG_REWRITE,
    d_bins->devData(),d_bins->devSize(), d_bins->mclFlags(),
    d_ibins->devData(),d_ibins->devSize(), d_ibins->mclFlags(),
    d_flag->devData(),d_flag->devSize(), d_flag->mclFlags(),
    &atoms_per_bin,sizeof(atoms_per_bin), MCL_ARG_SCALAR,
    &nall,sizeof(nall), MCL_ARG_SCALAR,
    &bininv,sizeof(bininv), MCL_ARG_SCALAR,
    &prd,sizeof(prd), MCL_ARG_SCALAR, 
    &mbinlo,sizeof(mbinlo), MCL_ARG_SCALAR,
    &nbin,sizeof(nbin), MCL_ARG_SCALAR, 
    &mbin,sizeof(mbin), MCL_ARG_SCALAR
    );
  //fprintf(stderr, "Kernel Launched.\n");

  return hdl;
}

mcl_handle* Neighbor::resize_and_bin(Atom &atom)
{
  int ibin,nlocal,nall;
  MMD_float3 *x;

  nlocal = atom.nlocal;
  nall = atom.nlocal + atom.nghost;
  x = atom.x;

  xprd = atom.box.xprd;
  yprd = atom.box.yprd;
  zprd = atom.box.zprd;

#if MDPREC == 2
    MMD_floatK3 bininv;
#else
    MMD_float3 bininv;
#endif
  bininv.x=bininvx;bininv.y=bininvy;bininv.z=bininvz;
  MMD_float3 prd;
  prd.x=xprd;prd.y=yprd;prd.z=zprd;
  cl_int3 mbinlo;
  mbinlo.x=mbinxlo;mbinlo.y=mbinylo;mbinlo.z=mbinzlo;
  cl_int3 nbin;
  nbin.x=nbinx;nbin.y=nbiny;nbin.z=nbinz;
  cl_int3 mbin;
  mbin.x=mbinx;mbin.y=mbiny;mbin.z=mbinz;

  d_flag->download();
  if(!d_flag->hostData()[0]) return NULL;
  mcl_unregister_buffer(d_bins->devData());
  delete d_bins;
  atoms_per_bin*=2;
  d_bins = new cMCLData<int,xx>(mcl, MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_DYNAMIC, mbins*atoms_per_bin);
  bins = d_bins->hostData();
	for(int i=0; i< mbins; i++) bincount[i]=0;
  d_flag->hostData()[0]=0;
  d_flag->upload();
  d_bincount->upload();
  atom.d_x->upload();
  mcl_handle* hdl = mcl->LaunchKernel("neighbor_kernel.h", "neighbor_bin",nall, 0, NULL, 12,
    atom.d_x->devData(),atom.d_x->devSize(), atom.d_x->mclFlags(),
    d_bincount->devData(),d_bincount->devSize(), d_bincount->mclFlags() | MCL_ARG_REWRITE,
    d_bins->devData(),d_bins->devSize(), d_bins->mclFlags(),
    d_ibins->devData(),d_ibins->devSize(), d_ibins->mclFlags(),
    d_flag->devData(),d_flag->devSize(), d_flag->mclFlags(),
    &atoms_per_bin, sizeof(atoms_per_bin), MCL_ARG_SCALAR,
     &nall,sizeof(nall), MCL_ARG_SCALAR,
    &bininv,sizeof(bininv), MCL_ARG_SCALAR,
     &prd,sizeof(prd), MCL_ARG_SCALAR, 
    &mbinlo,sizeof(mbinlo), MCL_ARG_SCALAR,
     &nbin,sizeof(nbin), MCL_ARG_SCALAR, 
    &mbin,sizeof(mbin), MCL_ARG_SCALAR
  );

  return hdl;
}



/* convert xyz atom coords into local bin #
   take special care to insure ghost atoms with
   coord >= prd or coord < 0.0 are put in correct bins */

int Neighbor::coord2bin(MMD_float x, MMD_float y, MMD_float z)
{
  int ix,iy,iz;
  if (x >= xprd)
    ix = (int) ((x-xprd)*bininvx) + nbinx - mbinxlo;
  else if (x >= 0.0)
    ix = (int) (x*bininvx) - mbinxlo;
  else
    ix = (int) (x*bininvx) - mbinxlo - 1;
  
  if (y >= yprd)
    iy = (int) ((y-yprd)*bininvy) + nbiny - mbinylo;
  else if (y >= 0.0)
    iy = (int) (y*bininvy) - mbinylo;
  else
    iy = (int) (y*bininvy) - mbinylo - 1;
  
  if (z >= zprd)
    iz = (int) ((z-zprd)*bininvz) + nbinz - mbinzlo;
  else if (z >= 0.0)
    iz = (int) (z*bininvz) - mbinzlo;
  else
    iz = (int) (z*bininvz) - mbinzlo - 1;
 // printf("%i",iz*mbiny*mbinx + iy*mbinx + ix + 1);
  return (iz*mbiny*mbinx + iy*mbinx + ix + 1);
}


/*
setup neighbor binning parameters
bin numbering is global: 0 = 0.0 to binsize
                         1 = binsize to 2*binsize
                         nbin-1 = prd-binsize to binsize
                         nbin = prd to prd+binsize
                         -1 = -binsize to 0.0
coord = lowest and highest values of ghost atom coords I will have
        add in "small" for round-off safety
mbinlo = lowest global bin any of my ghost atoms could fall into
mbinhi = highest global bin any of my ghost atoms could fall into
mbin = number of bins I need in a dimension
stencil() = bin offsets in 1-d sense for stencil of surrounding bins
*/

int Neighbor::setup(Atom &atom)
{
  int i,j,k,nmax;
  MMD_float coord;
  int mbinxhi,mbinyhi,mbinzhi;
  int nextx,nexty,nextz;
 
  cutneighsq = cutneigh*cutneigh;
 
  xprd = atom.box.xprd;
  yprd = atom.box.yprd;
  zprd = atom.box.zprd;

  /*
c bins must evenly divide into box size, 
c   becoming larger than cutneigh if necessary
c binsize = 1/2 of cutoff is near optimal

  if (flag == 0) {
    nbinx = 2.0 * xprd / cutneigh;
    nbiny = 2.0 * yprd / cutneigh;
    nbinz = 2.0 * zprd / cutneigh;
    if (nbinx == 0) nbinx = 1;
    if (nbiny == 0) nbiny = 1;
    if (nbinz == 0) nbinz = 1;
  }
  */

  binsizex = xprd/nbinx;
  binsizey = yprd/nbiny;
  binsizez = zprd/nbinz;
  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;

  coord = atom.box.xlo - cutneigh - SMALL*xprd;
  mbinxlo = static_cast<int>(coord*bininvx);
  if (coord < 0.0) mbinxlo = mbinxlo - 1;
  coord = atom.box.xhi + cutneigh + SMALL*xprd;
  mbinxhi = static_cast<int>(coord*bininvx);

  coord = atom.box.ylo - cutneigh - SMALL*yprd;
  mbinylo = static_cast<int>(coord*bininvy);
  if (coord < 0.0) mbinylo = mbinylo - 1;
  coord = atom.box.yhi + cutneigh + SMALL*yprd;
  mbinyhi = static_cast<int>(coord*bininvy);

  coord = atom.box.zlo - cutneigh - SMALL*zprd;
  mbinzlo = static_cast<int>(coord*bininvz);
  if (coord < 0.0) mbinzlo = mbinzlo - 1;
  coord = atom.box.zhi + cutneigh + SMALL*zprd;
  mbinzhi = static_cast<int>(coord*bininvz);

/* extend bins by 1 in each direction to insure stencil coverage */

  mbinxlo = mbinxlo - 1;
  mbinxhi = mbinxhi + 1;
  mbinx = mbinxhi - mbinxlo + 1;

  mbinylo = mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  mbiny = mbinyhi - mbinylo + 1;

  mbinzlo = mbinzlo - 1;
  mbinzhi = mbinzhi + 1;
  mbinz = mbinzhi - mbinzlo + 1;

  /*
compute bin stencil of all bins whose closest corner to central bin
  is within neighbor cutoff
for partial Newton (newton = 0),
  stencil is all surrounding bins including self
for full Newton (newton = 1),
  stencil is bins to the "upper right" of central bin, does NOT include self
next(xyz) = how far the stencil could possibly extend
factor < 1.0 for special case of LJ benchmark so code will create
  correct-size stencil when there are 3 bins for every 5 lattice spacings
  */

  nextx = static_cast<int>(cutneigh*bininvx);
  if (nextx*binsizex < FACTOR*cutneigh) nextx++;
  nexty = static_cast<int>(cutneigh*bininvy);
  if (nexty*binsizey < FACTOR*cutneigh) nexty++;
  nextz = static_cast<int>(cutneigh*bininvz);
  if (nextz*binsizez < FACTOR*cutneigh) nextz++;

  nmax = (2*nextz+1) * (2*nexty+1) * (2*nextx+1);
  delete d_stencil;
  d_stencil = new cMCLData<int,xx>(mcl, MCL_ARG_INPUT | MCL_ARG_RESIDENT | MCL_ARG_BUFFER | MCL_ARG_RDONLY, nmax);
  stencil = d_stencil->hostData();


  nstencil = 0;
  for (k = -nextz; k <= nextz; k++) {
    for (j = -nexty; j <= nexty; j++) {
      for (i = -nextx; i <= nextx; i++) {
	    if (bindist(i,j,k) < cutneighsq) {
	      stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
	    }
      }
    }
  }
  d_stencil->upload();

  mbins = mbinx*mbiny*mbinz;
  d_bincount = new cMCLData<int,xx>(mcl, MCL_ARG_RESIDENT | MCL_ARG_INPUT | MCL_ARG_DYNAMIC | MCL_ARG_BUFFER, mbins);
  bincount = d_bincount->hostData();
  d_bins = new cMCLData<int,xx>(mcl, MCL_ARG_RESIDENT | MCL_ARG_DYNAMIC | MCL_ARG_BUFFER, mbins*atoms_per_bin);
  bins = d_bins->hostData();
  d_flag = new cMCLData<int,xx>(mcl, MCL_ARG_BUFFER | MCL_ARG_INPUT | MCL_ARG_OUTPUT | MCL_ARG_RESIDENT | MCL_ARG_REWRITE, 1);
  return 0;
}
      
/* compute closest distance between central bin (0,0,0) and bin (i,j,k) */

MMD_float Neighbor::bindist(int i, int j, int k)
{
  MMD_float delx,dely,delz;

  if (i > 0)
    delx = (i-1)*binsizex;
  else if (i == 0)
    delx = 0.0;
  else
    delx = (i+1)*binsizex;

  if (j > 0)
    dely = (j-1)*binsizey;
  else if (j == 0)
    dely = 0.0;
  else
    dely = (j+1)*binsizey;

  if (k > 0)
    delz = (k-1)*binsizez;
  else if (k == 0)
    delz = 0.0;
  else
    delz = (k+1)*binsizez;
  
  return (delx*delx + dely*dely + delz*delz);
}

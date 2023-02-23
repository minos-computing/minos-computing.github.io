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

#include <algorithm>
#include <iterator>
#include "stdio.h"
#include "stdlib.h"
#include "comm.h"

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 100
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))


void Comm::get_my_loc(int my_loc[], int id){
  my_loc[0] = id / (procgrid[1] * procgrid[2]);
  my_loc[1] = (id % (procgrid[1] * procgrid[2])) / procgrid[2];
  my_loc[2] = (id % (procgrid[1] * procgrid[2])) % procgrid[2];
}

int Comm::neighbor(int my_loc[], int dim, int dir){
  int n_loc[3];
  n_loc[0] = my_loc[0];
  n_loc[1] = my_loc[1];
  n_loc[2] = my_loc[2];
  n_loc[dim] += dir;
  if(n_loc[dim] < 0) n_loc[dim] += procgrid[dim];
  if(n_loc[dim] >= procgrid[dim]) n_loc[dim] %= procgrid[dim];
  return (n_loc[0] * procgrid[1] * procgrid[2]) + (n_loc[1] * procgrid[2]) + n_loc[2];
}

Comm::Comm()
{
  maxsend = NULL;

}

Comm::~Comm()
{
}

/* setup spatial-decomposition communication patterns */

int Comm::setup(MMD_float cutneigh, Atom atom[], int nparts)
{
  int i;
  int periods[3];
  MMD_float prd[3];
  int myloc[3];
  double lo,hi;
  int ineed,idim,nbox;
  
  npatitions = nparts;
  
  prd[0] = atom[0].box.xprd;
  prd[1] = atom[0].box.yprd;
  prd[2] = atom[0].box.zprd;

  /* setup 3-d grid of procs */

  MMD_float area[3];
  area[0] = prd[0] * prd[1];
  area[1] = prd[0] * prd[2];
  area[2] = prd[1] * prd[2];

  MMD_float bestsurf = 2.0 * (area[0]+area[1]+area[2]);

  // loop thru all possible factorizations of nprocs
  // surf = surface area of a proc sub-domain
  // for 2d, insure ipz = 1

  int ipx,ipy,ipz,nremain;
  MMD_float surf;

  ipx = 1;
  while (ipx <= nparts) {
    if (nparts % ipx == 0) {
      nremain = nparts/ipx;
      ipy = 1;
      while (ipy <= nremain) {
        if (nremain % ipy == 0) {
          ipz = nremain/ipy;
          surf = area[0]/ipx/ipy + area[1]/ipx/ipz + area[2]/ipy/ipz;
          if (surf < bestsurf) {
            bestsurf = surf;
            procgrid[0] = ipx;
            procgrid[1] = ipy;
            procgrid[2] = ipz;
          }
        }
        ipy++;
      }
    }
    ipx++;
  }
  if (procgrid[0]*procgrid[1]*procgrid[2] != nparts) {
    printf("ERROR: Bad grid of processors\n");
    return 1;
  }

  //fprintf(stderr, "Process Grid: %d x %d x %d\n", procgrid[0], procgrid[1], procgrid[2]);

  /* determine where I am and my neighboring procs in 3d grid of procs */
  /* lo/hi = my local box bounds */

  for(i = 0; i < nparts; i++){
    get_my_loc(myloc, i);
    atom[i].box.xlo = myloc[0] * prd[0] / procgrid[0];
    atom[i].box.xhi = (myloc[0]+1) * prd[0] / procgrid[0];
    atom[i].box.ylo = myloc[1] * prd[1] / procgrid[1];
    atom[i].box.yhi = (myloc[1]+1) * prd[1] / procgrid[1];
    atom[i].box.zlo = myloc[2] * prd[2] / procgrid[2];
    atom[i].box.zhi = (myloc[2]+1) * prd[2] / procgrid[2];
  }
  

  /* need = # of boxes I need atoms from in each dimension */

  need[0] = static_cast<int>(cutneigh * procgrid[0] / prd[0] + 1);
  need[1] = static_cast<int>(cutneigh * procgrid[1] / prd[1] + 1);
  need[2] = static_cast<int>(cutneigh * procgrid[2] / prd[2] + 1);
 
  /* alloc comm memory */

  maxswap = 2 * (need[0]+need[1]+need[2]);
  

  slablo = (MMD_float *) malloc(nparts*maxswap*sizeof(MMD_float));
  slabhi = (MMD_float *) malloc(nparts*maxswap*sizeof(MMD_float));
  pbc_any = (int *) malloc(nparts*maxswap*sizeof(int));
  pbc_flagx = (int *) malloc(nparts*maxswap*sizeof(int));
  pbc_flagy = (int *) malloc(nparts*maxswap*sizeof(int));
  pbc_flagz = (int *) malloc(nparts*maxswap*sizeof(int));
  recvproc = (int *) malloc(nparts*maxswap*sizeof(int));

  sendnum = (int *) malloc(nparts*maxswap*sizeof(int));
  recvnum = (int *) malloc(nparts*maxswap*sizeof(int));
  firstrecv = (int *) malloc(nparts*maxswap*sizeof(int));


  d_sendlist = (cMCLData<int, xy>**)malloc(sizeof(cMCLData<int, xy>*) * npatitions);
  sendlist = (int ***) malloc(npatitions*sizeof(int**));
  for(int j = 0; j < npatitions; j++){
    sendlist[j] = (int **) malloc(maxswap*sizeof(int*));
    for (i = 0; i < maxswap; i++)
      sendlist[j][i] = (int *) malloc(BUFMIN*sizeof(int));

    d_sendlist[j] = new cMCLData<int,xy>(mcl,(int*)sendlist[j], MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_DYNAMIC | MCL_ARG_INPUT,maxswap,BUFMIN, 0);
  }
  maxsendlist = (int *) malloc(nparts*maxswap*sizeof(int));
  for (i = 0; i < nparts*maxswap; i++) maxsendlist[i] = BUFMIN;
  

  /* setup 4 parameters for each exchange: (spart,rpart,slablo,slabhi)
     recvproc(nswap) = proc to recv from at each swap
     slablo/slabhi(nswap) = slab boundaries (in correct dimension) of atoms
                            to send at each swap
     1st part of if statement is sending to the west/south/down
     2nd part of if statement is sending to the east/north/up
     nbox = atoms I send originated in this box */
  
  /* set commflag if atoms are being exchanged across a box boundary
     commflag(idim,nswap) =  0 -> not across a boundary
                          =  1 -> add box-length to position when sending
                          = -1 -> subtract box-length from pos when sending */

  temp_buffers = new cMCLData<MMD_float, xx>**[nparts];
  maxsend = new int[nparts];
  for(i = 0; i < nparts; i++){
    maxsend[i] = BUFMIN;
    get_my_loc(myloc, i);

    nswap = 0;
    for (idim = 0; idim < 3; idim++) {
      for (ineed = 0; ineed < 2*need[idim]; ineed++) {
        pbc_any[(i*maxswap) + nswap] = 0;
        pbc_flagx[(i*maxswap) + nswap] = 0;
        pbc_flagy[(i*maxswap) + nswap] = 0;
        pbc_flagz[(i*maxswap) + nswap] = 0;

        if (ineed % 2 == 0) {
          recvproc[(i*maxswap) + nswap] = neighbor(myloc, idim, 1);

          nbox = myloc[idim] + ineed/2;
          lo = nbox * prd[idim] / procgrid[idim];
          if (idim == 0) hi = atom[i].box.xlo + cutneigh;
          if (idim == 1) hi = atom[i].box.ylo + cutneigh;
          if (idim == 2) hi = atom[i].box.zlo + cutneigh;
          hi = MIN(hi,(nbox+1) * prd[idim] / procgrid[idim]);
          if (myloc[idim] == 0) {
            pbc_any[(i*maxswap) + nswap] = 1;
            if (idim == 0) pbc_flagx[(i*maxswap) + nswap] = 1;
            if (idim == 1) pbc_flagy[(i*maxswap) + nswap] = 1;
            if (idim == 2) pbc_flagz[(i*maxswap) + nswap] = 1;
          }
        } else {
          recvproc[(i*maxswap) + nswap] = neighbor(myloc, idim, -1);
          
          nbox = myloc[idim] - ineed/2;
          hi = (nbox+1) * prd[idim] / procgrid[idim];
          if (idim == 0) lo = atom[i].box.xhi - cutneigh;
          if (idim == 1) lo = atom[i].box.yhi - cutneigh;
          if (idim == 2) lo = atom[i].box.zhi - cutneigh;
          lo = MAX(lo,nbox * prd[idim] / procgrid[idim]);
          if (myloc[idim] == procgrid[idim]-1) {
            pbc_any[(i*maxswap) + nswap] = 1;
            if (idim == 0) pbc_flagx[(i*maxswap) + nswap] = -1;
            if (idim == 1) pbc_flagy[(i*maxswap) + nswap] = -1;
            if (idim == 2) pbc_flagz[(i*maxswap) + nswap] = -1;
          }
        }

        slablo[(i*maxswap)+nswap] = lo;
        slabhi[(i*maxswap)+nswap] = hi;
        nswap++;
      }
    }
    temp_buffers[i] = new cMCLData<MMD_float, xx>*[nswap];
    for(int j = 0; j < nswap; j++)
      temp_buffers[i][j] = new cMCLData<MMD_float, xx>(mcl, MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_DYNAMIC, maxsend[i], 0, 0);
  }



  return 0;
}

/* communication of atom info every timestep */

mcl_handle** Comm::communicate(Atom atom[], int i, mcl_handle** waitlist)
{
  int partition, iswap;
  int pbc_flags[4];
  MMD_float *buf;
  mcl_handle** hdls = new mcl_handle*[npatitions * nswap];
  mcl_handle** hdls_2 = new mcl_handle*[npatitions * nswap];
  uint64_t rewrite = i == 0 ? MCL_ARG_REWRITE : 0;
  //fprintf(stderr, "Starting communicate..."); 
  for(partition = 0; partition < npatitions; partition++){
    for (iswap = 0; iswap < nswap; iswap++) {
      /* pack buffer */
      pbc_flags[0] = pbc_any[(partition*maxswap) + iswap];
      pbc_flags[1] = pbc_flagx[(partition*maxswap) + iswap];
      pbc_flags[2] = pbc_flagy[(partition*maxswap) + iswap];
      pbc_flags[3] = pbc_flagz[(partition*maxswap) + iswap];
      
      MMD_float3 pbc;
      pbc.x=atom[partition].box.xprd*pbc_flagx[(partition*maxswap) + iswap];
      pbc.y=atom[partition].box.yprd*pbc_flagy[(partition*maxswap) + iswap];
      pbc.z=atom[partition].box.zprd*pbc_flagz[(partition*maxswap) + iswap];
      int offset = iswap * maxsendlist[(partition*maxswap)];
      //atom.d_x->download();
      //atom.pack_comm(sendnum[iswap],sesndlist[iswap],buf_send,pbc_flags);

      /* exchange with another proc
        if self, set recv buffer to send buffer */

      if (recvproc[(partition * maxswap) + iswap] != partition) {
        hdls[(partition * maxswap) + iswap] = mcl->LaunchKernel("atom_kernel.h", "atom_pack_comm", sendnum[(partition * maxswap) + iswap], 1, &waitlist[partition], 6,
            atom[partition].d_x->devData(),atom[partition].d_x->devSize(),atom[partition].d_x->mclFlags(),
            temp_buffers[partition][iswap]->devData(),temp_buffers[partition][iswap]->devSize(),temp_buffers[partition][iswap]->mclFlags(),
            d_sendlist[partition]->devData(),d_sendlist[partition]->devSize(),d_sendlist[partition]->mclFlags() | rewrite,
            &offset,sizeof(offset), MCL_ARG_SCALAR,
            &pbc,sizeof(pbc), MCL_ARG_SCALAR, 
            &sendnum[(partition * maxswap) + iswap],sizeof(sendnum[(partition * maxswap) + iswap]), MCL_ARG_SCALAR
        );
      } else {
        //atom[partition].cpu_comm_self(d_sendlist[partition]->devData(), offset, pbc, firstrecv[(partition*maxswap) + iswap], sendnum[(partition*maxswap) + iswap]);
        hdls[(partition * maxswap) + iswap] = mcl->LaunchKernel("atom_kernel.h", "atom_comm_self", sendnum[(partition * maxswap) + iswap], 1, &waitlist[partition], 6,
            atom[partition].d_x->devData(),atom[partition].d_x->devSize(),atom[partition].d_x->mclFlags(),
            d_sendlist[partition]->devData(),d_sendlist[partition]->devSize(),d_sendlist[partition]->mclFlags() | rewrite,
            &offset,sizeof(offset), MCL_ARG_SCALAR,
            &pbc,sizeof(pbc), MCL_ARG_SCALAR, 
            &firstrecv[(partition * maxswap) + iswap],sizeof(firstrecv[(partition * maxswap) + iswap]), MCL_ARG_SCALAR,
            &sendnum[(partition * maxswap) + iswap],sizeof(sendnum[(partition * maxswap) + iswap]), MCL_ARG_SCALAR);
      }
    }
  }

  for (iswap = 0; iswap < nswap; iswap++) {
    for(partition = 0; partition < npatitions; partition++){
      int recv = recvproc[(partition * maxswap) + iswap];
      if (recv != partition) {
        mcl_handle** wait = &hdls[(recv * maxswap) + iswap];
        hdls_2[(partition * maxswap) + iswap] =  mcl->LaunchKernel("atom_kernel.h", "atom_unpack_comm", recvnum[(partition * maxswap) + iswap], 1, wait, 4,
            atom[partition].d_x->devData(),atom[partition].d_x->devSize(),atom[partition].d_x->mclFlags(),
            temp_buffers[recv][iswap]->devData(),temp_buffers[recv][iswap]->devSize(),temp_buffers[recv][iswap]->mclFlags(),
            &firstrecv[(partition * maxswap) + iswap],sizeof(firstrecv[(partition * maxswap) + iswap]),MCL_ARG_SCALAR,
            &recvnum[(partition * maxswap) + iswap],sizeof(recvnum[(partition * maxswap) + iswap]),MCL_ARG_SCALAR
        );
      } else {
        hdls_2[(partition * maxswap) + iswap] =  hdls[(recv * maxswap) + iswap];
        hdls[(recv * maxswap) + iswap] = NULL;
      }
    }
  }
  for (iswap = 0; iswap < nswap; iswap++) {
    for(partition = 0; partition < npatitions; partition++) {
      if(hdls[(partition * maxswap) + iswap])
        to_free.push_back(hdls[(partition * maxswap) + iswap]);
      to_free.push_back(hdls_2[(partition * maxswap) + iswap]);  
    }
  }
  //fprintf(stderr, "done.\n");
  delete[] hdls;

  return hdls_2;
}

/* reverse communication of atom info every timestep */
      
void Comm::reverse_communicate(Atom atom[])
{
  // Unimplemented: Needed for half neighbor
  throw "Not Implemented";
}

/* exchange:
   move atoms to correct proc boxes
   send out atoms that have left my box, receive ones entering my box
   this routine called before every reneighboring
   atoms exchanged with all 6 stencil neighbors
*/

void Comm::exchange(Atom atom[])
{
  int i,j,m,n,idim,nrecv,nlocal;
  int* nsend = new int[npatitions];
  int myloc[3];
  MMD_float lo,hi,value;
  MMD_float3 *x;

  /* enforce PBC */

  for(i = 0; i < npatitions; i++){
    atom[i].pbc();
  }
  
  /* loop over dimensions */

  for (idim = 0; idim < 3; idim++) {

    /* only exchange if more than one proc in this dimension */
    //fprintf(stderr, "dim: %d, size: %d\n", idim, procgrid[idim]);

    if (procgrid[idim] == 1) continue;

    /* fill buffer with atoms leaving my box
       when atom is deleted, fill it in with last atom */

    for(j = 0; j < npatitions; j++) nsend[j] = 0;

    for(j = 0; j < npatitions; j++){
      i = 0;
      MMD_float* buf_send = temp_buffers[j][0]->hostData();

      if (idim == 0) {
        lo = atom[j].box.xlo;
        hi = atom[j].box.xhi;
      } else if (idim == 1) {
        lo = atom[j].box.ylo;
        hi = atom[j].box.yhi;
      } else {
        lo = atom[j].box.zlo;
        hi = atom[j].box.zhi;
      }

      x = atom[j].x;

      nlocal = atom[j].nlocal;

      MMD_float xdim;

      while (i < nlocal) {
        if(idim==0) xdim=x[i].x;
        if(idim==1) xdim=x[i].y;
        if(idim==2) xdim=x[i].z;

        if (xdim < lo || xdim >= hi) {
          if (nsend[j] + 6 >= maxsend[j]) buf_send = growsend(nsend[j], j, 0);
          nsend[j] += atom[j].pack_exchange(i,&buf_send[nsend[j]]);
          atom[j].copy(nlocal-1,i);
          nlocal--;
        } else i++;
      }
      atom[j].nlocal = nlocal;
    }
          
    /* check incoming atoms to see if they are in my box
       if they are, add to my list */

    for(j = 0; j < npatitions; j++) {
      get_my_loc(myloc, j);
      n = atom[j].nlocal;

      if (idim == 0) {
        lo = atom[j].box.xlo;
        hi = atom[j].box.xhi;
      } else if (idim == 1) {
        lo = atom[j].box.ylo;
        hi = atom[j].box.yhi;
      } else {
        lo = atom[j].box.zlo;
        hi = atom[j].box.zhi;
      }

      m = 0;

      nrecv = nsend[neighbor(myloc, idim, -1)];
      MMD_float* buf_recv =  temp_buffers[neighbor(myloc, idim, -1)][0]->hostData();
      int num_zeros = 0;
      while (m < nrecv) {
        //fprintf(stderr, "Recieving atom from: %d\n", m);
        value = buf_recv[m+idim];
        if(value == 0){
          num_zeros++;
        }
        if (value >= lo && value < hi)
          m += atom[j].unpack_exchange(n++,&buf_recv[m]);
        else m += atom[j].skip_exchange(&buf_recv[m]);
      }

      if(procgrid[idim] > 2){
        m = 0;
        nrecv = nsend[neighbor(myloc, idim, 1)];
        buf_recv =  temp_buffers[neighbor(myloc, idim, 1)][0]->hostData();
        while (m < nrecv) {
          value = buf_recv[m+idim];
          if (value >= lo && value < hi)
            m += atom[j].unpack_exchange(n++,&buf_recv[m]);
          else m += atom[j].skip_exchange(&buf_recv[m]);
        }
      }
      atom[j].nlocal = n;
    }
  }

  delete[] nsend;
}

/* borders:
   make lists of nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a communicate (so don't need to explicitly
     call communicate routine on reneighboring timestep)
   this routine is called before every reneighboring
*/

void Comm::borders(Atom* atom)
{
  int j,i,m,n,iswap,idim,ineed,nrecv;
  int* nsend = new int[npatitions];
  int* nfirst = new int[npatitions];
  int* nlast = new int[npatitions];
  MMD_float lo,hi;
  int pbc_flags[4];
  MMD_float3 *x;
  MMD_float *buf;

  /* erase all ghost atoms */

  for(j = 0; j < npatitions; j++) atom[j].nghost = 0;
  for(j = 0; j < npatitions; j++) nfirst[j] = 0;

  /* do swaps over all 3 dimensions */

  iswap = 0;
  

  for (idim = 0; idim < 3; idim++) {
    for(j = 0; j < npatitions; j++) nlast[j] = 0;
    for (ineed = 0; ineed < 2*need[idim]; ineed++) {
      /* find all atoms (own & ghost) within slab boundaries lo/hi
      store atom indices in list for use in future timesteps */

      for(j = 0; j < npatitions; j++){
        lo = slablo[(j*maxswap) + iswap];
        hi = slabhi[(j*maxswap) + iswap];
        pbc_flags[0] = pbc_any[(j*maxswap) + iswap];
        pbc_flags[1] = pbc_flagx[(j*maxswap) + iswap];
        pbc_flags[2] = pbc_flagy[(j*maxswap) + iswap];
        pbc_flags[3] = pbc_flagz[(j*maxswap) + iswap];


        x = atom[j].x;

        if (ineed % 2 == 0) {
          nfirst[j] = nlast[j];
          nlast[j] = atom[j].nlocal + atom[j].nghost;
        }

        nsend[j] = 0;
        m = 0;

        MMD_float xdim;
        MMD_float* buf_send = temp_buffers[j][iswap]->hostData();
        for (i = nfirst[j]; i < nlast[j]; i++) {
          if(idim==0) xdim=x[i].x;
          if(idim==1) xdim=x[i].y;
          if(idim==2) xdim=x[i].z;
          if (xdim >= lo && xdim < hi) {
            if (m + 3 >= maxsend[j]) buf_send = growsend(m, j, iswap);
            m += atom[j].pack_border(i,&buf_send[m],pbc_flags);
            if (nsend[j] >= maxsendlist[(j*maxswap) + iswap]) growlist(iswap,nsend[j],j);
            sendlist[j][iswap][nsend[j]++] = i;
          }
        }
      }

      for(j = 0; j < npatitions; j++) {
        buf = temp_buffers[recvproc[(j*maxswap) + iswap]][iswap]->hostData();
        nrecv = nsend[recvproc[(j*maxswap) + iswap]];

        /* unpack buffer */
        n = atom[j].nlocal + atom[j].nghost;
        m = 0;
        for (i = 0; i < nrecv; i++) {
          m += atom[j].unpack_border(n++,&buf[m]);
        }
          

        /* set all pointers & counters */
        sendnum[(j*maxswap) + iswap] = nsend[j];
        recvnum[(j*maxswap) + iswap] = nrecv;
        firstrecv[(j*maxswap) + iswap] = atom[j].nlocal + atom[j].nghost;
        atom[j].nghost += nrecv;
      }
      
      iswap++;
    }
  }

  for(j = 0; j < npatitions; j++){
    d_sendlist[j]->upload();
  }

  delete[] nlast;
  delete[] nfirst;
  delete[] nsend;
}

/* realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA */

MMD_float* Comm::growsend(int n, int partition, int iswap)
{
  int old_size = maxsend[partition];
  maxsend[partition] = static_cast<int>(BUFFACTOR * n);

  for(int i = 0; i < nswap; i++) {
    mcl_unregister_buffer(temp_buffers[partition][i]->devData());
    cMCLData<MMD_float, xx>* temp = new cMCLData<MMD_float, xx>(mcl, MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_DYNAMIC, maxsend[partition], 0, 0);
    std::copy(temp_buffers[partition][i]->hostData(), temp_buffers[partition][i]->hostData() + old_size, temp->hostData());
    delete temp_buffers[partition][i];
    temp_buffers[partition][i] = temp;
  }
  return temp_buffers[partition][iswap]->hostData();
}

/* realloc the size of the iswap sendlist as needed with BUFFACTOR */

int** Comm::growlist(int iswapa, int n, int partition)
{
	unsigned int* dim=d_sendlist[partition]->getDim();
	unsigned int maxswap=dim[0];

	for(int iswap=0;iswap<dim[0];iswap++)
	{
    maxsendlist[(partition * maxswap) + iswap] = static_cast<int>(BUFFACTOR * n);
    sendlist[partition][iswap] = 
      (int *) realloc(sendlist[partition][iswap],maxsendlist[(partition * maxswap) + iswap]*sizeof(int));
	}
  mcl_unregister_buffer(d_sendlist[partition]->devData());
	delete d_sendlist[partition];
  d_sendlist[partition] = new cMCLData<int,xy>(mcl,(int*)sendlist[partition], MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_DYNAMIC | MCL_ARG_INPUT , maxswap,maxsendlist[(partition * maxswap)], 0);
  return sendlist[partition];
}

void Comm::free()
{
    for(auto hdl : to_free) mcl_hdl_free(hdl);
}

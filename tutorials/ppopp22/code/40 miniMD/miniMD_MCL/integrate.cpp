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
#include "integrate.h"
#include <minos.h>
#include <queue>
#include <cstring>
#include <chrono>
#include <thread>

#define NUM_SHARED_BUF 100
using namespace std;

Integrate::Integrate() {}
Integrate::~Integrate() {}

void Integrate::setup(int partitions)
{
    dtforce = 0.5 * dt;
    integrate_init_hdls = new mcl_handle *[partitions];
    force_hdls = new mcl_handle *[partitions];
    integrate_final_hdls = new mcl_handle *[partitions];
    neighbor_hdls = new mcl_handle *[partitions];
    
    for (int j = 0; j < partitions; j++)
    {
        integrate_init_hdls[j] = NULL;
        force_hdls[j] = NULL;
        integrate_final_hdls[j] = NULL;
        neighbor_hdls[j] = NULL;
    }
}

void Integrate::run(Atom atom[], Force &force, Neighbor neighbor[],
                    Comm &comm, Thermo &thermo, Timer &timer, int partitions, int share)
{
    mcl_handle** comm_hdls;

    MMD_float3** shared_mem;
    mcl_handle** share_hdls;
    char* shared_buf_name;
    if(share){
        shared_mem = new MMD_float3*[NUM_SHARED_BUF];
        share_hdls = new mcl_handle*[ntimes * partitions];
        shared_buf_name = (char*)malloc(strlen("mcl_pos_buffer_") + 7);
    }
    int nwait = 0;
    mcl_handle **waitlist = NULL;
    for (int n = 0; n < ntimes; n += neighbor[0].every)
    {
        for (int i = 0; i < neighbor[0].every - 1; i++)
        {
            //fprintf(stderr, "Starting iteration %d:%d\n", n, i);
            for (int j = 0; j < partitions; j++)
            {
                if(share && (n > 0) && (i > 0))
                {
                    nwait = 1;
                    waitlist = &share_hdls[((n + i - 1) * partitions) + j];

                }
                else if ((n > 0) && (i > 0))
                {
                    nwait = 1;
                    waitlist = &integrate_final_hdls[j];
                }
                else
                {
                    nwait = 0;
                    waitlist = NULL;
                }
                integrate_init_hdls[j] = mcl->LaunchKernel("integrate_kernel.h", "integrate_initial", atom[j].nlocal, nwait, waitlist, 7,
                                                           atom[j].d_x->devData(), atom[j].d_x->devSize(), atom[j].d_x->mclFlags(),
                                                           atom[j].d_v->devData(), atom[j].d_v->devSize(), atom[j].d_v->mclFlags(),
                                                           atom[j].d_f->devData(), atom[j].d_f->devSize(), atom[j].d_f->mclFlags(),
                                                           &atom[j].nlocal, sizeof(atom[j].nlocal), MCL_ARG_SCALAR,
                                                           &dt, sizeof(dt), MCL_ARG_SCALAR,
                                                           &dtforce, sizeof(dtforce), MCL_ARG_SCALAR,
                                                           &atom[j].nmax, sizeof(atom[j].nmax), MCL_ARG_SCALAR);
            }

            timer.stamp();
            comm_hdls = comm.communicate(atom, i, integrate_init_hdls);
            timer.stamp(TIME_COMM);

            for (int j = 0; j < partitions; j++)
            {
                //mcl_hdl_free(integrate_init_hdls[j]);
                //integrate_init_hdls[j] = NULL;
                //if (force_hdls[j])
                //{
                //    mcl_hdl_free(force_hdls[j]);
                //    force_hdls[j] = NULL;
                //    mcl_hdl_free(integrate_final_hdls[j]);
                //    integrate_final_hdls[j] = NULL;
                //}

                force_hdls[j] = force.compute(atom[j], neighbor[j], comm.nswap, &comm_hdls[j * comm.maxswap]);
            }
            delete[] comm_hdls;

            for (int j = 0; j < partitions; j++)
            {
                nwait = 1;
                waitlist = &force_hdls[j];
                integrate_final_hdls[j] = mcl->LaunchKernel("integrate_kernel.h", "integrate_final", atom[j].nlocal, nwait, waitlist, 5,
                                                            atom[j].d_v->devData(), atom[j].d_v->devSize(), atom[j].d_v->mclFlags(),
                                                            atom[j].d_f->devData(), atom[j].d_f->devSize(), atom[j].d_f->mclFlags(),
                                                            &atom[j].nlocal, sizeof(atom[j].nlocal), MCL_ARG_SCALAR,
                                                            &dtforce, sizeof(dtforce), MCL_ARG_SCALAR,
                                                            &atom[j].nmax, sizeof(atom[j].nmax), MCL_ARG_SCALAR);
            }
            // if(thermo.nstat) {
            //   thermo.compute(n + i,atom,neighbor,force,timer,comm);
            // }
            if (share)
            {
                for (int j = 0; j < partitions; j++)
                {
                    int natoms = min(65536, atom[j].nlocal);
                    int idx = ((n+i) * partitions) + j;
                    int buf_idx = idx % NUM_SHARED_BUF;
                    if(idx < NUM_SHARED_BUF) {
                        sprintf(shared_buf_name, "mcl_pos_buffer_%05d", buf_idx);
                        shared_mem[idx]  = (MMD_float3*)mcl_get_shared_buffer(shared_buf_name, 65536 * 3 * sizeof(float), 
                                            MCL_ARG_DYNAMIC | MCL_SHARED_MEM_NEW | MCL_SHARED_MEM_DEL_OLD);
                    }
                    nwait = 1;
                    waitlist = &integrate_final_hdls[j];
                    share_hdls[idx] = mcl->LaunchKernelShared("share_kernel.h", "copy_atoms", atom[j].nlocal, nwait, waitlist, 3,
                                                            atom[j].d_x->devData(), atom[j].d_x->devSize(), atom[j].d_x->mclFlags(),
                                                            shared_mem[idx], natoms * 3 * sizeof(float),  MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_SHARED | MCL_ARG_DYNAMIC,
                                                            &natoms, sizeof(natoms), MCL_ARG_SCALAR);
                }
            }

        }
        //mcl_wait_all();
        //fprintf(stderr, "Starting iteration %d:%d\n", n, neighbor[0].every - 1);
        for (int j = 0; j < partitions; j++)
        {
            nwait = 1;
            waitlist = &integrate_final_hdls[j];
            integrate_init_hdls[j] = mcl->LaunchKernel("integrate_kernel.h", "integrate_initial", atom[j].nlocal, nwait, waitlist, 7,
                                                       atom[j].d_x->devData(), atom[j].d_x->devSize(), atom[j].d_x->mclFlags() | MCL_ARG_OUTPUT,
                                                       atom[j].d_v->devData(), atom[j].d_v->devSize(), atom[j].d_v->mclFlags() | MCL_ARG_OUTPUT,
                                                       atom[j].d_f->devData(), atom[j].d_f->devSize(), atom[j].d_f->mclFlags(),
                                                       &atom[j].nlocal, sizeof(atom[j].nlocal), MCL_ARG_SCALAR,
                                                       &dt, sizeof(dt), MCL_ARG_SCALAR,
                                                       &dtforce, sizeof(dtforce), MCL_ARG_SCALAR,
                                                       &atom[j].nmax, sizeof(atom[j].nmax), MCL_ARG_SCALAR);
        }
        //fprintf(stderr, "Finished enqueing tasks.\n");
        mcl_wait_all();
        timer.stamp();
        //fprintf(stderr, "Finished all integrate initial tasks.\n");

        for (int j = 0; j < partitions; j++)
        {
            mcl_hdl_free(integrate_init_hdls[j]);
            integrate_init_hdls[j] = NULL;
            mcl_hdl_free(force_hdls[j]);
            force_hdls[j] = NULL;
            mcl_hdl_free(integrate_final_hdls[j]);
            integrate_final_hdls[j] = NULL;
        }
        //fprintf(stderr, "Freed all handles.\n");

        for (int j = 0; j < partitions; j++)
        {
            atom[j].d_x->download();
            atom[j].d_v->download();
        }

        comm.exchange(atom);
        comm.borders(atom);
        timer.stamp(TIME_COMM);

        //fprintf(stderr, "Binning atoms.\n");
        for (int j = 0; j < partitions; j++)
        {
            atom[j].d_x->upload();
            atom[j].d_v->upload();
            neighbor[j].resize_buffers(atom[j]);
            neighbor_hdls[j] = neighbor[j].binatoms(atom[j]);
            pending_list.push(j);
        }

        while (!pending_list.empty())
        {
            int j = pending_list.front();
            pending_list.pop();
            mcl_wait(neighbor_hdls[j]);
            mcl_hdl_free(neighbor_hdls[j]);
            neighbor_hdls[j] = neighbor[j].resize_and_bin(atom[j]);
            if (neighbor_hdls[j])
                pending_list.push(j);
            else
            {
                neighbor_hdls[j] = neighbor[j].build(atom[j]);
            }
        }

        for (int j = 0; j < partitions; j++)
        {
            mcl_wait(neighbor_hdls[j]);
            mcl_hdl_free(neighbor_hdls[j]);
            neighbor_hdls[j] = neighbor[j].reneigh(atom[j]);
            if (neighbor_hdls[j])
                pending_list.push(j);
            else
            {
                force_hdls[j] = force.compute(atom[j], neighbor[j], 0, NULL);
            }
        }

        while (!pending_list.empty())
        {
            int j = pending_list.front();
            pending_list.pop();
            mcl_wait(neighbor_hdls[j]);
            mcl_hdl_free(neighbor_hdls[j]);
            neighbor_hdls[j] = neighbor[j].reneigh(atom[j]);
            if (neighbor_hdls[j])
                pending_list.push(j);
            else
            {
                force_hdls[j] = force.compute(atom[j], neighbor[j], 0, NULL);
            }
        }
        //fprintf(stderr, "All force handles enqueued.\n");
        timer.stamp(TIME_NEIGH);

        for (int j = 0; j < partitions; j++)
        {
            uint64_t output = n + 1 >= ntimes ? MCL_ARG_OUTPUT : 0;
            integrate_final_hdls[j] = mcl->LaunchKernel("integrate_kernel.h", "integrate_final", atom[j].nlocal, 1, &force_hdls[j], 5,
                                                        atom[j].d_v->devData(), atom[j].d_v->devSize(), atom[j].d_v->mclFlags() | MCL_ARG_REWRITE | output,
                                                        atom[j].d_f->devData(), atom[j].d_f->devSize(), atom[j].d_f->mclFlags() | output,
                                                        &atom[j].nlocal, sizeof(atom[j].nlocal), MCL_ARG_SCALAR,
                                                        &dtforce, sizeof(dtforce), MCL_ARG_SCALAR,
                                                        &atom[j].nmax, sizeof(atom[j].nmax), MCL_ARG_SCALAR);
        }
        timer.stamp(TIME_FORCE);

        if (share)
        {
            for (int j = 0; j < partitions; j++)
            {
                int natoms = min(65536, atom[j].nlocal);
                int idx = ((n+neighbor[0].every-1) * partitions) + j;
                int buf_idx = idx % NUM_SHARED_BUF;
                if(idx < NUM_SHARED_BUF) {
                    sprintf(shared_buf_name, "mcl_pos_buffer_%05d", buf_idx);
                    shared_mem[idx]  = (MMD_float3*)mcl_get_shared_buffer(shared_buf_name, 65536 * 3 * sizeof(float), 
                                        MCL_ARG_DYNAMIC | MCL_SHARED_MEM_NEW | MCL_SHARED_MEM_DEL_OLD);
                }
                nwait = 1;
                waitlist = &integrate_final_hdls[j];
                share_hdls[idx] = mcl->LaunchKernelShared("share_kernel.h", "copy_atoms", atom[j].nlocal, nwait, waitlist, 3,
                                                        atom[j].d_x->devData(), atom[j].d_x->devSize(), atom[j].d_x->mclFlags(),
                                                        shared_mem[idx], natoms * 3 * sizeof(float),  MCL_ARG_BUFFER | MCL_ARG_RESIDENT | MCL_ARG_SHARED | MCL_ARG_DYNAMIC,
                                                        &natoms, sizeof(natoms), MCL_ARG_SCALAR);
            }
        }

        // if(thermo.nstat) {
        //   thermo.compute(n + neighbor->every - 1,atom,neighbor,force,timer,comm);
        // }
    }
    mcl_wait_all();

    if(share)
    {
        for(int i = 0; i < ntimes * partitions; i++){
            mcl_free_shared_buffer(shared_mem[i]);
        }
        this_thread::sleep_for(chrono::seconds(1));
    }
    

    for (int j = 0; j < partitions; j++)
    {
        atom[j].d_x->download();
        atom[j].d_v->download();
        atom[j].d_f->download();
    }
}

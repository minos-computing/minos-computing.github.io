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
//#define CUDA_PROFILER_ENABLED
//#include <cuda_profiler_api.h>
#include "stdio.h"
#include "stdlib.h"
#include <cstring>
#include <queue>
#include "ljs.h"
#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "integrate.h"
#include "thermo.h"
#include "comm.h"
#include "timer.h"
#include "variant.h"
#include "mcl_wrapper.h"
#include "mcl_data.h"
#include "precision.h"
#include <unistd.h>

#define MAXLINE 256

int input(In &, const char*);
void create_box(Atom &, int, int, int, double);
int create_atoms(Atom &, int, int, int, double);
void create_velocity(double, Atom*, Thermo &, int);
void output(In &, Atom*, Force&, Neighbor*, Comm &,
            Thermo &, Integrate &, Timer &, int, int);
int read_lammps_data(MCLWrapper* mcl, Atom &atom, Comm &comm, Neighbor &neighbor, Integrate &integrate, Thermo &thermo, char* file, int units, int nparts);
void mcl_verify(int res, timespec start);

int main(int argc, char **argv)
{
  //Common miniMD settings
  In in;
  in.datafile = NULL;
  int nparts=16;                 //number partitions
  int num_threads=32;		      //number of Threads per Block threads
  int num_steps=-1;             //number of timesteps (if -1 use value from lj.in)
  int system_size=-1;           //size of the system (if -1 use value from lj.in)
  int check_safeexchange=0;     //if 1 complain if atom moves further than 1 subdomain length between exchanges
  int do_safeexchange=0;        //if 1 use safe exchange mode [allows exchange over multiple subdomains]
  int use_sse=0;                //setting for SSE variant of miniMD only
  int screen_yaml=0;            //print yaml output to screen also
  int yaml_output=0;            //print yaml output
  int halfneigh=0;              //1: use half neighborlist; 0: use full neighborlist; -1: use original miniMD version half neighborlist force
  char* input_file = NULL;
  int ghost_newton = 0;
  int skip_gpu = 999;
  int neighbor_size = -1;
  int workers = 1;
  int share = 0;

  //MCL specific
  int use_tex = 0;
  int threads_per_atom = 1;

  for(int i = 0; i < argc; i++) {
    if((strcmp(argv[i], "-i") == 0) || (strcmp(argv[i], "--input_file") == 0)) {
      input_file = argv[++i];
      continue;
    }
    if((strcmp(argv[i],"-np")==0)||(strcmp(argv[i],"--nparts")==0)) {nparts=atoi(argv[++i]); continue;}
    if((strcmp(argv[i],"-w")==0)||(strcmp(argv[i],"--workers")==0)) {workers=atoi(argv[++i]); continue;}
  }

  MCLWrapper* mcl = new MCLWrapper;
  mcl->Init(argc,argv,workers);

  int error = 0;
  if(input_file == NULL)
    error = input(in, "in.lj.miniMD");
  else
	error = input(in, input_file);

  if (error)
  {
	  exit(0);
  }

  for(int i=0;i<argc;i++)
  {
     if((strcmp(argv[i],"-t")==0)||(strcmp(argv[i],"--num_threads")==0)) {num_threads=atoi(argv[++i]); continue;}
     if((strcmp(argv[i],"-n")==0)||(strcmp(argv[i],"--nsteps")==0))  {num_steps=atoi(argv[++i]); continue;}
     if((strcmp(argv[i],"-s")==0)||(strcmp(argv[i],"--size")==0))  {system_size=atoi(argv[++i]); continue;}
     if((strcmp(argv[i],"--share")==0))  {share=1; continue;}
     if((strcmp(argv[i],"--half_neigh")==0))  {halfneigh=atoi(argv[++i]); continue;}
     if((strcmp(argv[i],"-sse")==0))  {use_sse=atoi(argv[++i]); continue;}
     if((strcmp(argv[i],"--check_exchange")==0))  {check_safeexchange=1; continue;}
     if((strcmp(argv[i],"-o")==0)||(strcmp(argv[i],"--yaml_output")==0))  {yaml_output=atoi(argv[++i]); continue;}
     if((strcmp(argv[i],"--yaml_screen")==0))  {screen_yaml=1; continue;}
     if((strcmp(argv[i], "-f") == 0) || (strcmp(argv[i], "--data_file") == 0)) {
       if(in.datafile == NULL) in.datafile = new char[1000];

       strcpy(in.datafile, argv[++i]);
       continue;
     }
     if((strcmp(argv[i], "-u") == 0) || (strcmp(argv[i], "--units") == 0)) {
       in.units = strcmp(argv[++i], "metal") == 0 ? 1 : 0;
       continue;
     }

     if((strcmp(argv[i], "-p") == 0) || (strcmp(argv[i], "--force") == 0)) {
       in.forcetype = strcmp(argv[++i], "eam") == 0 ? FORCEEAM : FORCELJ;
       continue;
     }
     if((strcmp(argv[i], "-gn") == 0) || (strcmp(argv[i], "--ghost_newton") == 0)) {
       ghost_newton = atoi(argv[++i]);
       continue;
     }
     if((strcmp(argv[i], "--skip_gpu") == 0)) {
       skip_gpu = atoi(argv[++i]);
       continue;
     }
     if((strcmp(argv[i], "-b") == 0) || (strcmp(argv[i], "--neigh_bins") == 0))  {
       neighbor_size = atoi(argv[++i]);
       continue;
     }
	 if((strcmp(argv[i],"-tex")==0)||(strcmp(argv[i],"--texture")==0)) {use_tex=atoi(argv[++i]); continue;}
     if((strcmp(argv[i],"-tpa")==0)) {threads_per_atom=atoi(argv[++i]); continue;}
     if((strcmp(argv[i],"-h")==0)||(strcmp(argv[i],"--help")==0))
     {
        printf("\n---------------------------------------------------------\n");
        printf("-------------" VARIANT_STRING "------------\n");
        printf("---------------------------------------------------------\n\n");

        printf("miniMD is a simple, parallel molecular dynamics (MD) code,\n"
               "which is part of the Mantevo project at Sandia National\n"
   	           "Laboratories ( http://www.mantevo.org ).\n"
	           "The original authors of MPI based miniMD are Steve Plimpton (sjplimp@sandia.gov) ,\n"
               "Paul Crozier (pscrozi@sandia.gov) with current versions \n"
               "written by Christian Trott (crtrott@sandia.gov).\n\n");
        printf("Commandline Options:\n");
        printf("\n  Execution configuration:\n");
        printf("\t-t / --num_threads <threads>: set number of threads per block (default 32)\n");
        printf("\t--half_neigh <int>:           use half neighborlists (default 0)\n"
               "\t                                0: full neighborlist\n"
               "\t                                1: half neighborlist (not supported in OpenCL variant)\n"
               "\t                               -1: original miniMD half neighborlist force \n"
               "\t                                   (not supported in OpenCL variant)\n");
        printf("\t-np / --nparts:               partition problem into grid of size nparts (default:1)\n");
        printf("\t-w  / --workers:              number of MCL workers to use (default:1)\n");
        printf("\t-sse <sse_version>:           use explicit sse intrinsics (use miniMD-SSE variant)\n");
        printf("\t-gn / --ghost_newton <int>:   set usage of newtons third law for ghost atoms\n"
               "\t                              (only applicable with half neighborlists)\n");
        printf("\n  Simulation setup:\n");
        printf("\t-i / --input_file <string>:   set input file to be used (default: in.lj.miniMD)\n");
        printf("\t-n / --nsteps <nsteps>:       set number of timesteps for simulation\n");
        printf("\t-s / --size <size>:           set linear dimension of systembox and neighbor bins\n");
        printf("\t-b / --neigh_bins <int>:      set linear dimension of neighbor bin grid\n");
        printf("\t-u / --units <string>:        set units (lj or metal), see LAMMPS documentation\n");
        printf("\t-p / --force <string>:        set interaction model (lj or eam)\n");
        printf("\t-f / --data_file <string>:    read configuration from LAMMPS data file\n");

        printf("\n  Miscelaneous:\n");
        printf("\t--check_exchange:             check whether atoms moved further than subdomain width\n");
        printf("\t--safe_exchange:              perform exchange communication with all MPI processes\n"
	           "\t                              within rcut_neighbor (outer force cutoff)\n");
        printf("\t--yaml_output <int>:          level of yaml output (default 0)\n");
        printf("\t--yaml_screen:                write yaml output also to screen\n");
        printf("\t-tex / --texture <int>:       use texture cache in force kernel (default 0)\n");
        printf("\t-h / --help:                  display this help message\n\n");
        printf("---------------------------------------------------------\n\n");

        exit(0);
     }
  }

  Atom* atom = new Atom[nparts];
  Neighbor* neighbor = new Neighbor[nparts];
  Force force;
  Integrate integrate;
  Thermo thermo;
  Comm comm;
  Timer timer;

  if(in.forcetype == FORCEEAM) {
	  printf("ERROR: " VARIANT_STRING " does not yet support EAM simulations. Exiting.\n");
	  exit(0);
  }
  if(ghost_newton!=0)
  {
    printf("ERROR: -ghost_newton %i is not supported in " VARIANT_STRING ". Exiting.\n",ghost_newton);
    exit(0);
  }
  if(halfneigh!=0)
  {
    printf("ERROR: -half_neigh %i is not supported in " VARIANT_STRING ". Exiting.\n",halfneigh);
    exit(0);
  }
  if(use_tex!=0)
  {
    printf("ERROR: -tex %i is currently broken. Exiting.\n",use_tex);
    exit(0);
  }
  if(use_sse)
  {
    #ifndef VARIANT_SSE
    printf("ERROR: Trying to run with -sse with miniMD reference version. Use SSE variant instead. Exiting.\n");
    exit(0);
    #endif
  }

  for(int i = 0; i < nparts; i++){
    atom[i].threads_per_atom = threads_per_atom;
    atom[i].use_tex = use_tex;
    atom[i].mcl = mcl;

    neighbor[i].halfneigh=halfneigh;
    neighbor[i].mcl = mcl;

    if(neighbor_size > 0) {
      neighbor[i].nbinx = neighbor_size;
      neighbor[i].nbiny = neighbor_size;
      neighbor[i].nbinz = neighbor_size;
    }

    if(neighbor_size < 0 && in.datafile == NULL) {
      MMD_float neighscale = 5.0 / 6.0;
      neighbor[i].nbinx = neighscale * in.nx;
      neighbor[i].nbiny = neighscale * in.ny;
      neighbor[i].nbinz = neighscale * in.ny;
    }

    if(neighbor_size < 0 && in.datafile)
      neighbor[i].nbinx = -1;

    if(neighbor[i].nbinx == 0) neighbor[i].nbinx = 1;

    if(neighbor[i].nbiny == 0) neighbor[i].nbiny = 1;

    if(neighbor[i].nbinz == 0) neighbor[i].nbinz = 1;

    neighbor[i].every = in.neigh_every;
    neighbor[i].cutneigh = in.neigh_cut;
  }

  mcl->blockdim = num_threads;
  comm.do_safeexchange=do_safeexchange;
  force.use_sse=use_sse;
  

  integrate.mcl = mcl;
  force.mcl = mcl;
  comm.mcl = mcl;

  if(num_steps > 0) in.ntimes = num_steps;

  if(system_size > 0) {
    in.nx = system_size;
    in.ny = system_size;
    in.nz = system_size;
  }

  integrate.ntimes = in.ntimes;
  integrate.dt = in.dt;
  force.cutforce = in.force_cut;
  thermo.nstat = in.thermo_nstat;

  printf("# Create System:\n");

  if(in.datafile && 0) {
    // read_lammps_data(atom, comm, neighbor, integrate, thermo, in.datafile, in.units);
    // MMD_float volume = atom[0].box.xprd * atom[0].box.yprd * atom[0].box.zprd;
    // in.rho = 1.0 * atom.natoms / volume;
    // force.setup();

  } else {
    for(int i = 0; i < nparts; i++){
      create_box(atom[i], in.nx, in.ny, in.nz, in.rho);
    }
    comm.setup(neighbor[0].cutneigh, atom, nparts);

    for(int i = 0; i < nparts; i++){
      neighbor[i].setup(atom[i]);
    }
    
    integrate.setup(nparts);

    force.setup();

    for(int i = 0; i < nparts; i++){
      create_atoms(atom[i], in.nx, in.ny, in.nz, in.rho);
    }
    thermo.setup(mcl, in.rho, integrate, atom[0], in.units, nparts);
    
    create_velocity(in.t_request, atom, thermo, nparts);
  }
  printf("# Done .... \n");

  fprintf(stdout, "# " VARIANT_STRING " output ...\n");
  fprintf(stdout, "# Systemparameters: \n");
  fprintf(stdout, "\t# Inputfile: %s\n", input_file == 0 ? "in.lj.miniMD" : input_file);
  fprintf(stdout, "\t# Datafile: %s\n", in.datafile ? in.datafile : "None");
  fprintf(stdout, "\t# ForceStyle: %s\n", in.forcetype == FORCELJ ? "LJ" : "EAM");
  fprintf(stdout, "\t# Units: %s\n", in.units == 0 ? "LJ" : "METAL");
  fprintf(stdout, "\t# Atoms: %i\n", atom[0].natoms);
  fprintf(stdout, "\t# System size: %2.2lf %2.2lf %2.2lf (unit cells: %i %i %i)\n", atom[0].box.xprd, atom[0].box.yprd, atom[0].box.zprd, in.nx, in.ny, in.nz);
  fprintf(stdout, "\t# Density: %lf\n", in.rho);
  fprintf(stdout, "\t# Force cutoff: %lf\n", force.cutforce);
  fprintf(stdout, "\t# Neigh cutoff: %lf\n", neighbor[0].cutneigh);
  fprintf(stdout, "\t# Half neighborlists: %i\n", neighbor[0].halfneigh);
  fprintf(stdout, "\t# Neighbor bins: %i %i %i\n", neighbor[0].nbinx, neighbor[0].nbiny, neighbor[0].nbinz);
  fprintf(stdout, "\t# Neighbor frequency: %i\n", neighbor[0].every);
  fprintf(stdout, "\t# Timestep size: %lf\n", integrate.dt);
  fprintf(stdout, "\t# Thermo frequency: %i\n", thermo.nstat);
  fprintf(stdout, "\t# Ghost Newton: %i\n", ghost_newton);
  fprintf(stdout, "\t# Use SSE intrinsics: %i\n", force.use_sse);
  fprintf(stdout, "\t# Do safe exchange: %i\n", comm.do_safeexchange);
  fprintf(stdout, "\t# Size of float: %li\n\n",sizeof(MMD_float));

  
  comm.exchange(atom);
  comm.borders(atom);

  // int count = 0;
  // for(int j = 0; j < nparts; j++){
  //   for(int i = 0; i < comm.nswap; i++) {
  //     int offset = i * comm.maxsendlist[(j*comm.maxswap)];
  //     int first = comm.firstrecv[(j * comm.maxswap) + i];
  //     for(int k = 0; k < comm.sendnum[(j*comm.maxswap) + i]; k++) {
  //       int idx = comm.d_sendlist[j]->devData()[offset + k];
  //       if(idx > atom[j].nlocal + atom[j].nghost){
  //         fprintf(stderr, "Error in send list!\n");
  //         return -1;
  //       }
  //       if(atom[j].d_x->devData()[idx].x == 0 || atom[j].d_x->devData()[first + k].x == 0){
  //         count += 1;
  //       }
  //     }
  //   }
  // }
  // fprintf(stderr, "Send list verified. Count: %d!\n", count);
  // return 0;

  mcl_handle** hdls = new mcl_handle*[nparts];
  std::queue<int> pending_list;
  for(int j = 0; j < nparts; j++){
    atom[j].d_x->upload();
    atom[j].d_v->upload();
    atom[j].d_vold->upload();
    neighbor[j].resize_buffers(atom[j]);
    hdls[j] = neighbor[j].binatoms(atom[j]);
    pending_list.push(j);
  }
  
  while(!pending_list.empty()){
    int j = pending_list.front();
    pending_list.pop();
    mcl_wait(hdls[j]);
    mcl_hdl_free(hdls[j]);
    hdls[j] = neighbor[j].resize_and_bin(atom[j]);
    if(hdls[j]) pending_list.push(j);
    else {
      hdls[j] = neighbor[j].build(atom[j]);
    } 
  }

  for(int j = 0; j < nparts; j++){
    mcl_wait(hdls[j]);
    mcl_hdl_free(hdls[j]);
    hdls[j] = neighbor[j].reneigh(atom[j]);
    if(hdls[j]) pending_list.push(j);
  }

  while(!pending_list.empty()){
    int j = pending_list.front();
    pending_list.pop();
    mcl_wait(hdls[j]);
    mcl_hdl_free(hdls[j]);
    hdls[j] = neighbor[j].reneigh(atom[j]);
    if(hdls[j]) pending_list.push(j);
  }

  printf("# Starting dynamics ...\n");
  printf("# Timestep T U P Time\n");
  
  //fprintf(stderr, "Starting thermo compute...");
  //thermo.compute(0,atom,neighbor,force,timer,comm);
  //fprintf(stderr, "Done.\n");
  
  for(int j = 0; j < nparts; j++){
    hdls[j] = force.compute(atom[j], neighbor[j], 0, NULL);
  }
  mcl_wait_all();

  for(int j = 0; j < nparts; j++){
    mcl_hdl_free(hdls[j]);
  }
  
  //cudaProfilerStart();
  timespec start;
  clock_gettime(CLOCK_MONOTONIC, &start);

  timer.start(TIME_TOTAL);
  integrate.run(atom,force,neighbor,comm,thermo,timer,nparts,share);
  timer.stop(TIME_TOTAL);

  mcl_verify(0, start);
  comm.free();
  //cudaProfilerStop();

  //thermo.compute(-1,atom,neighbor,force,timer,comm);

  int natoms = 0;
  for(int j = 0; j < nparts; j++){
    natoms += atom[j].nlocal;
  }
  double time_other=timer.array[TIME_TOTAL]-timer.array[TIME_FORCE]-timer.array[TIME_NEIGH]-timer.array[TIME_COMM];
  //printf("\n\n");
  //printf("# Performance Summary:\n");
  //printf("# block_size nsteps natoms t_total t_force t_neigh t_comm t_other performance grep_string t_extra\n");
  //printf("%i %i %i %f %lf %lf %lf %lf %lf %lf PERF_SUMMARY \n\n\n",
  //    num_threads,integrate.ntimes,natoms,
  //    timer.array[TIME_TOTAL],timer.array[TIME_FORCE],timer.array[TIME_NEIGH],timer.array[TIME_COMM],time_other,
  //    1.0*natoms*integrate.ntimes/timer.array[TIME_TOTAL],timer.array[TIME_TEST]);

  if(yaml_output)
  output(in,atom,force,neighbor,comm,thermo,integrate,timer,screen_yaml,nparts);

  delete mcl;
  return 0;
}

#define TEST_SUCCESS "\x1B[32mSUCCESS \x1B[0m"
#define TEST_FAILED   "\x1B[31mFAILED \x1B[0m"
#define BILLION 1000000000ULL
#define tdiff(end,start) BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec

void mcl_verify(int res, timespec start)
{
	time_t t;
    timespec end;

	clock_gettime(CLOCK_MONOTONIC, &end);
	t = time(NULL);
	printf("-------------------------------------------\n");
	printf("End time:       %s", ctime(&t));
	printf("Result:         %s\n",res? "Failure": "Success");
	printf("Start Time: %lld.%.9ld\n", (long long) start.tv_sec, start.tv_nsec);
	printf("End Time: %lld.%.9ld\n", (long long) end.tv_sec, end.tv_nsec);
	printf("Execution time: %f seconds\n", ((float)tdiff(end,start))/BILLION);
	printf("===========================================\n");
}

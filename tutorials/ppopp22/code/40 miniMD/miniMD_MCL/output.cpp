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
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "ljs.h"
#include "atom.h"
#include "integrate.h"
#include "force.h"
#include "neighbor.h"
#include "comm.h"
#include "thermo.h"
#include "timer.h"
#include <time.h>
#include "variant.h"

void stats(int, double*, double*, double*, double*, int, int*);

void output(In &in, Atom atom[], Force& force, Neighbor neighbor[], Comm &comm,
            Thermo &thermo, Integrate &integrate, Timer &timer, int nparts, int screen_yaml)
{
  int i, n;
  int histo[10];
  double tmp, ave, max, min, total;
  FILE* fp;

  /* enforce PBC, then check for lost atoms */
  int natoms = 0, nlost = 0;
  for(int j = 0; j < nparts; j++){
    atom[i].pbc();
    natoms += atom[i].nlocal;

    for(i = 0; i < atom[i].nlocal; i++) {
      if(atom[i].x[i].x < 0.0 || atom[i].x[i].x >= atom[i].box.xprd ||
          atom[i].x[i].y < 0.0 || atom[i].x[i].y >= atom[i].box.yprd ||
          atom[i].x[i].z < 0.0 || atom[i].x[i].z >= atom[i].box.zprd) 
            nlost++;
    }
  }

  if(natoms != atom[0].natoms || nlost > 0) {
    printf("Atom counts = %d %d %d\n", nlost, natoms, atom[0].natoms);
    printf("ERROR: Incorrect number of atoms\n");
    return;
  }

  /* long-range energy and pressure corrections Whats this???*/

  double engcorr = 8.0 * 3.1415926 * in.rho *
                   (1.0 / (9.0 * pow(force.cutforce, MMD_float(9.0))) - 1.0 / (3.0 * pow(force.cutforce, MMD_float(3.0))));
  double prscorr = 8.0 * 3.1415926 * in.rho * in.rho *
                   (4.0 / (9.0 * pow(force.cutforce, MMD_float(9.0))) - 2.0 / (3.0 * pow(force.cutforce, MMD_float(3.0))));

  /* thermo output */

  double conserve;
  time_t general_time = time(NULL);
  struct tm local_time = *localtime(&general_time);
  char filename[256];

  sprintf(filename, "miniMD-%4d-%02d-%02d-%02d-%02d-%02d.yaml",
          local_time.tm_year + 1900, local_time.tm_mon + 1, local_time.tm_mday,
          local_time.tm_hour, local_time.tm_min, local_time.tm_sec);

  fp = fopen(filename, "w");

  if(screen_yaml) {
    fprintf(stdout, "run_configuration: \n");
    fprintf(stdout, "  variant: " VARIANT_STRING "\n");
    fprintf(stdout, "  datafile: %s\n", in.datafile ? in.datafile : "None");
    fprintf(stdout, "  units: %s\n", in.units == 0 ? "LJ" : "METAL");
    fprintf(stdout, "  atoms: %i\n", atom[0].natoms);
    fprintf(stdout, "  system_size: %2.2lf %2.2lf %2.2lf\n", atom[0].box.xprd, atom[0].box.yprd, atom[0].box.zprd);
    fprintf(stdout, "  unit_cells: %i %i %i\n", in.nx, in.ny, in.nz);
    fprintf(stdout, "  density: %lf\n", in.rho);
    fprintf(stdout, "  force_type: %s\n", in.forcetype == FORCELJ ? "LJ" : "EAM");
    fprintf(stdout, "  force_cutoff: %lf\n", force.cutforce);
    fprintf(stdout, "  neighbor_cutoff: %lf\n", neighbor[0].cutneigh);
    fprintf(stdout, "  neighbor_type: %i\n", neighbor[0].halfneigh);
    fprintf(stdout, "  neighbor_bins: %i %i %i\n", neighbor[0].nbinx, neighbor[0].nbiny, neighbor[0].nbinz);
    fprintf(stdout, "  neighbor_frequency: %i\n", neighbor[0].every);
    fprintf(stdout, "  timestep_size: %lf\n", integrate.dt);
    fprintf(stdout, "  thermo_frequency: %i\n", thermo.nstat);
    fprintf(stdout, "  ghost_newton: %i\n", neighbor[0].ghost_newton);
    fprintf(stdout, "  sse_intrinsics: %i\n", force.use_sse);
    fprintf(stdout, "  safe_exchange: %i\n", comm.do_safeexchange);
    fprintf(stdout, "  float_size: %li\n\n",sizeof(MMD_float));
  }

  fprintf(fp, "run_configuration: \n");
  fprintf(fp, "  variant: " VARIANT_STRING "\n");
  fprintf(fp, "  datafile: %s\n", in.datafile ? in.datafile : "None");
  fprintf(fp, "  units: %s\n", in.units == 0 ? "LJ" : "METAL");
  fprintf(fp, "  atoms: %i\n", atom[0].natoms);
  fprintf(fp, "  system_size: %2.2lf %2.2lf %2.2lf\n", atom[0].box.xprd, atom[0].box.yprd, atom[0].box.zprd);
  fprintf(fp, "  unit_cells: %i %i %i\n", in.nx, in.ny, in.nz);
  fprintf(fp, "  density: %lf\n", in.rho);
  fprintf(fp, "  force_type: %s\n", in.forcetype == FORCELJ ? "LJ" : "EAM");
  fprintf(fp, "  force_cutoff: %lf\n", force.cutforce);
  fprintf(fp, "  neighbor_cutoff: %lf\n", neighbor[0].cutneigh);
  fprintf(fp, "  neighbor_type: %i\n", neighbor[0].halfneigh);
  fprintf(fp, "  neighbor_bins: %i %i %i\n", neighbor[0].nbinx, neighbor[0].nbiny, neighbor[0].nbinz);
  fprintf(fp, "  neighbor_frequency: %i\n", neighbor[0].every);
  fprintf(fp, "  timestep_size: %lf\n", integrate.dt);
  fprintf(fp, "  thermo_frequency: %i\n", thermo.nstat);
  fprintf(fp, "  ghost_newton: %i\n", neighbor[0].ghost_newton);
  fprintf(fp, "  sse_intrinsics: %i\n", force.use_sse);
  fprintf(fp, "  safe_exchange: %i\n", comm.do_safeexchange);
  fprintf(fp, "  float_size: %li\n\n",sizeof(MMD_float));

  if(screen_yaml)
    fprintf(stdout, "\n\nthermodynamic_output:\n");

  fprintf(fp, "\n\nthermodynamic_output:\n");

  for(i = 0; i < thermo.mstat; i++) {
    conserve = (1.5 * thermo.tmparr[i] + thermo.engarr[i]) /
                (1.5 * thermo.tmparr[0] + thermo.engarr[0]);

    if(screen_yaml) {
      fprintf(stdout, "  timestep: %d \n", thermo.steparr[i]);
      fprintf(stdout, "      T*:           %15.10g \n", thermo.tmparr[i]);
      fprintf(stdout, "      U*:           %15.10g \n", thermo.engarr[i]);
      fprintf(stdout, "      P*:           %15.10g \n", thermo.prsarr[i]);
      fprintf(stdout, "      Conservation: %15.10g \n", conserve);
    }

    fprintf(fp    , "  timestep: %d \n", thermo.steparr[i]);
    fprintf(fp    , "      T*:           %15.10g \n", thermo.tmparr[i]);
    fprintf(fp    , "      U*:           %15.10g \n", thermo.engarr[i]);
    fprintf(fp    , "      P*:           %15.10g \n", thermo.prsarr[i]);
    fprintf(fp    , "      Conservation: %15.10g \n", conserve);
  }
  /* performance output */

  fprintf(stdout, "\n\n");
  fprintf(fp, "\n\n");

  double time_total = timer.array[TIME_TOTAL];
  double mflops = 4.0 / 3.0 * 3.1415926 *
                  pow(force.cutforce, MMD_float(3.0)) * in.rho * 0.5 *
                  23 * natoms * integrate.ntimes / time_total / 1000000.0;

  if(screen_yaml) {
    fprintf(stdout, "time:\n");
    fprintf(stdout, "  total:\n");
    fprintf(stdout, "    time: %g \n", time_total);
    fprintf(stdout, "    performance: %10.5e \n", natoms * integrate.ntimes / time_total);
  }

  fprintf(fp,    "time:\n");
  fprintf(fp,    "  total:\n");
  fprintf(fp,    "    time: %g \n", time_total);
  fprintf(fp,    "    performance: %10.5e \n", natoms * integrate.ntimes / time_total);

  double time_force = timer.array[TIME_FORCE];
  if(screen_yaml)
    fprintf(stdout, "  force: %g\n", time_force);
  fprintf(fp,    "  force: %g\n", time_force);

  double time_neigh = timer.array[TIME_NEIGH];
  if(screen_yaml)
    fprintf(stdout, "  neigh: %g\n", time_neigh);
  fprintf(fp,    "  neigh: %g\n", time_neigh);

  double time_comm = timer.array[TIME_COMM];
  if(screen_yaml)
    fprintf(stdout, "  comm:  %g\n", time_comm);
  fprintf(fp,    "  comm:  %g\n", time_comm);


  double time_other = time_total - (time_force + time_neigh + time_comm);
  if(screen_yaml)
    fprintf(stdout, "  other: %g\n", time_other);
  fprintf(fp,    "  other: %g\n", time_other);

  if(screen_yaml)
    fprintf(stdout, "\n");
  fprintf(fp, "\n");

  fclose(fp);
}
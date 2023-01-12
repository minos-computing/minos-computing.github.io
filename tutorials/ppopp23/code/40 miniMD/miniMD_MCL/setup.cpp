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

#include <cstdio>
#include <cmath>
#include "atom.h"
#include "thermo.h"
#include "precision.h"
#include "integrate.h"
#include "neighbor.h"

#include <cstring>
#include <cstdio>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

double random(int*);

#define NSECTIONS 3
#define MAXLINE 255
char line[MAXLINE];
char keyword[MAXLINE];
FILE* fp;

void read_lammps_parse_keyword(int first)
{
  int eof = 0;
  char buffer[MAXLINE];

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if(!first) {
    if(fgets(line, MAXLINE, fp) == NULL) eof = 1;
  }

  while(eof == 0 && strspn(line, " \t\n\r") == strlen(line)) {
    if(fgets(line, MAXLINE, fp) == NULL) eof = 1;
  }

  if(fgets(buffer, MAXLINE, fp) == NULL) eof = 1;

  // if eof, set keyword empty and return

  if(eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs


  // copy non-whitespace portion of line into keyword

  int start = strspn(line, " \t\n\r");
  int stop = strlen(line) - 1;

  while(line[stop] == ' ' || line[stop] == '\t'
        || line[stop] == '\n' || line[stop] == '\r') stop--;

  line[stop + 1] = '\0';
  strcpy(keyword, &line[start]);
}

void read_lammps_header(Atom& atom)
{
  int n;
  char* ptr;

  // customize for new sections

  const char* section_keywords[NSECTIONS] =
  {"Atoms", "Velocities", "Masses"};

  // skip 1st line of file

  char* eof = fgets(line, MAXLINE, fp);

  // customize for new header lines
  int ntypes = 0;

  while(1) {

    if(fgets(line, MAXLINE, fp) == NULL) n = 0;
    else n = strlen(line) + 1;

    if(n == 0) {
      line[0] = '\0';
      return;
    }

    // trim anything from '#' onward
    // if line is blank, continue

    double xlo, xhi, ylo, yhi, zlo, zhi;

    if(ptr = strchr(line, '#')) * ptr = '\0';

    if(strspn(line, " \t\n\r") == strlen(line)) continue;

    // search line for header keyword and set corresponding variable

    if(strstr(line, "atoms")) sscanf(line, "%i", &atom.natoms);
    else if(strstr(line, "atom types")) sscanf(line, "%i", &ntypes);

    // check for these first
    // otherwise "triangles" will be matched as "angles"

    else if(strstr(line, "xlo xhi")) {
      sscanf(line, "%lg %lg", &xlo, &xhi);
      atom.box.xprd = xhi - xlo;
    } else if(strstr(line, "ylo yhi")) {
      sscanf(line, "%lg %lg", &ylo, &yhi);
      atom.box.yprd = yhi - ylo;
    } else if(strstr(line, "zlo zhi")) {
      sscanf(line, "%lg %lg", &zlo, &zhi);
      atom.box.zprd = zhi - zlo;
    } else break;
  }

  // error check on total system size


  // check that exiting string is a valid section keyword

  read_lammps_parse_keyword(1);

  for(n = 0; n < NSECTIONS; n++)
    if(strcmp(keyword, section_keywords[n]) == 0) break;

  if(n == NSECTIONS) {
    char str[36 + 255 + 1];
    sprintf(str, "Unknown identifier in data file: %s", keyword);
  }

  // error check on consistency of header values
}

void read_lammps_atoms(Atom &atom, MMD_float3* x)
{
  int i;

  int nread = 0;
  int natoms = atom.natoms;
  atom.nlocal = 0;

  int type;
  double xx, xy, xz;

  while(nread < natoms) {
    fgets(line, MAXLINE, fp);
    sscanf(line, "%i %i %lg %lg %lg", &i, &type, &xx, &xy, &xz);
    i--;
    x[i].x = xx;
    x[i].y = xy;
    x[i].z = xz;
    nread++;
  }

}

void read_lammps_velocities(Atom &atom, MMD_float3* v)
{
  int i;

  int nread = 0;
  int natoms = atom.natoms;

  double x, y, z;

  while(nread < natoms) {
    fgets(line, MAXLINE, fp);
    sscanf(line, "%i %lg %lg %lg", &i, &x, &y, &z);
    i--;
    v[i].x = x;
    v[i].y = y;
    v[i].z = z;
    nread++;
  }

  // check that all atoms were assigned correctly

}

int read_lammps_data(MCLWrapper* mcl, Atom* atom, Comm &comm, Neighbor* neighbor, Integrate &integrate, Thermo &thermo, char* file, int units, int nparts)
{
  fp = fopen(file, "r");

  if(fp == NULL) {
    char str[128];
    sprintf(str, "Cannot open file %s", file);
  }

  read_lammps_header(atom[0]);
  for(int j = 0; j < nparts; j++){
    atom[j].box.xprd = atom[0].box.xprd;
    atom[j].box.yprd = atom[0].box.zprd;
    atom[j].box.zprd = atom[0].box.yprd;
    atom[j].natoms = atom[0].natoms;
  }

  comm.setup(neighbor[0].cutneigh, atom, nparts);

  for(int j = 0; j < nparts; j++){
    if(neighbor[j].nbinx < 0) {
      MMD_float volume = atom[j].box.xprd * atom[j].box.yprd * atom[j].box.zprd;
      MMD_float rho = 1.0 * atom[j].natoms / volume;
      MMD_float neigh_bin_size = pow(rho * 16, MMD_float(1.0 / 3.0));
      neighbor[j].nbinx = atom[j].box.xprd / neigh_bin_size;
      neighbor[j].nbiny = atom[j].box.yprd / neigh_bin_size;
      neighbor[j].nbinz = atom[j].box.zprd / neigh_bin_size;
    }
    if(neighbor[j].nbinx == 0) neighbor[j].nbinx = 1;

    if(neighbor[j].nbiny == 0) neighbor[j].nbiny = 1;

    if(neighbor[j].nbinz == 0) neighbor[j].nbinz = 1;

    neighbor[j].setup(atom[j]);
  }

  integrate.setup(nparts);

  //force->setup();

  thermo.setup(mcl, atom[0].box.xprd * atom[0].box.yprd * atom[0].box.zprd / atom[0].natoms, integrate, atom[0], units, nparts);

  MMD_float3* x = new MMD_float3[atom[0].natoms];
  MMD_float3* v = new MMD_float3[atom[0].natoms];

  int atomflag = 0;
  int tmp;

  while(strlen(keyword)) {
    if(strcmp(keyword, "Atoms") == 0) {
      read_lammps_atoms(atom[0], x);
      atomflag = 1;
    } else if(strcmp(keyword, "Velocities") == 0) {
      if(atomflag == 0) printf("Must read Atoms before Velocities\n");

      read_lammps_velocities(atom[0], v);
    } else if(strcmp(keyword, "Masses") == 0) {
      fgets(line, MAXLINE, fp);
      if(sizeof(MMD_float) == 4)
        sscanf(line, "%i %g", &tmp, &atom[0].mass);
      else
        sscanf(line, "%i %g", &tmp, &atom[0].mass);

      for(int j = 0; j < nparts; j++) atom[j].mass = atom[0].mass;
    }

    read_lammps_parse_keyword(0);
  }

  for(int i = 0; i < atom[0].natoms; i++) {
    for(int j = 0; j < nparts; j++){
      if(x[i].x >= atom[j].box.xlo && x[i].x < atom[j].box.xhi &&
        x[i].y >= atom[j].box.ylo && x[i].y < atom[j].box.yhi &&
        x[i].z >= atom[j].box.zlo && x[i].z < atom[j].box.zhi)
          atom[j].addatom(x[i].x, x[i].y, x[i].z, v[i].x, v[i].y, v[i].z);
    }
    
  }

  // check that all atoms were assigned correctly
  return 0;
}

/* create simulation box */

void create_box(Atom &atom, int nx, int ny, int nz, double rho)
{
  double lattice = pow((4.0 / rho), (1.0 / 3.0));
  atom.box.xprd = nx * lattice;
  atom.box.yprd = ny * lattice;
  atom.box.zprd = nz * lattice;
}

/* initialize atoms on fcc lattice in parallel fashion */

int create_atoms(Atom &atom, int nx, int ny, int nz, double rho)
{
  /* total # of atoms */

  atom.natoms = 4 * nx * ny * nz;
  atom.nlocal = 0;

  /* determine loop bounds of lattice subsection that overlaps my sub-box
     insure loop bounds do not exceed nx,ny,nz */

  double alat = pow((4.0 / rho), (1.0 / 3.0));
  int ilo = static_cast<int>(atom.box.xlo / (0.5 * alat) - 1);
  int ihi = static_cast<int>(atom.box.xhi / (0.5 * alat) + 1);
  int jlo = static_cast<int>(atom.box.ylo / (0.5 * alat) - 1);
  int jhi = static_cast<int>(atom.box.yhi / (0.5 * alat) + 1);
  int klo = static_cast<int>(atom.box.zlo / (0.5 * alat) - 1);
  int khi = static_cast<int>(atom.box.zhi / (0.5 * alat) + 1);

  ilo = MAX(ilo, 0);
  ihi = MIN(ihi, 2 * nx - 1);
  jlo = MAX(jlo, 0);
  jhi = MIN(jhi, 2 * ny - 1);
  klo = MAX(klo, 0);
  khi = MIN(khi, 2 * nz - 1);

  //fprintf(stderr, "Boundaries: x: (%e, %e), y: (%e, %e), z: (%e, %e)\n", atom.box.xlo, atom.box.xhi, atom.box.ylo, atom.box.yhi, atom.box.zlo, atom.box.zhi);
  /* each proc generates positions and velocities of atoms on fcc sublattice
       that overlaps its box
     only store atoms that fall in my box
     use atom # (generated from lattice coords) as unique seed to generate a
       unique velocity
     exercise RNG between calls to avoid correlations in adjacent atoms */

  double xtmp, ytmp, ztmp, vx, vy, vz;
  int i, j, k, m, n;
  int sx = 0;
  int sy = 0;
  int sz = 0;
  int ox = 0;
  int oy = 0;
  int oz = 0;
  int subboxdim = 8;

  int iflag = 0;

  int num_zeros = 0;
  while(oz * subboxdim <= khi) {
    k = oz * subboxdim + sz;
    j = oy * subboxdim + sy;
    i = ox * subboxdim + sx;

    if(iflag) continue;

    if(((i + j + k) % 2 == 0) &&
        (i >= ilo) && (i <= ihi) &&
        (j >= jlo) && (j <= jhi) &&
        (k >= klo) && (k <= khi)) {

      xtmp = 0.5 * alat * i;
      ytmp = 0.5 * alat * j;
      ztmp = 0.5 * alat * k;

      if(xtmp >= atom.box.xlo && xtmp < atom.box.xhi &&
          ytmp >= atom.box.ylo && ytmp < atom.box.yhi &&
          ztmp >= atom.box.zlo && ztmp < atom.box.zhi) {
        n = k * (2 * ny) * (2 * nx) + j * (2 * nx) + i + 1;

        for(m = 0; m < 5; m++) random(&n);

        vx = random(&n);

        for(m = 0; m < 5; m++) random(&n);

        vy = random(&n);

        for(m = 0; m < 5; m++) random(&n);

        vz = random(&n);

        if(xtmp == 0 || ytmp == 0 || ztmp == 0){
          num_zeros++;
        }
        atom.addatom(xtmp, ytmp, ztmp, vx, vy, vz);
      }
    }
    

    sx++;

    if(sx == subboxdim) {
      sx = 0;
      sy++;
    }

    if(sy == subboxdim) {
      sy = 0;
      sz++;
    }

    if(sz == subboxdim) {
      sz = 0;
      ox++;
    }

    if(ox * subboxdim > ihi) {
      ox = 0;
      oy++;
    }

    if(oy * subboxdim > jhi) {
      oy = 0;
      oz++;
    }
  }
  return 0;
}

/* adjust initial velocities to give desired temperature */

void create_velocity(double t_request, Atom* atom, Thermo &thermo, int nparts)
{
  int i, j;

  /* zero center-of-mass motion */

  double vxtot = 0.0;
  double vytot = 0.0;
  double vztot = 0.0;

  for(j = 0; j < nparts; j++){
    for(i = 0; i < atom[j].nlocal; i++) {
      vxtot += atom[j].v[i].x;
      vytot += atom[j].v[i].y;
      vztot += atom[j].v[i].z;
    }

    vxtot /= atom[j].natoms;
    vytot /= atom[j].natoms;
    vztot /= atom[j].natoms;

    for(i = 0; i < atom[j].nlocal; i++) {
      atom[j].v[i].x -= vxtot;
      atom[j].v[i].y -= vytot;
      atom[j].v[i].z -= vztot;
    }
  }

  /* rescale velocities, including old ones */
  thermo.t_act = 0;
  double t = thermo.temperature(atom, nparts);
  double factor = sqrt(t_request / t);

  for(j = 0; j < nparts; j++){
    for(i = 0; i < atom[j].nlocal; i++) {
      atom[j].v[i].x *= factor;
      atom[j].v[i].y *= factor;
      atom[j].v[i].z *= factor;
    }
  }
}

/* Park/Miller RNG w/out MASKING, so as to be like f90s version */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double random(int* idum)
{
  int k;
  double ans;

  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;

  if(*idum < 0) *idum += IM;

  ans = AM * (*idum);
  return ans;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

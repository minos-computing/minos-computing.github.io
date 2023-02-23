#include "precision.h"

__kernel void copy_atoms(__global MMD_floatK3* x, __global float* x_copy, int nlocal)
{
  int i = get_global_id(0);
  if(i<nlocal)
  {
    x_copy[(3*i)] = x[i].x;
    x_copy[(3*i)+1] = x[i].y;
    x_copy[(3*i)+2] = x[i].z;
  }

}
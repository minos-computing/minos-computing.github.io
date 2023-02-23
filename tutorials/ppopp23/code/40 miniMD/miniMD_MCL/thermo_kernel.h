#include "precision.h"

void warp_reduce_k3(__local MMD_floatK3* sdata, int tid){
    //sdata[tid] += sdata[tid + 32];
    sdata[tid] += sdata[tid + 16];
    sdata[tid] += sdata[tid + 8];
    sdata[tid] += sdata[tid + 4];
    sdata[tid] += sdata[tid + 2];
    sdata[tid] += sdata[tid + 1];
}

void warp_reduce_k2(__local float2* sdata, int tid){
    //sdata[tid] += sdata[tid + 32];
    sdata[tid] += sdata[tid + 16];
    sdata[tid] += sdata[tid + 8];
    sdata[tid] += sdata[tid + 4];
    sdata[tid] += sdata[tid + 2];
    sdata[tid] += sdata[tid + 1];
}

__kernel void energy_virial(__global MMD_floatK3* x, __global int* numneigh, __global int* neighbors, 
                            __global float2* sum, __local float2* temp, MMD_float cutforcesq, int maxneighs, int nlocal) 
{
    MMD_float sr2, sr6, phi, pair, rsq;
    MMD_floatK3 xi, delx;
    float2 ei;
    int i = get_global_id(0);

    int tid = get_local_id(0);
    int block_id = get_group_id(0);
    int block_dim = get_local_size(0);
    temp[tid] = (0.0f,0.0f);

    if(i<nlocal)
    {
        __global int* neighs = neighbors + i;
        xi = x[i];
        ei = (float2)(0.0f,0.0f);

        for (int k = 0; k < numneigh[i]; k++) {
            int j = neighs[k * nlocal];
            delx = (xi - x[j]);
            rsq = delx.x*delx.x + delx.y*delx.y + delx.z*delx.z;
            if (rsq < cutforcesq) {
                sr2 = 1.0f/rsq;
                sr6 = sr2*sr2*sr2;
                phi = sr6*(sr6-1.0f);
                pair = 48.0f * sr6 * (sr6 - 0.5f) * sr2;
                ei += (float2)(4.0f*phi, rsq * pair);
            }
        }
        temp[tid] = ei;
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    for(int s=block_dim/2; s>16; s = s>>1) {
        if (tid < s) {
            temp[tid] += temp[tid + s];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if(tid < 16) warp_reduce_k2(temp, tid);
    if (tid == 0) sum[block_id] = temp[0];
}

__kernel void temperature(__global const MMD_floatK3* v, __global MMD_float* sum, __local MMD_floatK3* temp, int nlocal) 
{
    int tid = get_local_id(0);
    int block_id = get_group_id(0);
    int block_dim = get_local_size(0);
    int grid_dim = get_global_size(0) * 2;
    int i = (block_id * block_dim * 2) + tid;

    temp[tid] = (MMD_floatK3)(0.0f, 0.0f, 0.0f);

    while(i < nlocal){
        temp[tid] += (v[i] * v[i]);
        temp[tid] += (v[i + block_dim] * v[i + block_dim]);
        i += grid_dim;
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    for(int s=block_dim/2; s>16; s=s>>1) {
        if (tid < s) {
            temp[tid] += temp[tid + s];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if(tid < 16) warp_reduce_k3(temp, tid);
    if (tid == 0) {
        sum[block_id] = temp[0].x + temp[0].y + temp[0].z;
    }
}
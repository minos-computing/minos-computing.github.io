__kernel void vec_add (__global float* out, __global float* x, __global float* y, int n)
{
        const int i = get_global_id(0);

        if (i < n) 
           out[i] = x[i] + y[i];
}
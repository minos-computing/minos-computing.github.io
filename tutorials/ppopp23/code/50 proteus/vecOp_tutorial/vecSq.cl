__kernel void vec_square (__global float* out, __global float* x, int n)
{
        const int i = get_global_id(0);

        if (i < n) 
           out[i] = x[i] * x[i];
}
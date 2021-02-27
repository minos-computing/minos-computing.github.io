#include <string.h>
#include <assert.h>

#include <minos.h>
#include "utils.h"

#define BLOCK_DIM 16

struct tc_struct_hdl{
	mcl_handle* hdl[4];
};

inline void transpose(mcl_handle** hdl, float* in, float* out, size_t n)
{
	int ret;
	size_t bsize = n * n * sizeof(float);
	int offset = 0;
	size_t szGlobalWorkSize[3] = { n, n, 1};
	size_t szLocalWorkSize[3]  = {BLOCK_DIM, BLOCK_DIM, 1};

	*hdl = mcl_task_create();
	assert(*hdl);
	
	ret = mcl_task_set_kernel(*hdl, "./transpose.cl", "transpose", 6, "", 0x0);
	assert(!ret);

	ret  = mcl_task_set_arg(*hdl, 0, (void*) out, bsize, MCL_ARG_OUTPUT | MCL_ARG_BUFFER);
	ret |= mcl_task_set_arg(*hdl, 1, (void*) in,  bsize, MCL_ARG_INPUT | MCL_ARG_BUFFER);
	ret |= mcl_task_set_arg(*hdl, 2, (void*) &offset, sizeof(int), MCL_ARG_INPUT | MCL_ARG_SCALAR);
	ret |= mcl_task_set_arg(*hdl, 3, (void*) &n, sizeof(int), MCL_ARG_INPUT | MCL_ARG_SCALAR);
	ret |= mcl_task_set_arg(*hdl, 4, (void*) &n, sizeof(int), MCL_ARG_INPUT | MCL_ARG_SCALAR);
	ret |= mcl_task_set_arg(*hdl, 5, NULL, (BLOCK_DIM + 1) * BLOCK_DIM * sizeof(float), MCL_ARG_LOCAL);
	assert(!ret);

	ret = mcl_exec(*hdl, szGlobalWorkSize, szLocalWorkSize, MCL_TASK_ANY);
	assert(!ret);

}

inline void gemm(mcl_handle** hdl, float* C, float* A, float* B, size_t n)
{
	int ret;
	size_t bsize = n * n * sizeof(float);
	size_t szGlobalWorkSize[3] = { n, n, 1};
	size_t szLocalWorkSize[3]  = {BLOCK_DIM, BLOCK_DIM, 1};

	*hdl = mcl_task_create();
	assert(*hdl);
	
	ret = mcl_task_set_kernel(*hdl, "./matrixMul.cl", "matrixMul", 8, "", 0x0);
	assert(!ret);

	ret  = mcl_task_set_arg(*hdl, 0, (void*) C, bsize, MCL_ARG_OUTPUT | MCL_ARG_BUFFER);
	ret |= mcl_task_set_arg(*hdl, 1, (void*) A, bsize, MCL_ARG_INPUT | MCL_ARG_BUFFER);
	ret |= mcl_task_set_arg(*hdl, 2, (void*) B, bsize, MCL_ARG_INPUT | MCL_ARG_BUFFER);
	ret |= mcl_task_set_arg(*hdl, 3, NULL, sizeof(float) * BLOCK_DIM *BLOCK_DIM, MCL_ARG_LOCAL);
	ret |= mcl_task_set_arg(*hdl, 4, NULL, sizeof(float) * BLOCK_DIM *BLOCK_DIM, MCL_ARG_LOCAL);
	ret |= mcl_task_set_arg(*hdl, 5, (void*) &n, sizeof(int), MCL_ARG_INPUT | MCL_ARG_SCALAR);
	ret |= mcl_task_set_arg(*hdl, 6, (void*) &n, sizeof(int), MCL_ARG_INPUT | MCL_ARG_SCALAR);
	ret |= mcl_task_set_arg(*hdl, 7, (void*) &n, sizeof(int), MCL_ARG_INPUT | MCL_ARG_SCALAR);
	assert(!ret);

	ret = mcl_exec(*hdl, szGlobalWorkSize, szLocalWorkSize, MCL_TASK_ANY);
	assert(!ret);

}

int main(int argc, char** argv)
{
	float *A, *B, *C;
	float *AT, *BT, *CT;
	unsigned long i, j;
	int ret;
	struct tc_struct_hdl* tc_hdl;
	
	mcl_banner("Tensor Contraction Skeleton");
	parse_global_opts(argc, argv);
	
	mcl_init(1, 0x0);

	A  = (float*) malloc(size * size * sizeof(float));
	AT = (float*) malloc(size * size * sizeof(float));
	B  = (float*) malloc(size * size * sizeof(float));
	BT = (float*) malloc(size * size * sizeof(float));
	C  = (float*) malloc(size * size * sizeof(float));
	CT = (float*) malloc(size * size * sizeof(float));
	tc_hdl = (struct tc_struct_hdl*) malloc (rep * sizeof(struct tc_struct_hdl));
	
	if(!A || !B || !C || !AT || !BT || !CT || !tc_hdl){
		printf("Error allocating vectors. Aborting.");
		ret = -1;
		goto err;
	}

	srand48(13579862);
	for(i=0; i<size; ++i){
		for(j=0; j<size; ++j){
			A[i*size+j] = (float)(0.5 + drand48()*1.5);
		}
	}
	
	for(i=0; i<size; ++i){
		for(j=0; j<size; ++j){
			B[i*size+j] = (float)(0.5 + drand48()*1.5);
		}
	}

	memset((char*) C, 0x0, size*size*sizeof(float));

	printf("-------------------------------------------\n");
	printf("\t Launching transposes...\n");
	for(i=0; i<rep; i++){
		transpose(&(tc_hdl[i].hdl[0]), A, AT, size);
		transpose(&(tc_hdl[i].hdl[1]), B, BT, size);
	}

	for(i=0; i<rep; i++){
		mcl_wait(tc_hdl[i].hdl[0]);
		mcl_wait(tc_hdl[i].hdl[1]);
		gemm(&(tc_hdl[i].hdl[2]),CT, AT, BT, size);
	}

	for(i=0; i<rep; i++){
		mcl_wait(tc_hdl[i].hdl[2]);
		transpose(&(tc_hdl[i].hdl[3]), CT, C, size);
	}

	
	mcl_wait_all();
	
	mcl_verify(0);

	free(tc_hdl);
	free(AT);
	free(BT);
	free(CT);
	free(C);
	free(B);
	free(A);

	mcl_finit();
	return 0;
 err:
	return -1;
}


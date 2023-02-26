#include <string.h>
#include <assert.h>

#include <minos.h>
#include "utils.h"

#define BLOCK_DIM 16

struct tc_struct_hdl{
	mcl_handle* hdl[4];
};

void square (mcl_handle** hdl, float* in, float* out, size_t n)
{
	int ret;
	size_t bsize = n * sizeof(float);
	int offset = 0;
	size_t szGlobalWorkSize[3] = { n, 1, 1};
	size_t szLocalWorkSize[3]  = {BLOCK_DIM, 1, 1};

	*hdl = mcl_task_create();
	assert(*hdl);
	
	ret = mcl_task_set_kernel(*hdl, "vec_square", 3);
	assert(!ret);

	ret  = mcl_task_set_arg(*hdl, 0, (void*) out, bsize, MCL_ARG_OUTPUT | MCL_ARG_BUFFER);
	ret |= mcl_task_set_arg(*hdl, 1, (void*) in,  bsize, MCL_ARG_INPUT | MCL_ARG_BUFFER);
	ret |= mcl_task_set_arg(*hdl, 2, (void*) &n, sizeof(int), MCL_ARG_SCALAR);
	assert(!ret);

	ret = mcl_exec(*hdl, szGlobalWorkSize, szLocalWorkSize, MCL_TASK_PROTEUS);
        //ret = mcl_exec(*hdl, szGlobalWorkSize, szLocalWorkSize, MCL_TASK_CPU);
        //ret = mcl_exec(*hdl, szGlobalWorkSize, szLocalWorkSize, MCL_TASK_ANY);
	assert(!ret);

}

void add (mcl_handle** hdl, float* C, float* A, float* B, size_t n)
{
	int ret;
	size_t bsize = n * sizeof(float);
	size_t szGlobalWorkSize[3] = { n, 1, 1};
	size_t szLocalWorkSize[3]  = {BLOCK_DIM, 1, 1};

	*hdl = mcl_task_create();
	assert(*hdl);
	
	ret = mcl_task_set_kernel(*hdl, "vec_add", 4);
	assert(!ret);

	ret  = mcl_task_set_arg(*hdl, 0, (void*) C, bsize, MCL_ARG_OUTPUT | MCL_ARG_BUFFER);
	ret |= mcl_task_set_arg(*hdl, 1, (void*) A, bsize, MCL_ARG_INPUT | MCL_ARG_BUFFER);
	ret |= mcl_task_set_arg(*hdl, 2, (void*) B, bsize, MCL_ARG_INPUT | MCL_ARG_BUFFER);
	ret |= mcl_task_set_arg(*hdl, 3, (void*) &n, sizeof(int), MCL_ARG_SCALAR);
	assert(!ret);

        ret = mcl_exec(*hdl, szGlobalWorkSize, szLocalWorkSize, MCL_TASK_PROTEUS);
	assert(!ret);

}

int main(int argc, char** argv)
{
	float *A, *B, *C;
	float *AT, *BT, *CT;
	unsigned long i;
	int ret;
	struct tc_struct_hdl* tc_hdl;
	
	mcl_banner("Vector Operator Skeleton");
	parse_global_opts(argc, argv);
	
	mcl_init(5, 0x0);

	mcl_prg_load("./vecSq.cl", "", MCL_PRG_SRC);
	mcl_prg_load("./vecAdd.cl", "", MCL_PRG_SRC);
	mcl_prg_load("vec_square", "", MCL_PRG_BUILTIN);   // the path should be name of builtin kernel 
	mcl_prg_load("vec_add", "", MCL_PRG_BUILTIN);

	A  = (float*) malloc(size * sizeof(float));
	AT = (float*) malloc(size * sizeof(float));
	B  = (float*) malloc(size * sizeof(float));
	BT = (float*) malloc(size * sizeof(float));
	C  = (float*) malloc(size * sizeof(float));
	CT = (float*) malloc(size * sizeof(float));
	tc_hdl = (struct tc_struct_hdl*) malloc (rep * sizeof(struct tc_struct_hdl));
	
	if(!A || !B || !C || !AT || !BT || !CT || !tc_hdl){
		printf("Error allocating vectors. Aborting.");
		ret = -1;
		goto err;
	}

	srand48(13579862);
	for(i=0; i<size; ++i){
		A[i] = (float)(0.5 + drand48()*1.5);
	}
	
	for(i=0; i<size; ++i){
		B[i] = (float)(0.5 + drand48()*1.5);
	}

	memset((char*) C, 0x0, size*sizeof(float));

	printf("-------------------------------------------\n");
	printf("\t Launching squares...\n");
	for(i=0; i<rep; i++){
		square(&(tc_hdl[i].hdl[0]), A, AT, size);
		square(&(tc_hdl[i].hdl[1]), B, BT, size);
	}

	for(i=0; i<rep; i++){
		mcl_wait(tc_hdl[i].hdl[0]);
		mcl_wait(tc_hdl[i].hdl[1]);
		add(&(tc_hdl[i].hdl[2]),CT, AT, BT, size);
	}

	for(i=0; i<rep; i++){
		mcl_wait(tc_hdl[i].hdl[2]);
		square(&(tc_hdl[i].hdl[3]), CT, C, size);
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


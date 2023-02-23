#include <unistd.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>

#include <minos.h>
#include "utils.h"



int test_mcl(double* A, double* B, double* C, size_t N)
{
	struct timespec start, end;
	mcl_handle*     hdl = NULL;
	uint64_t        pes[MCL_DEV_DIMS] = {N, N, 1};
	const size_t    msize    = N * N * sizeof(double);
	uint64_t        i;
	unsigned int    errs = 0;
	float           rtime;
	int             ret;
	
	printf("MCL Test...");

	clock_gettime(CLOCK_MONOTONIC,&start);

        if(mcl_prg_load("./gemmN.cl", "-DSINGLE_PRECISION", MCL_PRG_SRC)){
                printf("Error loading program %s. Aborting.\n", "./gemmN.cl");
                goto err;
        }

	hdl = mcl_task_create();
	if(!hdl){
	  printf("Error creating MCL task. Aborting.\n");
	  goto err;
	}

	if(mcl_task_set_kernel(hdl, "gemmN", 4)){
	  
	  printf("Error setting %s kernel. Aborting.\n", "gemmN");
	  goto err;
	}
	
	if(mcl_task_set_arg(hdl, 0, (void*) A, msize, MCL_ARG_INPUT|MCL_ARG_BUFFER|MCL_ARG_RESIDENT)){
	  printf("Error setting up task input A. Aborting.\n");
	  goto err;
	}
	
	if(mcl_task_set_arg(hdl, 1, (void*) B, msize, MCL_ARG_INPUT|MCL_ARG_BUFFER|MCL_ARG_RESIDENT)){
	  printf("Error setting up task input B. Aborting.\n");
	  goto err;
	}
	
	if(mcl_task_set_arg(hdl, 2, (void*) &N, sizeof(int), MCL_ARG_INPUT|MCL_ARG_SCALAR)){
	  printf("Error setting up task input N. Aborting.\n");
	  goto err;
	}
	
	if(mcl_task_set_arg(hdl, 3, (void*) C, msize, MCL_ARG_OUTPUT|MCL_ARG_BUFFER)){
	  printf("Error setting up task output. Aborting.\n");
	  goto err;
	}
	
	if((ret = mcl_exec(hdl, pes, NULL, MCL_TASK_ANY))){
	  printf("Error submitting task (%d)! Aborting.\n", ret);
	  goto err;
	}
	
	if(mcl_wait(hdl)){
	  printf("Request timed out!\n");
	  goto err;
	}
	
	clock_gettime(CLOCK_MONOTONIC, &end);

	if(hdl->ret == MCL_RET_ERROR){
	  printf("Error executing task %"PRIu64"!\n", i);
	  errs++;
	}
	if(errs)
		printf("Detected %u errors!\n",errs);
	else{
		rtime = ((float)tdiff(end,start))/BILLION;
		printf("Done.\n  Test time : %f seconds\n", rtime);
		printf("  Throughput: %f tasks/s\n", ((float)rep)/rtime);
	}
	
	mcl_hdl_free(hdl);
	
	return errs;

 err:
	return -1;
}

int main(int argc, char** argv)
{
	double          *A, *B, *C;
	int             i, j, ret = -1;
	
	mcl_banner("GEMM N Test");
	parse_global_opts(argc, argv);
	
	mcl_init(1,0x0);
	
	A = (double*) malloc(size * size * sizeof(double));
	B = (double*) malloc(size * size * sizeof(double));
	C = (double*) malloc(size * size * sizeof(double));
	
	if(!A || !B || !C){
		printf("Error allocating vectors. Aborting.");
		goto err;
	}

	srand48(13579862);
	for(i=0; i<size; ++i){
		for(j=0; j<size; ++j){
			A[i*size+j] = (double)(0.5 + drand48()*1.5);
		}
	}
	
	for(i=0; i<size; ++i){
		for(j=0; j<size; ++j){
			B[i*size+j] = (double)(0.5 + drand48()*1.5);
		}
	}
	

	ret = test_mcl(A,B,C,size);

	mcl_finit();

	free(A);
	free(B);
	free(C);
 err:
	return ret;
}

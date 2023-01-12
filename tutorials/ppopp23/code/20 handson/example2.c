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
	mcl_handle**    hdl = NULL;
	uint64_t        pes[MCL_DEV_DIMS] = {N, N, 1};
	const size_t    msize    = N * N * sizeof(double);
	uint64_t        i;
	unsigned int    errs = 0;
	float           rtime;
	int             ret;
	
	printf("MCL Test...");

	hdl = (mcl_handle**) malloc(sizeof(mcl_handle*) * rep);
	if(!hdl){
		printf("Error allocating memmory for MCL hanlders. Aborting.\n");
		goto err;
	}

        if(mcl_prg_load("./gemmN.cl", "-DSINGLE_PRECISION", MCL_PRG_SRC)){
                printf("Error loading program %s. Aborting.\n", "./gemmN.cl");
                goto err;
        }

	clock_gettime(CLOCK_MONOTONIC,&start);
	for(i=0; i<rep; i++){

		hdl[i] = mcl_task_create();
		if(!hdl[i]){
			printf("Error creating MCL task. Aborting.\n");
			continue;
		}
		if(mcl_task_set_kernel(hdl[i], "gemmN", 4)){
                        printf("Error setting %s kernel. Aborting.\n", "gemmN");
			continue;
		}
		
		if(mcl_task_set_arg(hdl[i], 0, (void*) A, msize, MCL_ARG_INPUT|MCL_ARG_BUFFER|MCL_ARG_RESIDENT)){
			printf("Error setting up task input A. Aborting.\n");
			continue;
		}
		
		if(mcl_task_set_arg(hdl[i], 1, (void*) B, msize, MCL_ARG_INPUT|MCL_ARG_BUFFER|MCL_ARG_RESIDENT)){
			printf("Error setting up task input B. Aborting.\n");
			continue;
		}

		if(mcl_task_set_arg(hdl[i], 2, (void*) &N, sizeof(int), MCL_ARG_INPUT|MCL_ARG_SCALAR)){
			printf("Error setting up task input N. Aborting.\n");
			continue;
		}

		if(mcl_task_set_arg(hdl[i], 3, (void*) C, msize, MCL_ARG_OUTPUT|MCL_ARG_BUFFER)){
			printf("Error setting up task output. Aborting.\n");
			continue;
		}
		
		if((ret = mcl_exec(hdl[i], pes, NULL, MCL_TASK_ANY))){
		  printf("Error submitting task (%d)! Aborting.\n", ret);
			continue;
		}

		if(mcl_wait(hdl[i])){
		       printf("Request timed out!\n");
		       continue;
		}
	}
	
	clock_gettime(CLOCK_MONOTONIC, &end);

	for(i=0; i<rep; i++)
		if(hdl[i]->ret == MCL_RET_ERROR){
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
	
	for(i=0; i<rep; i++)
		mcl_hdl_free(hdl[i]);
	
	free(hdl);
	
	return errs;

 err_hdl:
	free(hdl);
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

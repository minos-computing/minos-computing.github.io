#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <string.h>

#include <minos.h>

#define BILLION 1000000000ULL
#define tdiff(end,start) BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec

#define IMGSIZE 28*28

int readPGMFile(char* path, float* buffer, size_t buf_elem){
	FILE *in_file ;
	char ch;
	in_file = fopen(path, "r");
	if(in_file == NULL){
		fprintf(stderr,"\n Error opening pgm. Aborting.\n");
		return -1;
	}
	
	ch = fgetc(in_file);
	if(ch != 'P'){
		fprintf(stderr,"Error invalid pgm file (only accepts P5)\n");
		fclose(in_file);
		return -1;
	}

	ch = getc(in_file);
	int type = ch - 48; //convert char to single digit int;
	if (type != 5){
		fprintf(stderr,"Error invalid pgm file (only accepts P5)\n");
		fclose(in_file);
		return -1;
	}

	while(getc(in_file) != '\n');
	while(getc(in_file) == '#'){
		while(getc(in_file) != '\n');
	}
	fseek(in_file, -1, SEEK_CUR);             /* backup one character */
	int h,w,max __attribute__((unused)), ret __attribute__((unused));
	ret = fscanf(in_file,"%d",&h);
	ret = fscanf(in_file,"%d",&w);
	ret = fscanf(in_file,"%d",&max);


	if ((size_t) h*w != buf_elem){
		fprintf(stderr,"Error mismatched buffer sizes (pgm -- h: %d x w: %d = %d, supplied -- %ld)\n",h,w,h*w,buf_elem);
		fclose(in_file);
		return -1;
	}

	while (getc(in_file) != '\n');
	uint32_t i;
    printf("h: %d w: %d max: %d\n",h,w,max);
		
	for (i =0; i<buf_elem; i++){
		buffer[i] = (float)getc(in_file);
		printf("%c%c"," .:-=+*#%@"[(int)buffer[i] / 26 ],((i + 1) % 28) ? '\0' : '\n');
	}
	printf("\n");
	fclose(in_file);
	return 0;

}

int main(int argc, char** argv)
{
	uint64_t            i,r;
	mcl_handle**        hdls;
	struct timespec     start, end;
	float            	**in, **out;
	int                 ret = 0;
	uint64_t            pes[MCL_DEV_DIMS] = {1,1,1};
	unsigned int        errs, submitted;
	uint64_t 			num_inferences =2;
	uint64_t 			workers = 2;

	char* 				dla_bin = "mnist/mnist.nvdla";
	char*				image_dir = "mnist/";	
	

	hdls     = (mcl_handle**) malloc(num_inferences * sizeof(mcl_handle*));
	in       = (float**) malloc(10 * sizeof(float*));
	out      = (float**) malloc(num_inferences * sizeof(float*));
	
	if(!in || !out ){
		printf("Error allocating memory. Aborting.\n");
		goto err;
	}
	
	mcl_init(workers,0x0);
	int ndevs = mcl_get_ndev();

        printf("Current systems has total of %u devices.\n", ndevs);

	for(i=0;i<10;i++){
		in[i] = (float*) malloc(sizeof(float)*IMGSIZE);
		if(!in[i] ){
			printf("Error allocating memory. Aborting.\n");
			goto err;
		}
		char img_name[18];
		sprintf(img_name,"%s%lu.pgm",image_dir,i);
		if (readPGMFile(img_name,in[i],IMGSIZE) == -1){
			goto err;
		}
	}

	for(i=0;i<num_inferences;i++){
		out[i] = (float*) calloc(1,sizeof(float)* 10);
		if(!out[i] ){
			printf("Error allocating memory. Aborting.\n");
			goto err;
		}
	}

	printf("Synchronous Test...\n");
	clock_gettime(CLOCK_MONOTONIC,&start);
	submitted=0;
	errs=0;
	for (i=0;i<num_inferences;++i){		
		hdls[i] = mcl_task_create();
		
		if(!hdls[i]){
			printf("Error creating task %" PRIu64 ".\n",i);
			continue;
		}

		if(mcl_task_set_kernel(hdls[i], dla_bin, "DLA_MNIST", 2, "", MCL_KERNEL_BIN)){
			printf("Error setting task kernel %s for request %"PRIu64".\n","DLA_MNIST", i);
			continue;
		}
		int digit = rand()%10;
		printf("infering %d\n",digit);
		if(mcl_task_set_arg(hdls[i], 0, (void*) in[digit], sizeof(float)*IMGSIZE,
					MCL_ARG_INPUT | MCL_ARG_BUFFER)){
			printf("Error setting argument for task %"PRIu64".\n",i);
			continue;
		}		

		if(mcl_task_set_arg(hdls[i], 1, (void*) out[i], sizeof(float)*10,
					MCL_ARG_OUTPUT | MCL_ARG_BUFFER)){
			printf("Error setting output for task %"PRIu64".\n",i);
			continue;
		}		

		int ret = mcl_exec(hdls[i], pes, NULL, MCL_TASK_NVDLA);
		if(ret){
			printf("Error (%d) executing task %"PRIu64".",ret, i);
			continue;
		}
		
		submitted++;

		if(mcl_wait(hdls[i])){
			printf("Request %" PRIu64 " timed out!",i);
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &end);
	for (i=0;i<num_inferences;++i){
		printf("Prediction for request %lu\n",i);
		uint64_t j =0;
		for(j=0;j<10;j++){
			if (out[i][j] > 0.001){
				printf("%ld\n",j);
			}
			out[i][j] =0; //reset for async test
		}
		printf("--------------------------------------\n");
	}
	for(i=0; i<num_inferences; i++){
		if(hdls[i]->status != MCL_REQ_COMPLETED ){
			printf("Request %"PRIu64" status=%"PRIx64" return=%x"PRIx64"\n",
			       i,  hdls[i]->status, hdls[i]->ret);
			ret = 1;
			errs++;
		}
		mcl_hdl_free(hdls[i]);
	}
	
	if(!errs) {
		printf("Done.\n  Test time: %f seconds\n", ((float)tdiff(end,start))/BILLION);
	} else {
		printf("Detected %u errors!\n",errs);
	}

	printf("=============================================================================\n");

	printf("Asynchronous Test...\n");
	clock_gettime(CLOCK_MONOTONIC,&start);
	errs=0;
	submitted=0;
	for (i=0;i<num_inferences;++i){	
		hdls[i] = mcl_task_create();
		
		if(!hdls[i]){
			printf("Error creating task %" PRIu64 ".\n",i);
			continue;
		}

		if(mcl_task_set_kernel(hdls[i], dla_bin, "DLA_MNIST", 2, "", MCL_KERNEL_BIN)){
			printf("Error setting task kernel %s for request %"PRIu64".\n","DLA_MNIST", i);
			continue;
		}
		int digit = rand()%10;
                printf("infering %d\n",digit);
		if(mcl_task_set_arg(hdls[i], 0, (void*) in[digit], sizeof(float)*IMGSIZE,
					MCL_ARG_INPUT | MCL_ARG_BUFFER)){
			printf("Error setting argument for task %"PRIu64".\n",i);
			continue;
		}		

		if(mcl_task_set_arg(hdls[i], 1, (void*) out[i], sizeof(float)*10,
					MCL_ARG_OUTPUT | MCL_ARG_BUFFER)){
			printf("Error setting output for task %"PRIu64".\n",i);
			continue;
		}		

		int ret = mcl_exec(hdls[i], pes, NULL, MCL_TASK_NVDLA);
		if(ret){
			printf("Error (%d) executing task %"PRIu64".",ret, i);
			continue;
		}
		
		submitted++;
	}
	

	if(submitted == num_inferences){
	  	mcl_wait_all();
	} else{
	  	printf("Not all requests have been successfully submitted! (%u/%"PRIu64")\n", submitted,num_inferences);
	}
	
	clock_gettime(CLOCK_MONOTONIC, &end);
	for (i=0;i<num_inferences;++i){
		printf("Prediction for request %lu\n",i);
		uint64_t j =0;
		for(j=0;j<10;j++){
			if (out[i][j] > 0.001){
				printf("%ld\n",j);
			}
			out[i][j] =0; //reset for async test
		}
		printf("--------------------------------------\n");
	}
	for(i=0; i<num_inferences; i++){
		if(hdls[i]->status != MCL_REQ_COMPLETED ){
			printf("Request %"PRIu64" status=%"PRIx64" retrun=%x"PRIx64"\n",
			       i, hdls[i]->status, hdls[i]->ret);
			ret = 1;
			errs++;
		}
	}
	
	if(!errs){
		printf("Done.\n  Test time: %f seconds\n", ((float)tdiff(end,start))/BILLION);
	} else{
		printf("Detected %u errors!\n",errs);
	}
	printf("=============================================================================\n");
		
	for(i=0; i<num_inferences; i++){
		mcl_hdl_free(hdls[i]);
		free(out[i]);
	}
	for(i=0; i<10;i++){
		free(in[i]);
		
	}

	mcl_finit();
	
	free(out);
	free(in);
	free(hdls);
	
	return 0;
	
 err:
	exit(EXIT_FAILURE);
}

//-----------------example used in slides------------------------//

/* handle command line arguments, declare variables, etc*/

// char* dla_bin = "path to .nvdla binary"

// float **in; //an array of digit images convert to 1-D arrays of floats
// float **out; //an array of 10 element arrays specifying the predicted digit

// /* input/output/ omitted for slides, please see full source for details */

// mcl_handle** hdls =(mcl_handle**) malloc(num_inferences * size_of(mcl_handle*));
// mcl_init(workers,0x0);

// for (i=0;i<num_inferences;i++){
// 	hdls[i]=mcl_task_create();
// 	ret = mcl_task_set_kernel(hdls[i], dla_bin, "DLA_MNIST", 2, "", MCL_KERNEL_BIN));
// 	ret = mcl_task_set_arg(hdls[i], 0, (void*) in[i], sizeof(float)*IMGSIZE, MCL_ARG_INPUT | MCL_ARG_BUFFER));
// 	ret = mcl_task_set_arg(hdls[i], 1, (void*) out[i], sizeof(float)*10, MCL_ARG_OUTPUT | MCL_ARG_BUFFER));	
// 	ret =  mcl_exec(hdls[h_idx], pes, NULL, MCL_TASK_NVDLA);
// 	// ret = mcl_wait(hdls[i]); //uncomment for synchronous execution
// }
// mcl_wait_all();

// //print output predictions

// for(i=0;i<num_inferences; i++){
// 	mcl_hdl_free(hdls[i]);
// }

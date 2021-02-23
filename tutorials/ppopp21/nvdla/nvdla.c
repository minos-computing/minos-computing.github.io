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
	uint64_t 			num_digits =2;
	uint64_t 			workers = 2;

	char* 				dla_bin = "mnist/mnist.nvdla";
	char*				image_dir = "mnist/";	


	hdls     = (mcl_handle**) malloc(num_digits * sizeof(mcl_handle*));
	in       = (float**) malloc(num_digits * sizeof(float*));
	out      = (float**) malloc(num_digits * sizeof(float*));
	
	if(!in || !out ){
		printf("Error allocating memory. Aborting.\n");
		goto err;
	}
	
	mcl_init(workers,0x0);

	for(i=0;i<num_digits;i++){
		in[i] = (float*) malloc(sizeof(float)*IMGSIZE);
		out[i] = (float*) malloc(sizeof(float)* 10);
		if(!in[i] || !out[i] ){
			printf("Error allocating memory. Aborting.\n");
			goto err;
		}
		char img_name[18];
		sprintf(img_name,"%s%lu.pgm",image_dir,i);
		if (readPGMFile(img_name,in[i],IMGSIZE) == -1){
			goto err;
		}
		uint64_t j;
		for (j = 0; j<10; j++){
			out[i][j] = 0;
		}
	}

	printf("Synchronous Test...\n");
	clock_gettime(CLOCK_MONOTONIC,&start);
	submitted=0;
	errs=0;
	for (i=0;i<num_digits;++i){		
		uint64_t h_idx=i;
		hdls[h_idx] = mcl_task_create();
		
		if(!hdls[h_idx]){
			printf("Error creating task %" PRIu64 ".\n",i);
			continue;
		}

		if(mcl_task_set_kernel(hdls[h_idx], dla_bin, "DLA_MNIST", 2, "", MCL_KERNEL_BIN)){
			printf("Error setting task kernel %s for request %"PRIu64".\n","DLA_MNIST", h_idx);
			continue;
		}

		if(mcl_task_set_arg(hdls[h_idx], 0, (void*) in[i], sizeof(float)*IMGSIZE,
					MCL_ARG_INPUT | MCL_ARG_BUFFER)){
			printf("Error setting argument for task %"PRIu64".\n",h_idx);
			continue;
		}		

		if(mcl_task_set_arg(hdls[h_idx], 1, (void*) out[i], sizeof(float)*10,
					MCL_ARG_OUTPUT | MCL_ARG_BUFFER)){
			printf("Error setting output for task %"PRIu64".\n",h_idx);
			continue;
		}		

		int ret = mcl_exec(hdls[h_idx], pes, NULL, MCL_TASK_NVDLA);
		if(ret){
			printf("Error (%d) executing task %"PRIu64".",ret, h_idx);
			continue;
		}
		
		submitted++;

		if(mcl_wait(hdls[h_idx])){
			printf("Request %" PRIu64 " timed out!",h_idx);
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &end);
	for (i=0;i<num_digits;++i){
		printf("Predictions (prob > 0.001) for %lu\n",i);
		uint64_t j =0;
		for(j=0;j<10;j++){
			if (out[i][j] > 0.001){
				printf("%ld %f\n",j,out[i][j]);
			}
			out[i][j] =0; //reset for async test
		}
		printf("--------------------------------------\n");
	}
	for(i=0; i<num_digits; i++){
		printf("i: %"PRIx64" %p\n",i,hdls[i]);
		if(hdls[i]->status != MCL_REQ_COMPLETED ){
			printf("Request %"PRIu64" status=%"PRIx64" return=%x"PRIx64"\n",
			       i,  hdls[i]->status, hdls[i]->ret);
			ret = 1;
			errs++;
		}
	}
	
	if(!errs) {
		printf("Done.\n  Test time: %f seconds\n", ((float)tdiff(end,start))/BILLION);
	} else {
		printf("Detected %u errors!\n",errs);
	}

	for(i=0; i<num_digits; i++){
		mcl_hdl_free(hdls[i]);
	}

	printf("=============================================================================\n");

	printf("Asynchronous Test...\n");
	clock_gettime(CLOCK_MONOTONIC,&start);
	errs=0;
	submitted=0;
	for (i=0;i<num_digits;++i){		
		uint64_t h_idx=r*num_digits + i;
		hdls[h_idx] = mcl_task_create();
		
		if(!hdls[h_idx]){
			printf("Error creating task %" PRIu64 ".\n",i);
			continue;
		}

		if(mcl_task_set_kernel(hdls[h_idx], dla_bin, "DLA_MNIST", 2, "", MCL_KERNEL_BIN)){
			printf("Error setting task kernel %s for request %"PRIu64".\n","DLA_MNIST", h_idx);
			continue;
		}

		if(mcl_task_set_arg(hdls[h_idx], 0, (void*) in[i], sizeof(float)*IMGSIZE,
					MCL_ARG_INPUT | MCL_ARG_BUFFER)){
			printf("Error setting argument for task %"PRIu64".\n",h_idx);
			continue;
		}		

		if(mcl_task_set_arg(hdls[h_idx], 1, (void*) out[i], sizeof(float)*10,
					MCL_ARG_OUTPUT | MCL_ARG_BUFFER)){
			printf("Error setting output for task %"PRIu64".\n",h_idx);
			continue;
		}		

		int ret = mcl_exec(hdls[h_idx], pes, NULL, MCL_TASK_NVDLA);
		if(ret){
			printf("Error (%d) executing task %"PRIu64".",ret, h_idx);
			continue;
		}
		
		submitted++;
	}
	

	if(submitted == num_digits){
	  	mcl_wait_all();
	} else{
	  	printf("Not all requests have been successfully submitted! (%u/%"PRIu64")\n", submitted,num_digits);
	}
	
	clock_gettime(CLOCK_MONOTONIC, &end);
	for (i=0;i<num_digits;++i){
		printf("Predictions (prob > 0.001) for %lu\n",i);
		uint64_t j =0;
		for(j=0;j<10;j++){
			if (out[i][j] > 0.001){
				printf("%ld %f\n",j,out[i][j]);
			}
			out[i][j] =0; //reset for async test
		}
		printf("--------------------------------------\n");
	}
	for(i=0; i<num_digits; i++){
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
		
	for(i=0; i<num_digits; i++){
		mcl_hdl_free(hdls[i]);
		free(in[i]);
		free(out[i]);
	}

	mcl_finit();
	
	free(out);
	free(in);
	free(hdls);
	
	return 0;
	
 err:
	exit(EXIT_FAILURE);
}

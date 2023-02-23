#include <unistd.h> 
#include <iostream>
#include <cstdlib> 
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <iostream>
#include <fstream>
#include <iterator>

#include "Connection.h"
#include "Message.h"

#define ERROR_TAG "[NVDLA CLIENT]"

std::vector<uint8_t> loadBinary(char* file){
    std::ifstream infile;
    infile.open(file,std::ios::binary);
    infile.seekg(0,std::ios::end);
    int size = infile.tellg();
    infile.seekg(0,std::ios::beg);
    std::cout<<"loading binary image: "<<file<<" size: "<<size<<std::endl;
    std::vector<uint8_t> buffer(size);
    infile.read((char*)buffer.data(),size);
    infile.close();
    return buffer;
}

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
    std::cout<<"h: "<<h<<" w: "<<w<<" max: "<<max<<std::endl;
		
	for (i =0; i<buf_elem; i++){
		buffer[i] = (float)getc(in_file);
		printf("%c%c"," .:-=+*#%@"[(int)buffer[i] / 26 ],((i + 1) % 28) ? '\0' : '\n');
	}
	printf("\n");
	fclose(in_file);
	return 0;
}

// parsePGM(char* path,){
//     std::ifstream infile(path);
//     if(infile.fail()){
//         std::cerr<<"invalid file"<<std::endl;
//     }
//     int i =0;
//     std::string header[3];
//     while (i < 3){
//         std::getline(infile,header[i]);
//         if (header[i] != '#'){ //check for comments
//             i++;
//         }
//     }
//     if (header[i].compare("P5") != 0){
//         std::cerr<< "invalid pgm magic value: "<<header[0]<<std::endl;
//     }
//     uint32_t width;
//     uint32_t height;
//     std::stringstream ss(header[1]);
//     if (!(ss >> width >> height)){
//         std::cerr<< "invalid pgm: "<<header[1]<<std::endl;
//     }
//     ss.clear();

//     uint32_t maxVal;
//     ss.str(header[2]);
//     if (!(ss >> maxVal)){
//         std::cerr<< "invalid pgm: "<<header[2]<<std::endl;
//     }

//     std::cout<<"h: "<<height<<" w: "<<weight<<" max: "<<maxVal<<std::endl;




    
// }



int main(int argc, char *argv[]){
    if (argc < 2){
        std::cerr<<"ERROR enter nvdla binary image as argument"<<std::endl;
        exit(1);
    }

    Connection con("127.0.0.1",6667);
    con.lock();
    sendInitMsg(&con);
    int64_t num_devices = getNumDevices(&con);
    int64_t max_devices = getMaxDevices(&con);
    std::cout<<"max devices "<< max_devices<<" num_devices "<<num_devices<<std::endl;
    auto buf = loadBinary(argv[1]);
    std::cout<<"sending binary image: "<<argv[1]<<" size: "<<buf.size()<<std::endl;
    sendLoadBinary(&con,buf.data(),buf.size(),std::string(argv[1]));
    float* buffer = new float[28*28];
    readPGMFile(argv[2], buffer, 28*28);
    sendInputData(&con,(uint8_t*)buffer,sizeof(float)*28*28);
    std::cout<<"sent buffer"<<std::endl;
    float* outBuffer = new float[10];
    bool ret =sendInfer(&con,(uint8_t*)outBuffer,sizeof(float)*10);
    std::cout<<"sent infer "<<ret<<std::endl;
    con.unlock();

	for (int i=0;i<10; i++){
		std::cout<<"i: "<<i<<" "<<outBuffer[i]<<std::endl;
	}
    
    delete[] outBuffer;
    delete[] buffer;
    return 0; 
}


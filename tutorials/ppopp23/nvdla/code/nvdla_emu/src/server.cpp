#include <unistd.h> 
#include <iostream>
#include <cstdlib> 
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <cstring>
#include <cerrno>





#include "nvdla/IRuntime.h"
#include "nvdla_inf.h"
#include "nvdla_os_inf.h"
#include "half.h"

#include "DlaUtils.h"

#include "dlatypes.h"
#include "dlaerror.h"

#include "ErrorMacros.h"

#include "WorkPool.h"
#include "Connection.h"
#include "Message.h"

#define ERROR_TAG "[NVDLA SERVER ERROR] "

WorkPool<std::function<void()>> workPool(1); //since this is running is qemu keep single threaded for now



 NvDlaError init_nvdla(NvdlaContext *ctx){
    NvDlaError e = NvDlaSuccess;
    ctx->rt = nvdla::createRuntime();
    std::cout<<"created runtime "<<(void*) ctx->rt<<std::endl;
    if (ctx->rt == NULL){
        ORIGINATE_ERROR(NvDlaError_BadParameter, "createRuntime() failed");
        sendAckMsg(ctx->connection,FAIL_MSG);
    }
    sendAckMsg(ctx->connection,INIT);
    return e;
}


void send_device_info(NvdlaContext *ctx, msgType type){
    if (ctx->rt){
        sendDevices(ctx->connection,type,ctx->rt->getNumDevices());
    }
    else{
        std::cout <<"warning nvlda not instantiated"<<std::endl;
        sendDevices(ctx->connection,type,-1);
    }
}


void print_tensor_info(nvdla::IRuntime *rt, nvdla::IRuntime::NvDlaTensor &tensor, int i){
    std::cout<<"Tensor: "<<i<<" "<<tensor.name<<" buf size: "<<tensor.bufferSize<<std::endl;
    std::cout<<"dims: "<<tensor.dims.n<<" "<<tensor.dims.c<<" "<<tensor.dims.w<<" "<< tensor.dims.h<<std::endl;
    std::cout<<"Data Format "<<(uint32_t)tensor.dataFormat<<" type: "<<(uint32_t)tensor.dataType<<" category: "<<(uint32_t)tensor.dataCategory
            <<" pixel format: "<<(uint32_t)tensor.pixelFormat<<" mapping: "<<(uint32_t)tensor.pixelMapping<<std::endl;
    std::cout<<"stride: ";
    for (int j=0;j<NVDLA_RUNTIME_TENSOR_DESC_NUM_STRIDES;++j){
        std::cout<<j<<": "<<(uint32_t)tensor.stride[j]<<", ";
    }
    std::cout<<std::endl;
}
void load_binary_image(NvdlaContext *ctx, uint8_t* buff){
    dataMsg* msg = (dataMsg*) buff;
    std::string tag((char*)msg->data);
    uint8_t* data = msg->data +msg->tagSize;
    if (tag != ctx->network){
        if (ctx->network == ""){
            ctx->rt->unload();
        }
        if (ctx->rt->load(data,0)){
            ctx->rt->initEMU();
            ctx->network=tag;
            sendAckMsg(ctx->connection, LOAD_BINARY);
            std::cout<<"loaded binary image"<<std::endl;
            int num_tensors;
            ctx->rt->getNumInputTensors(&num_tensors);
            std::cout<<"num input tensors: "<<num_tensors<<std::endl;
            nvdla::IRuntime::NvDlaTensor tensor;
            for (int i =0; i< num_tensors; ++i){
                ctx->rt->getInputTensorDesc(i, &tensor);
                print_tensor_info(ctx->rt,tensor,i);
            }
            for (int i =0; i< num_tensors; ++i){
                ctx->rt->getOutputTensorDesc(i, &tensor);
                print_tensor_info(ctx->rt,tensor,i);
            }
        }
        else{
            sendAckMsg(ctx->connection,FAIL_MSG);
        }
    }
    else{
        sendAckMsg(ctx->connection, LOAD_BINARY);
    }
}

NvDlaError load_input_data(NvdlaContext *ctx, uint8_t* buff){
    NvDlaError e = NvDlaSuccess;
    dataMsg* msg = (dataMsg*) buff;
    uint8_t* data = msg->data+msg->tagSize;
    nvdla::IRuntime::NvDlaTensor tensor;
    ctx->rt->getInputTensorDesc(0, &tensor);
    e =  copyImageToInputTensor(ctx->rt, ctx->input,msg->dataSize, data, &tensor);
    sendAckMsg(ctx->connection, INPUT_DATA);
    return e;
}
void performInfer(NvdlaContext *ctx, uint8_t* buff){
    std::cout<<"received infer"<<std::endl;
    prepareOutputTensor(ctx->rt,ctx->output);
    inferMsg* msg = (inferMsg*) buff;
    if(ctx->output->handle && ctx->output->data){
        std::cout<<"submitting task"<<std::endl;
        ctx->rt->submit();  
        uint8_t* outBuf = new uint8_t[msg->outputSize];
        prepareOutputBuffer(ctx->rt, ctx->output, msg->outputSize,  outBuf);
        sendData(ctx->connection,OUTPUT_DATA,outBuf, msg->outputSize);
        delete[] outBuf;
        // sendAckMsg(connection, INFER);
    }
}


void pollConnection(NvdlaContext* ctx){
    // std::cout<<"ctx "<<(void*)ctx<<" "<<(void*)ctx->connection<<std::endl;
    // while (connection){
        uint32_t closeConCount =0;
        ctx->connection->lock();
        while(pollWrapper(ctx->connection) > 0){ //we have a message!
            uint8_t * buff = NULL;
            if ( pollRecWrapper(ctx->connection,&buff) >= (int64_t) sizeof(msgHeader) ){
                // printMsgHeader(buff);
                msgHeader *header = (msgHeader*)buff;
                switch(header->type) {
                    case INIT: {
                        init_nvdla(ctx);
                        break;
                    }
                    case NUM_DEVICES: {
                        send_device_info(ctx,NUM_DEVICES);
                        break;
                    }
                    case MAX_DEVICES: {
                        send_device_info(ctx,MAX_DEVICES);
                        break;
                    }
                    case LOAD_BINARY: {
                        load_binary_image(ctx,buff);
                        break;
                    }
                    case INPUT_DATA: {
                        load_input_data(ctx,buff);
                        break;
                    }
                    case INFER: {
                        performInfer(ctx, buff);
                        break;
                    }
                    case CLOSE_SERVER_MSG: {
                        // shutDownServer();
                        break;
                    }
                    case PING_MSG: {
                        std::cout<<"received ping sending ack"<<std::endl;
                        sendAckMsg(ctx->connection,PING_MSG);
                        break;
                    }
                    case CLOSE_CON_MSG: {

                    }
                    default: {//Close connection on incorrect message type
                        ctx->connection->closeSocket();
                    }
                }
                if (buff){
                    delete[] buff;
                }
            }
            else { //We were not able to read a message header... Must be an error... close socket
                printf("Recv failed for client %s\n", ctx->connection->addr().c_str());
                ctx->connection->closeSocket();
            }

        }
        if (ctx->connection->unlock()) { //this means we closed the connection
            delete ctx;
            ctx = NULL;
        }
        else{
            workPool.addTask([ctx] { pollConnection(ctx); });
        }
    // }
}

int main(){

    uint16_t port=6667;
    struct sockaddr_in address;
    memset((char*)&address,0,sizeof(address));
    address.sin_family = AF_INET; 
    address.sin_addr.s_addr = htons(INADDR_ANY); 
    address.sin_port = htons( port );


    int sock_fd = socket(AF_INET,SOCK_STREAM,0);
    if (sock_fd < 0 ){
        std::cerr<<ERROR_TAG<<"failed to get server socket"<<std::endl;
        exit(1);
    }
    int opt=1;
    if(setsockopt(sock_fd,SOL_SOCKET,SO_REUSEADDR,&opt,sizeof(opt)) < 0){
        std::cerr<<ERROR_TAG<<"failed to set sockoptions"<<std::endl;
        exit(1);
    }
    if (bind(sock_fd, (struct sockaddr *)&address, sizeof(address))<0) { 
        std::cerr<<ERROR_TAG<<"failed bind socket"<<std::endl;
        exit(1);
    } 

    if (listen(sock_fd,3) < 0){
        std::cerr<<ERROR_TAG<<"failed listen socket"<<std::endl;
        exit(1);
    }

    socklen_t cli_len;
    sockaddr_in cli_addr;
    while(true){
        int new_sock = accept(sock_fd, (struct sockaddr *)&cli_addr,&cli_len);
        if (new_sock < 0){
            std::cerr<<ERROR_TAG<<"failed accept socket "<<std::strerror(errno)<<" "<<inet_ntoa(cli_addr.sin_addr)<<":"<<ntohs(cli_addr.sin_port)<<std::endl;
            exit(1);
        }
        workPool.addThreadWithTask([=] { 
            Connection* con = new Connection(new_sock,inet_ntoa(cli_addr.sin_addr),cli_addr.sin_port);
            NvdlaContext* ctx = new NvdlaContext(con);
            pollConnection(ctx);
        });

    }


    
}

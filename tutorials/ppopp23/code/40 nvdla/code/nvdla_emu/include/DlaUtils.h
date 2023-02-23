/*
 * Copyright (c) 2017-2019, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef DLA_UTILS_H
#define DLA_UTILS_H

#include <stdbool.h>
#include <sstream>
#include <iostream>
#include "nvdla/IRuntime.h"
#include "DlaImage.h"

#include "dlaerror.h"
#include "dlatypes.h"
#include "Connection.h"



class Connection;


struct NvdlaData{
    void* handle;
    uint64_t size;
    uint8_t* data;
    NvDlaImage* image;

    NvdlaData(){
        handle =NULL;
        size = 0;
        data = NULL;
        image = NULL;
    }

    ~NvdlaData(){
        delete image;
    }
};

struct NvdlaContext{
    Connection* connection;
    std::string network;
    nvdla::IRuntime *rt;
    NvdlaData *input;
    NvdlaData *output;
    NvdlaContext(Connection* con){
        connection = con;
        network="";
        rt = NULL;
        input = new NvdlaData();
        output = new NvdlaData();
        // std::cout<<"NvdlaContext: "<<(void*)this<<" "<<(void*)connection<<std::endl;
    }

    ~NvdlaContext(){
        if(rt){
            rt->unload();
            rt->freeSystemMemory(input->handle,input->size);
            rt->freeSystemMemory(output->handle,output->size);
            if (input){
                delete input;
            }
            if (output){
                delete output;
            }
            rt->stopEMU();
            std::cout<<"destroying runtime"<<std::endl;
            nvdla::destroyRuntime(rt);
        }
        
        delete connection;
    }
};

NvDlaError copyImageToInputTensor(nvdla::IRuntime *rt, NvdlaData* input, uint64_t inputSize, uint8_t* inputBuf, nvdla::IRuntime::NvDlaTensor *tensorDesc);
NvDlaError prepareOutputTensor(nvdla::IRuntime* rt, NvdlaData* output);
NvDlaError prepareOutputBuffer(nvdla::IRuntime* rt, NvdlaData* output, uint64_t size, uint8_t* outBuf);
NvDlaError Tensor2DIMG(const nvdla::IRuntime::NvDlaTensor* pTDesc, NvDlaImage* image);
NvDlaError DlaBuffer2DIMG(void** pBuffer, NvDlaImage* image);
NvDlaError DIMG2DlaBuffer(const NvDlaImage* image, void** pBuffer);

#endif // NVDLA_UTILS_DLAIMAGE_H

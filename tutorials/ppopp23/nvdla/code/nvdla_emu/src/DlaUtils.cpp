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

#include <math.h>

#include "DlaUtils.h"
#include "dlatypes.h"
#include "dlaerror.h"

#include "ErrorMacros.h"
#include "nvdla_os_inf.h"
#include "half.h"

#include <stdio.h> // snprintf


NvDlaError copyImageToInputTensor(nvdla::IRuntime *rt, NvdlaData* input, uint64_t inputSize, uint8_t* inputBuf, nvdla::IRuntime::NvDlaTensor *tensorDesc){
    NvDlaError e = NvDlaSuccess;

    uint64_t num_elems =tensorDesc->dims.n*tensorDesc->dims.c*tensorDesc->dims.h*tensorDesc->dims.w;
    uint64_t expectedElemSize = inputSize/num_elems;
    std::cout<<"expectedElemSize "<<expectedElemSize <<" "<<inputSize<<"/"<<num_elems<<std::endl;
    // for (int i =0; i<num_elems; i++){
		
	// 	printf("%c%c"," .:-=+*#%@"[(int)inputBuf[i] / 26 ],((i + 1) % 28) ? '\0' : '\n');
	// }
    NvDlaImage* R8Image = new NvDlaImage();
    NvDlaImage* tensorImage = new NvDlaImage();
    
    R8Image->m_meta.surfaceFormat = NvDlaImage::T_R8;
    R8Image->m_meta.width = tensorDesc->dims.w;
    R8Image->m_meta.height =  tensorDesc->dims.h;
    R8Image->m_meta.channel = 1;

    R8Image->m_meta.lineStride = tensorDesc->stride[1];
    R8Image->m_meta.surfaceStride = tensorDesc->stride[2];;
    R8Image->m_meta.size = tensorDesc->bufferSize;
    R8Image->m_pData = NvDlaAlloc(R8Image->m_meta.size);
    void* buf = R8Image->m_pData;
    memset(buf, 0, R8Image->m_meta.size);
    for (NvU32 y=0; y<R8Image->m_meta.height; y++){
        uint8_t* linebuf = reinterpret_cast<uint8_t*>(buf) + (y * R8Image->m_meta.lineStride);
        for (NvU32 x=0; x<R8Image->m_meta.width; x++){
            // std::cout<<"("<<(uint64_t) reinterpret_cast<uint8_t*>(inputBuf)[y*R8Image->m_meta.width + x]<<" ";
            // std::cout<<(uint64_t) reinterpret_cast<uint16_t*>(inputBuf)[y*R8Image->m_meta.width + x]<<" ";
            // std::cout<<(uint64_t) reinterpret_cast<uint32_t*>(inputBuf)[y*R8Image->m_meta.width + x]<<" ";
            // std::cout<<(uint64_t) reinterpret_cast<uint64_t*>(inputBuf)[y*R8Image->m_meta.width + x]<<std::endl;
            switch(expectedElemSize){
                case 1:{
                    linebuf[x] = inputBuf[y*R8Image->m_meta.width + x];
                    break;
                }
                case 2:{
                    linebuf[x] = (uint8_t) reinterpret_cast<uint16_t*>(inputBuf)[y*R8Image->m_meta.width + x];
                    break;
                }
                case 4:{
                    if (tensorDesc->dataType == NVDLA_DATA_TYPE_HALF || tensorDesc->dataType == NVDLA_DATA_TYPE_FLOAT){
                        linebuf[x] = (uint8_t) reinterpret_cast<float*>(inputBuf)[y*R8Image->m_meta.width + x];
                    }
                    else{
                        linebuf[x] = (uint8_t) reinterpret_cast<uint32_t*>(inputBuf)[y*R8Image->m_meta.width + x];
                    }
                    break;
                }
                case 8:{
                    linebuf[x] = (uint8_t) reinterpret_cast<uint64_t*>(inputBuf)[y*R8Image->m_meta.width + x];
                    break;
                }
                default:{
                    std::cerr<<"unexpected input elem size: "<<expectedElemSize<<std::endl;
                    break;
                }
            }
            // printf("%c%c"," .:-=+*#%@"[(int)linebuf[x] / 26 ],((y*R8Image->m_meta.width + x + 1) % 28) ? '\0' : '\n');
        }

        // memcpy(linebuf,inputBuf+y*R8Image->m_meta.width,R8Image->m_meta.width);
    }
    R8Image->printInfo();
    R8Image->printBuffer(true);

    tensorImage->m_meta.width = tensorDesc->dims.w;
    tensorImage->m_meta.height =  tensorDesc->dims.h;
    tensorImage->m_meta.channel = tensorDesc->dims.c;

    switch(tensorDesc->pixelFormat){
    case NVDLA_PIXEL_FORMAT_R8:
        std::cout<<"NVDLA_PIXEL_FORMAT_R8"<<std::endl;
        tensorImage->m_meta.surfaceFormat = NvDlaImage::T_R8;
        tensorImage->m_meta.lineStride    = tensorDesc->stride[1];
        tensorImage->m_meta.surfaceStride = 0;
        tensorImage->m_meta.size          = tensorDesc->bufferSize;
        break;

    case NVDLA_PIXEL_FORMAT_FEATURE:
        std::cout<<"NVDLA_PIXEL_FORMAT_FEATURE"<<std::endl;
        if (tensorDesc->dataType == NVDLA_DATA_TYPE_HALF){
            std::cout<<"NVDLA_DATA_TYPE_HALF"<<std::endl;
            tensorImage->m_meta.surfaceFormat = NvDlaImage::D_F16_CxHWx_x16_F;
        }
        else if (tensorDesc->dataType == NVDLA_DATA_TYPE_INT8){
            std::cout<<"NVDLA_DATA_TYPE_INT8"<<std::endl;
            tensorImage->m_meta.surfaceFormat = NvDlaImage::D_F8_CxHWx_x32_I;
        }
        else{
            ORIGINATE_ERROR(NvDlaError_NotSupported, "Unsupported (pixel, data) combination: (%u, %u)", tensorDesc->pixelFormat, tensorDesc->dataType);
        }

        tensorImage->m_meta.lineStride    = tensorDesc->stride[1];
        tensorImage->m_meta.surfaceStride = tensorDesc->stride[2];
        tensorImage->m_meta.size          = tensorDesc->bufferSize;
        break;

    case NVDLA_PIXEL_FORMAT_FEATURE_X8:
        std::cout<<"NVDLA_PIXEL_FORMAT_FEATURE_X8"<<std::endl;
        if (tensorDesc->dataType == NVDLA_DATA_TYPE_INT8){
            std::cout<<"NVDLA_DATA_TYPE_INT8"<<std::endl;
            tensorImage->m_meta.surfaceFormat = NvDlaImage::D_F8_CxHWx_x8_I;
        }
        else{
            ORIGINATE_ERROR(NvDlaError_NotSupported, "Unsupported (pixel, data) combination: (%u, %u)", tensorDesc->pixelFormat, tensorDesc->dataType);
        }

        tensorImage->m_meta.lineStride    = tensorDesc->stride[1];
        tensorImage->m_meta.surfaceStride = tensorDesc->stride[2];
        tensorImage->m_meta.size          = tensorDesc->bufferSize;
        break;

    case NVDLA_PIXEL_FORMAT_A16B16G16R16_F:
        std::cout<<"NVDLA_PIXEL_FORMAT_A16B16G16R16_F"<<std::endl;
        tensorImage->m_meta.surfaceFormat = NvDlaImage::T_A16B16G16R16_F;
        tensorImage->m_meta.lineStride    = tensorDesc->stride[1];
        tensorImage->m_meta.surfaceStride = 0;
        tensorImage->m_meta.size          = tensorDesc->bufferSize;
        break;

    case NVDLA_PIXEL_FORMAT_A8B8G8R8:
        std::cout<<"NVDLA_PIXEL_FORMAT_A8B8G8R8"<<std::endl;
        tensorImage->m_meta.surfaceFormat = NvDlaImage::T_A8B8G8R8;
        tensorImage->m_meta.lineStride    = tensorDesc->stride[1];
        tensorImage->m_meta.surfaceStride = 0;
        tensorImage->m_meta.size          = tensorDesc->bufferSize;
        break;
    default:
        ORIGINATE_ERROR(NvDlaError_NotSupported, "Unsupported pixel format: %u", tensorDesc->pixelFormat);
    }

    tensorImage->m_pData = NvDlaAlloc(tensorImage->m_meta.size);
    memset(tensorImage->m_pData, 0, tensorImage->m_meta.size);
    uint8_t* ibuf = static_cast<uint8_t*>(R8Image->m_pData);
    uint8_t* obuf = static_cast<uint8_t*>(tensorImage->m_pData);

    for (NvU32 y=0; y < R8Image->m_meta.height; y++)
    {
        for (NvU32 x=0; x < R8Image->m_meta.width; x++)
        {
            for (NvU32 z=0; z < R8Image->m_meta.channel; z++)
            {
                NvS32 ioffset = R8Image->getAddrOffset(x, y, z);
                NvS32 ooffset = tensorImage->getAddrOffset(x, y, z);

                if (ioffset < 0)
                    ORIGINATE_ERROR(NvDlaError_BadParameter);
                if (ooffset < 0)
                    ORIGINATE_ERROR(NvDlaError_BadParameter);

                NvU8* inp = ibuf + ioffset;

                if (tensorDesc->dataType == NVDLA_DATA_TYPE_HALF)
                {
                    half_float::half* outp = reinterpret_cast<half_float::half*>(obuf + ooffset);
                    *outp = static_cast<half_float::half>(float(*inp));
                }
                else if (tensorDesc->dataType == NVDLA_DATA_TYPE_INT8)
                {
                    char* outp = reinterpret_cast<char*>(obuf + ooffset);
                    //*outp = char(*inp); // no normalization happens
                    // compress the image from [0-255] to [0-127]
                    *outp = static_cast<NvS8>(std::floor((*inp * 127.0/255.0) + 0.5f));
                }
            }
        }
    }

    if (input->handle == NULL){
        rt->allocateSystemMemory(&input->handle, tensorDesc->bufferSize, (void**)&input->data);
        if (!rt->bindInputTensor(0, input->handle)){
            ORIGINATE_ERROR_FAIL(NvDlaError_BadParameter, "runtime->bindInputTensor() failed");
        }
        input->size = tensorDesc->bufferSize;
    }
    memcpy(input->data,(uint8_t*)tensorImage->m_pData,tensorImage->m_meta.size);

    


fail:
    if (R8Image != NULL && R8Image->m_pData != NULL)
        NvDlaFree(R8Image->m_pData);
    delete R8Image;

    return e;
}

 NvDlaError prepareOutputTensor(nvdla::IRuntime* rt, NvdlaData* output){
    NvDlaError e = NvDlaSuccess;
    nvdla::IRuntime::NvDlaTensor tensorDesc;
    PROPAGATE_ERROR_FAIL(rt->getOutputTensorDesc(0, &tensorDesc));
    output->image = new NvDlaImage();
    if (output->handle == NULL){
        rt->allocateSystemMemory(&(output->handle), tensorDesc.bufferSize, (void**)&(output->data));
        rt->bindOutputTensor(0, output->handle);
        output->size = tensorDesc.bufferSize;
    }
    

    PROPAGATE_ERROR_FAIL(Tensor2DIMG(&tensorDesc, output->image));
    PROPAGATE_ERROR_FAIL(DIMG2DlaBuffer(output->image, (void**)&output->data));
    

fail:
    return e;
}

NvDlaError prepareOutputBuffer(nvdla::IRuntime* rt, NvdlaData* output, uint64_t size, uint8_t* outBuf){
    NvDlaError e = NvDlaSuccess;
    std::stringstream sstream;
    nvdla::IRuntime::NvDlaTensor tensorDesc;
    uint64_t num_elems;
    uint64_t expectedElemSize;
    PROPAGATE_ERROR_FAIL(DlaBuffer2DIMG((void**)&output->data, output->image));
    PROPAGATE_ERROR_FAIL(rt->getOutputTensorDesc(0, &tensorDesc));
    num_elems = tensorDesc.dims.c*tensorDesc.dims.h*tensorDesc.dims.w*tensorDesc.dims.n;
    expectedElemSize = size/num_elems;
    output->image->packData(outBuf,expectedElemSize);
fail:
    return e;
}

NvDlaError Tensor2DIMG(const nvdla::IRuntime::NvDlaTensor* tensor, NvDlaImage* image)
{
    if (!tensor || !image)
        ORIGINATE_ERROR(NvDlaError_BadParameter);

    // Fill surfaceFormat
    NvDlaImage::PixelFormat surfaceFormat;
    switch (tensor->pixelFormat)
    {
        case TENSOR_PIXEL_FORMAT_FEATURE:
        {
            if (tensor->dataType == TENSOR_DATA_TYPE_HALF) {
                surfaceFormat = NvDlaImage::D_F16_CxHWx_x16_F;
            } else if (tensor->dataType == TENSOR_DATA_TYPE_INT16) {
                surfaceFormat = NvDlaImage::D_F16_CxHWx_x16_I;
            } else if (tensor->dataType == TENSOR_DATA_TYPE_INT8) {
                surfaceFormat = NvDlaImage::D_F8_CxHWx_x32_I;
            } else {
                ORIGINATE_ERROR(NvDlaError_NotSupported);
            }
        }
        break;
        case TENSOR_PIXEL_FORMAT_FEATURE_X8:
        {
            if (tensor->dataType == TENSOR_DATA_TYPE_HALF) {
                surfaceFormat = NvDlaImage::D_F16_CxHWx_x16_F;
            } else if (tensor->dataType == TENSOR_DATA_TYPE_INT16) {
                surfaceFormat = NvDlaImage::D_F16_CxHWx_x16_I;
            } else if (tensor->dataType == TENSOR_DATA_TYPE_INT8) {
                surfaceFormat = NvDlaImage::D_F8_CxHWx_x8_I;
            } else {
                ORIGINATE_ERROR(NvDlaError_NotSupported);
            }
        }
        break;

        default:
            REPORT_ERROR(NvDlaError_NotSupported, "Unexpected surface format %u, defaulting to D_F16_CxHWx_x16_F", tensor->pixelFormat);
            surfaceFormat = NvDlaImage::D_F16_CxHWx_x16_F;
            break;
    }

    image->m_meta.surfaceFormat = surfaceFormat;
    image->m_meta.width = tensor->dims.w;
    image->m_meta.height = tensor->dims.h;
    image->m_meta.channel = tensor->dims.c;

    image->m_meta.lineStride = tensor->stride[1];
    image->m_meta.surfaceStride = tensor->stride[2];
    image->m_meta.size = tensor->bufferSize;

    if ( 0 ) {
        NvDlaDebugPrintf("tensor2dimg: image meta: format=%u width=%u height=%u lineStride=%u surfaceStride=%u size=%u\n",
                        image->m_meta.surfaceFormat,
                        image->m_meta.width,  image->m_meta.height,
                        image->m_meta.lineStride, image->m_meta.surfaceStride,
                        (NvU32)image->m_meta.size);
    }

    // Allocate the buffer
    image->m_pData = NvDlaAlloc(image->m_meta.size);
    if (!image->m_pData)
        ORIGINATE_ERROR(NvDlaError_InsufficientMemory);

    // Clear the data
    memset(image->m_pData, 0, image->m_meta.size);

    return NvDlaSuccess;
}

NvDlaError DlaBuffer2DIMG(void** pBuffer, NvDlaImage* image){
    if (!(*pBuffer) || !image)
        ORIGINATE_ERROR(NvDlaError_BadParameter);

    memcpy(image->m_pData, *pBuffer, image->m_meta.size);

    return NvDlaSuccess;
}

NvDlaError DIMG2DlaBuffer(const NvDlaImage* image, void** pBuffer){
    if (!image || !(*pBuffer))
        ORIGINATE_ERROR(NvDlaError_BadParameter);

    memcpy(*pBuffer, image->m_pData, image->m_meta.size);

    return NvDlaSuccess;
}

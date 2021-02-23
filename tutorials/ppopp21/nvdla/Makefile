LIBS=-lmcl -lOpenCL -lm -lpthread -lrt
ifeq ($(CC),)
    CC := gcc
endif

NVDLA_EMU_BASE=/usr/local/nvdla
NVDLA_EMU_COMPILER=$(NVDLA_EMU_BASE)/nvdla_compiler 

ifeq ($(caffe),)
	caffe=mnist/mnist.caffemodel
endif

ifeq ($(prototxt),)
	prototxt=mnist/mnist.prototxt
endif


nvdla_emu_bin:
	@[ ! -f "${NVDLA_EMU_COMPILER}" ] || ( echo "NVDLA_EMU_COMPILER does not exist"; ls -altr ${NVDLA_EMU_COMPILER}; exit 1)
	@[ "${caffe}" ] || (echo "nvdla caffe image not supplied e.g. caffe=model.caffe "; exit 1)
	@[ "${prototxt}" ] || (echo "nvdla prototxt not supplied e.g. prototxt=model.prototxt"; exit 1)
	LD_LIBRARY_PATH=$(NVDLA_EMU_BASE) $(NVDLA_EMU_COMPILER) --prototxt ${prototxt} --caffemodel ${caffe}
	mv fast-math.nvdla mnist/mnist.nvdla
	rm output.protobug
	rm -r wisdom.dir

nvdla_test: nvdla.c
	$(CC) -g -O3  nvdla.c -o nvdla_test $(LIBS)
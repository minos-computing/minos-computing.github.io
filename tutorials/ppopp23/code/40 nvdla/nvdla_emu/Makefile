
TOOLCHAIN_PREFIX ?=



ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

ifeq ($(CXX),)
    CXX := g++
endif




NVDLA_SW_UMD_INC=/usr/local/nvdla/sw/umd/core/include
NVDLA_SW_UMD_LIB=/usr/local/nvdla/
NVDLA_SW_EXT_INC=/usr/local/nvdla/sw/umd/external/include

COMMON := Connection Message
COMMON_CPP := $(addsuffix .cpp, $(COMMON))
COMMON_H := $(addsuffix .h, $(COMMON))

SERVER_OBJS := $(addprefix src/, $(addsuffix _server.o,$(basename $(COMMON_CPP))) DlaUtils_server.o DlaImage_server.o WorkPool_server.o)
CLIENT_OBJS := $(addprefix src/, $(addsuffix _client.o,$(basename $(COMMON_CPP))))


CXX_FLAGS:= -std=c++11 -DNVDLA_UTILS_ERROR_TAG="\"DLAEMU\"" -g -O3

all: nvdla_emu_server client

src/%_server.o: src/%.cpp include/%.h
	@[ "${TOOLCHAIN_PREFIX}" ] || ( echo "TOOLCHAIN_PREFIX is not set"; exit 1 )
	$(TOOLCHAIN_PREFIX)$(CXX) $(CXX_FLAGS)  -c $< -o $@ -Iinclude -I$(NVDLA_SW_UMD_INC) -I$(NVDLA_SW_EXT_INC) 

%_client.o: %.cpp 
	$(CXX) $(CXX_FLAGS) -fpic -c $< -o $@ -Iinclude

nvdla_emu_server: src/server.cpp $(SERVER_OBJS) 
	@[ "${TOOLCHAIN_PREFIX}" ] || ( echo "TOOLCHAIN_PREFIX is not set"; exit 1 )
	$(TOOLCHAIN_PREFIX)$(CXX) $(CXX_FLAGS) src/server.cpp -o nvdla_emu_server $(SERVER_OBJS) -Iinclude -I$(NVDLA_SW_UMD_INC) -I$(NVDLA_SW_EXT_INC) -L$(NVDLA_SW_UMD_LIB) -lnvdla_runtime -lpthread 

client: src/client_test.cpp  $(CLIENT_OBJS)
	$(CXX) $(CXX_FLAGS) -shared -o libnvdla_emu.so $(CLIENT_OBJS)
	$(CXX) $(CXX_FLAGS) src/client_test.cpp -o client -Iinclude  -L$(CURDIR) -lnvdla_emu

install: nvdla_emu_server client

	install -d $(DESTDIR)$(PREFIX)/bin/
	install -m 755 nvdla_emu_server $(DESTDIR)$(PREFIX)/bin/
	install -d $(DESTDIR)$(PREFIX)/include/nvdla_emu
	install -m 644 include/Message.h $(DESTDIR)$(PREFIX)/include/nvdla_emu
	install -m 644 include/Connection.h $(DESTDIR)$(PREFIX)/include/nvdla_emu
	install -d $(DESTDIR)$(PREFIX)/lib/
	install -m 755 libnvdla_emu.so $(DESTDIR)$(PREFIX)/lib/

uninstall:
	rm -r $(DESTDIR)$(PREFIX)/bin/nvdla_emu_server || true
	rm -r $(DESTDIR)$(PREFIX)/include/nvdla_emu || true
	rm -r   $(DESTDIR)$(PREFIX)/lib/libnvdla_emu.so || true



clean:
	rm src/*.o || true
	rm nvdla_emu_server || true
	rm libnvdla_emu.so || true
	rm client || true









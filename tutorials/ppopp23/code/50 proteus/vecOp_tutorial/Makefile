CC=gcc-8

CFLAGS=-I. -I/usr/include
LIBS=-L/usr/lib -lmcl -lOpenCL -pthread -lm -lrt

DEPS = utils.h
OBJ = utils.o vecOp.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

vecOp: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)  

.PHONY: clean

clean:
	rm *.o vecOp

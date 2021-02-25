# MCL + NVDLA Emulator Tutorial

## Welcome
This portion of the tutorial provides a working example of performing Digit recognition using the NVDLA Virtual Platform (http://nvdla.org/vp.html) to emulate an open source NVDLA.
We have developed a docker containing both an NVDLA installation and all the appropraite OpenCL, Pocl, and MCL libraries.
Within this example we will cover the following topics:
* Highlevel overview of NVDLA QEMU configuration file
* Launching of NVDLA QEMU image + mcl net attached device server interface
* Compiling a caffe model into a nvdla binary
* Example MCL application to perform digit recognition
* Launching MCL scheduler using NVDLA devices
* Running our test application on the NVDLA emulators

## TL;DR
```
docker pull minoscomputing/ppopp21
docker run --name ppopp21 -it minoscomputing/ppopp21 /bin/bash
cd nvdla
./start_nvdla_emulator.sh 2

#start new terminal
docker exec -it  ppopp21 /bin/bash
cd nvdla
make #compiles both nvdla model and mcl test application
./launch_mcl_sched.sh
./run_nvdla_test.sh
```


## Downloading and starting the Docker image
```
docker pull mcl/ppopp21
docker run --name ppopp21 -it mcl/ppopp21 /bin/bash
cd nvdla
```

## Create Qemu configs and Launch NVDLA emulators
We have provided a script 'start_nvdla_emulator.sh' that handles correctly setting QEMU configs, starting the images, initializing the NVDLA device driver, and starting the net attached server.
Please see the tutorial slides for further details.
This script takes an option argument specifying the number of virtual NVDLAs to launch, for the tutorial we will use 2.
```
./start_nvdla_emulator.sh 2
```
We are going to leave this terminal running for now, eventually we will see some new output once inference tasks get assigned to the NVDLAs.
## Compile NVDLA Mnist Digit Recognition Model
Start a new server and attach to the currently running docker image
```
docker exec -it  ppopp21 /bin/bash
cd /nvdla
```
First we compile the caffe model (located in the mnist folder) into a binary that can be used the the NVDLA emulators.
The nvdla_compiler is part of the open source NVDLA SW platform and currently only parses caffe models.
```
make prototxt=mnist/mnist.prototxt caffe=mninst/mnist.caffemodel nvdla_emu_bin # explicitly providing model paths for illustrative purposes
```
the resultant binary mnist.nvdla will be located in the mnist folder.

## Compile MCL Nvdla example
The source code for our MCL NVDLA example is contained in nvdla.c (see tutorial slides for additional discussion).
```
make nvdla_test
```
the resultant binary is called 'nvdla_test'

## Launch the MCL scheduler
Recall that MCL utilizes a scheduler process on each 'node' which manages the various devices present.
Before executing an MCL application we must ensure the scheduler is running.
```
./launch_mcl_sched.sh
```
This script ensures MCL_SCHED is launched with the appropriate POCL environment variables so the the scheduler can discover the NVDLA emulators.
e.g. POCL_DEVICES='nvdlaem' POCL_NVDLAEMU0_PARAMETERS='127.0.0.1:6000'  (the actual values will be set according to number of emulator instances launched)
The scheduler gets launched in the background (so no need for a new terminal)

## Execute the test
Finally we are ready to actually perform our MCL + NVDLA inference test
```
./run_nvdla_test.sh
```
Again this a wrapper script which appropriately sets the POCL environment variables to use the NVDLA emulators.

The test reads in all the pgm image files in the mnist folder and then performs synchronous inferences for two digits (chosen at random)  followed by asynchronous inferences for two additional digits (also chosen at random). Sometimes the same digit is chosen more than once, rerunning the test should produce different digits. The model itself does not have perfect accuracy (for example it mis identifies the provided 0 image as a 5). Given we are using two NVDLA instances in this example, one would expect the execution time of the asynchronous test to be half that of the synchronous test (i.e. each instance is performing an instance simultaneously), but due to the docker+QEMU architecture execution appears to be serialized.

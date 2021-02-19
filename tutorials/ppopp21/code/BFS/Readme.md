# BFS Tutorial

## Welcome
This section of the tutorial is intended to experiment with a more complex benchmark. There are 2 applications in this folder: BFS, and BFS_modified. If you are following along in the tutorial, we began by running BFS, then make the modifications necessary to use resident memory for breadth first search. To build the modified version, uncomment that section from the CMakeLists.txt file. This section of the tutorial was focused on multiple GPUs. If you don't have access to multiple GPUs, you may run these benchmarks on your CPU or single GPU, just specify a smaller size using the -s, -r and -n options. For a full description of the op

## Building the Code
First intall MCL following the instructions here: #TODO: Include link to public MCL
This part of the tutorial follows a typical cmake build pattern. From this folder, run the commands:
```
mkdir b && cd b
cmake ..
make
```

## To Run the Code
First you must run the MCL scheduler. In the tutoria, we began with the basic round robin scheduler and used the BFS benchmark to demonstrate the usefullness of the locality aware schedulers. For best performance, run:
```
mcl_sched -s fffs -p hybrid
```
You may adjust the hyperparameters of the locality aware scheduler using the environment variables `MCL_SCHED_MAX_ATTEMPTS` and `MCL_SCHED_COPY_FACTOR`. For a full explanation of the locality aware scheduler see the [tutorial](https://minos-computing.github.io/tutorials/ppopp21/ppopp21.html) or the [paper](https://ieeexplore.ieee.org/document/9307939). Once the scheduler is running, you may run the code with:
```
./bfs [options]
```
For a full list of options, run `./bfs -h`. Beware that the larger graphs (`-s 3` and `-s 4` and well as the more saerches `-r 4` take lots of time on CPU or with a non locality-aware scheduler.)

## Attribution and Thanks
The BFS benchmark as well as other supporting files comes from the [SHOC benchmark suite] (https://github.com/vetter/shoc). The files are moved here for convenience during the tutorial.  The paper describing the SHOC benchmark suite is below:
Anthony Danalis, Gabriel Marin, Collin McCurdy, Jeremy S. Meredith, Philip C. Roth, Kyle Spafford, Vinod Tipparaju, and Jeffrey S. Vetter. 2010. The Scalable Heterogeneous Computing (SHOC) benchmark suite. In Proceedings of the 3rd Workshop on General-Purpose Computation on Graphics Processing Units (GPGPU-3). Association for Computing Machinery, New York, NY, USA, 63–74. DOI:https://doi.org/10.1145/1735688.1735702 
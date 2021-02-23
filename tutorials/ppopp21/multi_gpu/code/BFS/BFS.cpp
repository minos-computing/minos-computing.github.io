#include<fstream>
#include<iostream>
#include<limits.h>
#include<math.h>
#include<set>
#include<string.h>
#include<sys/time.h>
#include<time.h>

#include "Graph.h"
#include "OptionParser.h"
#include "ResultDatabase.h"
#include "utils.h"

// ****************************************************************************
// Function: addBenchmarkSpecOptions
//
// Purpose:
//   Add benchmark specific options parsing
//
// Arguments:
//   op: the options parser / parameter database
//
// Returns:  nothing
//
// Programmer: Alok Kamatar
// Creation: June 16, 2020
//
// Modifications:
//
// ****************************************************************************
void addBenchmarkSpecOptions(OptionParser &op)
{
    op.addOption("graph_file", OPT_STRING, "random", "name of graph file");
    op.addOption("degree", OPT_INT, "2", "average degree of nodes");
    op.addOption("dump-pl", OPT_BOOL, "false",
            "enable dump of path lengths to file");
    op.addOption("source_vertex", OPT_INT, "0",
            "vertex to start the traversal from");
    op.addOption("verify", OPT_BOOL, "false",
            "vertex to start the traversal from");
}

// ****************************************************************************
// Function: GetWorkSize
//
// Purpose:
//   Get the kernel configuration so that we have one thread mapped to one
//   vertex in the frontier.
//
// Arguments:
//   gws: global work size
//   lws: local work size
//   maxThreadsPerCore: the max work group size for specified device
//   numVerts: number of vertices in the graph
//
// Returns:  nothing
//
// Programmer: Aditya Sarwade
// Creation: June 16, 2011
//
// Modifications:
//   Alok Kamatar (06/16/2020): Taken from OpenCL and copied to MCL
// ****************************************************************************
void GetWorkSize(size_t *gws, size_t *lws, size_t maxThreadsPerCore,
                 cl_uint numVerts)
{
    int temp;
    gws[0]=(size_t)numVerts;
    temp=(int)ceil((float)gws[0]/(float)maxThreadsPerCore);
    gws[0]=temp*maxThreadsPerCore;
    lws[0]=maxThreadsPerCore;
}

// ****************************************************************************
// Function: verify_results
//
// Purpose:
//  Verify BFS results by comparing the output path lengths from cpu and gpu
//  traversals
//
// Arguments:
//   cpu_cost: path lengths calculated on cpu
//   gpu_cost: path lengths calculated on gpu
//   numVerts: number of vertices in the given graph
//   out_path_lengths: specify if path lengths should be dumped to files
//
// Returns:  nothing
//
// Programmer: Aditya Sarwade
// Creation: June 16, 2011
//
// Modifications:
//   Alok Kamatar (06/16/2020): Taken from OpenCL and copied to MCL
// ****************************************************************************
cl_uint verify_results(cl_uint *cpu_cost, cl_uint *gpu_cost, cl_uint numVerts,
                       bool out_path_lengths)
{
    cl_uint unmatched_nodes=0;

    //check cpu against gpu results
    for (int i=0;i<numVerts;i++)
    {
        if (gpu_cost[i]!=cpu_cost[i])
        {
            unmatched_nodes++;
        }
    }

    //if user wants to write path lengths to file
    if (out_path_lengths)
    {
        std::ofstream bfs_out_cpu("bfs_out_cpu.txt");
        std::ofstream bfs_out_gpu("bfs_out_ocl.txt");
        for (int i=0;i<numVerts;i++)
        {
            if (cpu_cost[i]!=UINT_MAX)
                bfs_out_cpu<<cpu_cost[i]<<"\n";
            else
                bfs_out_cpu<<"0\n";

            if (gpu_cost[i]!=UINT_MAX)
                bfs_out_gpu<<gpu_cost[i]<<"\n";
            else
                bfs_out_gpu<<"0\n";
        }
        bfs_out_cpu.close();
        bfs_out_gpu.close();
    }
    return unmatched_nodes;
}

// ******************************************************************
// Function: RunKernel1
//
// Purpose:
//   Runs the bfs_iiit t kernel and returns the runtime
//
// Arguments:
//
// Returns: Runtime in seconds
// Author: Alok Kamatar
// Created June 16, 2020
//
// ******************************************************************
mcl_handle* RunKernel1(char* kernel_path, uint32_t* frontier, uint32_t numVerts,
        uint32_t* edge_offsets, uint32_t* edge_list, uint32_t adj_list_length, uint32_t W_SZ, uint32_t CHUNK_SZ, 
        uint32_t iters, int* flag, uint64_t* pes, uint64_t* lsize) {
    int err;
    mcl_handle* hdl = mcl_task_create();
    if(!hdl) {
        MCL_CHECK_ERROR(1);
    }
    err = mcl_task_set_kernel(hdl, kernel_path, "BFS_kernel_warp", 8, "", 0x0);
    MCL_CHECK_ERROR(err);
    
    err = mcl_task_set_arg(
        hdl, 0, (void*)frontier, sizeof(uint32_t) * numVerts,
        MCL_ARG_BUFFER | MCL_ARG_INPUT | MCL_ARG_OUTPUT
    );
    
    MCL_CHECK_ERROR(err);
    
    err = mcl_task_set_arg(
        hdl, 1, (void*)edge_offsets, sizeof(uint32_t) * (numVerts+1), 
        MCL_ARG_BUFFER | MCL_ARG_INPUT
    );
    MCL_CHECK_ERROR(err);
    
    err = mcl_task_set_arg(
        hdl, 2, (void*)edge_list, sizeof(uint32_t) * adj_list_length, 
        MCL_ARG_BUFFER | MCL_ARG_INPUT
    );
    MCL_CHECK_ERROR(err);

    err = mcl_task_set_arg(
        hdl, 3, (void*)&W_SZ, sizeof(uint32_t), MCL_ARG_SCALAR
    );
    MCL_CHECK_ERROR(err);

    err = mcl_task_set_arg(
        hdl, 4, (void*)&CHUNK_SZ, sizeof(uint32_t), MCL_ARG_SCALAR
    );
    MCL_CHECK_ERROR(err);

    err = mcl_task_set_arg(
        hdl, 5, (void*)&numVerts, sizeof(uint32_t), MCL_ARG_SCALAR
    );
    MCL_CHECK_ERROR(err);

    err = mcl_task_set_arg(
        hdl, 6, (void*)&iters, sizeof(uint32_t), MCL_ARG_SCALAR
    );
    MCL_CHECK_ERROR(err);

    err = mcl_task_set_arg(
        hdl, 7, (void*)flag, sizeof(int), 
        MCL_ARG_BUFFER | MCL_ARG_INPUT | MCL_ARG_OUTPUT
    );
    MCL_CHECK_ERROR(err);

    err = mcl_exec(hdl, pes, lsize, MCL_TASK_ANY);
    MCL_CHECK_ERROR(err);

    return hdl;
}

// ****************************************************************************
// Function: RunTest1
//
// Purpose:
//   Runs the BFS benchmark using method 1 (IIIT-BFS method)
//
// Arguments:
//   resultDB: the benchmark stores its results in this ResultDatabase
//   op: the options parser / parameter database
//   G: input graph
//
// Returns:  nothing
// Programmer: Alok Kamatar
// Creation: June 16, 2020
//
// Modifications:
//
// ****************************************************************************
void RunTest1(ResultDatabase& resultDB, OptionParser& op, Graph *G)
{
    srand48(13579862);
    int iter_ops[] = {16, 64, 256, 4096};
    int _iterations = op.getOptionInt("iterations");
    const uint32_t num_tasks = iter_ops[_iterations - 1];

    //get graph info
    cl_uint *edgeArray=G->GetEdgeOffsets();
    cl_uint *edgeArrayAux=G->GetEdgeList();
    cl_uint adj_list_length=G->GetAdjacencyListLength();
    cl_uint numVerts = G->GetNumVertices();
    cl_uint numEdges = G->GetNumEdges();

    int err;
    char* kernel_path = "bfs_iiit.cl";

    //initialize frontier
    uint32_t** frontier = new uint32_t*[num_tasks];
    uint32_t* source_vertex = new uint32_t[num_tasks];
    for(int i = 0; i < num_tasks; i++) {
        frontier[i] = new uint32_t[numVerts];
        for (int index=0; index < numVerts; index++) {
            frontier[i][index]=UINT_MAX;
        }
        source_vertex[i] = lrand48()%numVerts;
        frontier[i][source_vertex[i]] = 0;
    }
    int* flag = new int[num_tasks];
    int* completed = new int[num_tasks];
    int* iters = new int[num_tasks];
    for(int i = 0; i < num_tasks; i++){
        flag[i] = 0;
        completed[i] = 0;
        iters[i] = 0;
    }

    //Get the kernel configuration
    size_t global_work_size[MCL_DEV_DIMS] = {1, 1, 1};
    size_t local_work_size[MCL_DEV_DIMS] = {1, 1, 1};
    size_t maxWorkItemsPerGroup=mcl::findMaxWorkGroupSize();
    GetWorkSize(global_work_size, local_work_size, maxWorkItemsPerGroup, numVerts);

    //configurable parameters for processing nodes
    cl_int W_SZ=32;
    cl_int CHUNK_SZ=32;

    // Do not run for large graphs, will take forever
    uint32_t** cpu_cost; 
    if(op.getOptionBool("verify")){
        cpu_cost = new uint32_t*[num_tasks];
        //Perform cpu bfs traversal for verify results
        for(int i = 0; i < num_tasks; i++) {
            cpu_cost[i] = new uint32_t[numVerts];
            G->GetVertexLengths(cpu_cost[i],source_vertex[i]);
        }
    }
    
    //number of passes
    int passes = op.getOptionInt("passes");
    mcl_handle** tasks = new mcl_handle*[num_tasks];

    //Start the benchmark
    std::cout<<"Running Benchmark\n";
    for (int j=0;j<passes;j++) {
        if (j>0) {
            //Reset the arrays to perform BFS again
            for(int i = 0; i < num_tasks; i++) {
                for (int index=0; index < numVerts; index++) {
                    frontier[i][index]=UINT_MAX;
                }
                frontier[i][source_vertex[i]] = 0;
                flag[i] = 0;
                completed[i] = 0;
                iters[i] = 0;
            }
        }

        //iteration count
        uint32_t finished = 0;
        //start CPU Timer to measure total time taken to complete benchmark
        int cpu_bfs_timer = Timer::Start();

        for(int i = 0; i < num_tasks; i++){
            tasks[i] = RunKernel1(
                kernel_path, frontier[i], numVerts, edgeArray, edgeArrayAux, 
                adj_list_length, W_SZ, CHUNK_SZ, iters[i], flag+i, global_work_size,
                NULL
            );
            iters[i]++;
        }

        //while there are nodes to traverse
        //flag is set if nodes exist in frontier
        int i = 0;
        while (finished < num_tasks) {
            if(completed[i] != 1 && mcl_test(tasks[i]) == MCL_REQ_COMPLETED) {
                mcl_hdl_free(tasks[i]);
                if(flag[i] == 1) {
                    flag[i] = 0;
                    tasks[i] = RunKernel1(
                        kernel_path, frontier[i], numVerts, edgeArray, edgeArrayAux, 
                        adj_list_length, W_SZ, CHUNK_SZ, iters[i], flag+i, global_work_size,
                        NULL
                    );
                    iters[i]++;
                } else {
                    finished++;
                    completed[i] = 1;
                }
            }
            i = (i+1)%num_tasks;
        };

        //stop the CPU timer
        double result_time = Timer::Stop(cpu_bfs_timer, "cpu_bfs_timer");
        result_time /= num_tasks;

        //count number of visited vertices
        uint32_t numVisited=0;
        for(int i=0;i<num_tasks;i++) {
            for(int index = 0; index < numVerts; index++) {
                if(frontier[i][index]!=UINT_MAX)
                    numVisited++;
            }
        }

        //Verify Results against serial BFS
        uint32_t unmatched_verts = 0;
        bool dump_paths = op.getOptionBool("dump-pl");
        if(op.getOptionBool("verify")){
            for(int i = 0; i < num_tasks; i++) {
                unmatched_verts += verify_results(cpu_cost[i],frontier[i],numVerts,
                    dump_paths);
            }
        }
        
        

        //total size transferred
        float gbytes=
            sizeof(uint32_t)*numVerts*2+   //2 frontiers
            sizeof(uint32_t)*numVerts+         //frontiers out array
            sizeof(uint32_t)*(numVerts+1)+       //edgeArray
            sizeof(uint32_t)*adj_list_length;    //edgeArrayAux
        gbytes/=(1000. * 1000. * 1000.);

        //populate the result database
        char atts[1024];
        sprintf(atts,"v:%d_e:%d ",numVerts,adj_list_length);
        if(unmatched_verts==0) {
            resultDB.AddResult("BFS_total",atts,"s",result_time);
            resultDB.AddResult("BFS",atts,"GB/s",gbytes/result_time);
            resultDB.AddResult("BFS_teps",atts,"Edges/s", numEdges/result_time);
            resultDB.AddResult("BFS_PCIe",atts,"Gbs/s", gbytes/(result_time));
            resultDB.AddResult("BFS_visited_vertices", atts, "N",numVisited);

        } else {
            resultDB.AddResult("BFS_total",atts,"s",FLT_MAX);
            resultDB.AddResult("BFS_kernel_time",atts,"s",FLT_MAX);
            resultDB.AddResult("BFS",atts,"GB/s",FLT_MAX);
            resultDB.AddResult("BFS_teps",atts,"Edges/s",FLT_MAX);
            resultDB.AddResult("BFS_PCIe",atts,"Gbs/s", FLT_MAX);
            resultDB.AddResult("BFS_visited_vertices", atts, "N",FLT_MAX);
        }

        std::cout<<"Test ";
        if(unmatched_verts==0) {
            std::cout<<"Passed\n";
        } else {
            std::cout<<"Failed\n";
            break;
        }
    }

    //Clean up
    delete[] flag;
    delete[] completed;
    delete[] iters;
    delete[] source_vertex;
    delete[] tasks;

    for(int i = 0; i < num_tasks; i++) {
        //delete[] cpu_cost[i];
        delete[] frontier[i];
    }
    //delete[] cpu_cost;
    delete[] frontier;
    
    std::cout << "Finished BFS test" << std::endl;
}

// ****************************************************************************
// Function: RunBenchmark
//
// Purpose:
//   Executes the BFS benchmark
//
// Arguments:
//   devcpp: opencl device
//   ctxcpp: the opencl context
//   queuecpp: the opencl command queue
//   resultDB: results from the benchmark are stored in this db
//   op: the options parser / parameter database
//
// Returns:  nothing
// Programmer: Aditya Sarwade
// Creation: June 16, 2011
//
// Modifications:
//
// ****************************************************************************
void RunBenchmark(ResultDatabase &resultDB, OptionParser &op)
{
    //adjacency list variables
    cl_mem h_edge,h_edgeAux;
    //number of vertices and edges in graph
    cl_uint numVerts=0,numEdges=0;
    //variable to get error code
    cl_int err_code;
    //Get the graph filename
    string inFileName = op.getOptionString("graph_file");
    //max degree in graph
    unsigned int max_deg=0;

    //Create graph
    Graph *G=new Graph();

    //Get pointers to edge offsets and edge list
    uint32_t **edge_ptr1 = G->GetEdgeOffsetsPtr();
    uint32_t **edge_ptr2 = G->GetEdgeListPtr();

    //Load simple k-way tree or from a file
    if (inFileName == "random")
    {
        //Load simple k-way tree
        //prob size specifies number of vertices
        cl_uint prob_sizes[4] = { 100,10000,1000000,100000000 };
       
        //Check for a valid size option and exit if not found
        if((op.getOptionInt("size") > 4) || (op.getOptionInt("size") <= 0))
        {
          cout<<"Please use a size between 1 and 4"<<endl;
          cout<<"Exiting..."<<endl;
          return;
        }

        numVerts = prob_sizes[op.getOptionInt("size")-1];
        int avg_degree = op.getOptionInt("degree");
        if(avg_degree<1)
            avg_degree=1;

        *edge_ptr1 = (uint32_t*)malloc(sizeof(uint32_t) * (numVerts+1));

        *edge_ptr2 = (uint32_t*)malloc(sizeof(uint32_t) * (numVerts*(avg_degree+1)));

        //Generate simple tree
        G->GenerateSimpleKWayGraph(numVerts,avg_degree);
    }
    else
    {
        //Read number of vertices and edges from first line of graph
        FILE *fp=fopen(inFileName.c_str(),"r");
        if(fp==NULL)
        {
            std::cout<<"\nFile not found!!!";
            return;
        }
        const char delimiters[]=" \n";
        char charBuf[MAX_LINE_LENGTH];
        fgets(charBuf,MAX_LINE_LENGTH,fp);
        char *temp_token = strtok (charBuf, delimiters);

        while(temp_token[0]=='%')
        {
            fgets(charBuf,MAX_LINE_LENGTH,fp);
            temp_token = strtok (charBuf, delimiters);
        }
        numVerts=atoi(temp_token);
        temp_token = strtok (NULL, delimiters);
        numEdges=atoi(temp_token);
        fclose(fp);

        *edge_ptr1 = (cl_uint *)malloc(sizeof(uint32_t)*(numVerts+1));
        *edge_ptr2 = (uint32_t*)malloc(sizeof(uint32_t)*(numEdges*2));

        //Load metis graph
        G->LoadMetisGraph(inFileName.c_str());

    }

    std::cout<<"Vertices: "<<G->GetNumVertices() << endl;
    std::cout<<"Edges: "<<G->GetNumEdges() << endl;
    unsigned long workers = op.getOptionInt("workers");
    //Run the test according to specified method
    int err = mcl_init(workers, 0x0);
    MCL_CHECK_ERROR(err);
    RunTest1(resultDB, op, G);
    std::cout << "Calling MCL finit" << endl;
    err = mcl_finit();
    MCL_CHECK_ERROR(err);
    std::cout << "Returned from MCL finit" << endl;

    free(*edge_ptr1);
    free(*edge_ptr2);
    delete G;
}

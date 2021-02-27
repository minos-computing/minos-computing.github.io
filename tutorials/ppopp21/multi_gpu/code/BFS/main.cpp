
#include <iostream>
#include <cstdlib>
#include <minos.h>

#include "ResultDatabase.h"
#include "OptionParser.h"
#include "utils.h"

using namespace std;

// Forward Declarations
void addBenchmarkSpecOptions(OptionParser &op);
void RunBenchmark(ResultDatabase &resultDB, OptionParser &op);

// ****************************************************************************
// Function: EnumerateDevices
//
// Purpose:
//   This function queries MCL about the available devices in the system, prints
//   those results to standard out
//
// Arguments:
//
// Returns:  nothing
//
// Programmer: Alok Kamatar
// Creation: Fri May 29 2020
// ****************************************************************************
void EnumerateDevices(void)
{
    mcl_init(1UL, 0x0);   
    int deviceCount = mcl_get_ndev();
    cout << "Number of devices = " << deviceCount << "\n";

    for (int device = 0; device < deviceCount; ++device)
    {
        mcl_device_info deviceProp;
        mcl_get_dev(device, &deviceProp);
    
        cout << "Device " << device << ":\n";
        cout << "  name               = '" << deviceProp.name << "'"
                << endl;
        cout << "  totalGlobalMem     = " << HumanReadable(
                deviceProp.mem_size) << endl;
        cout << "  Work Item size           = " << deviceProp.wisize << endl;
        cout << "  work Group size           = " << deviceProp.wgsize << endl;
        cout << "  Processing Elements = " << deviceProp.pes
                << endl;
        cout << "  Type             = " << deviceProp.type   << endl;
    }
    mcl_finit();
}

// ****************************************************************************
// Function: main
//
// Purpose:
//   The main function takes care of initialization,  then
//   performs the benchmark and prints results.
//
// Arguments:
//
//
// Programmer: Alok Kamatar
// Creation: 6/16/20
// ****************************************************************************
int main(int argc, char *argv[])
{
    int ret = 0;
    bool noprompt = false;

    try
    {        // Get args
        OptionParser op;

        //Add shared options to the parser
        op.addOption("verbose", OPT_BOOL, "", "enable verbose output", 'v');
        op.addOption("passes", OPT_INT, "10", "specify number of passes", 'n');
        op.addOption("size", OPT_INT, "1", "specify problem size", 's');
        op.addOption("infoDevices", OPT_BOOL, "",
                "show info for available platforms and devices", 'i');
        op.addOption("double", OPT_INT, "1", "specify double, 0 - no support, 1 - K support, 2 - AMD support", 'd');
        op.addOption("iterations", OPT_INT, "2", "Iterations (1-4) per pass to submit asychronously", 'r');
        op.addOption("workers", OPT_INT, "1", "Number of MCL workers", 'w');
#ifdef _WIN32
        op.addOption("noprompt", OPT_BOOL, "", "don't wait for prompt at program exit");
#endif

        addBenchmarkSpecOptions(op);
        if (!op.parse(argc, argv))
        {
            op.usage();
            return (op.HelpRequested() ? 0 : 1);
        }

        bool verbose = op.getOptionBool("verbose");
        bool infoDev = op.getOptionBool("infoDevices");
#ifdef _WIN32
        noprompt = op.getOptionBool("noprompt");
#endif

        if( infoDev )
        {
            EnumerateDevices();
            return 0;
        }
        ResultDatabase resultDB;

        // Run the benchmark
        RunBenchmark(resultDB, op);
        resultDB.DumpDetailed(cout);
    }
    catch( std::exception& e )
    {
        std::cerr << e.what() << std::endl;
        ret = 1;
    }
    catch( ... )
    {
        ret = 1;
    }

#ifdef _WIN32
    if (!noprompt)
    {
        cout << "Press return to exit\n";
        cin.get();
    }
#endif

    return ret;
}

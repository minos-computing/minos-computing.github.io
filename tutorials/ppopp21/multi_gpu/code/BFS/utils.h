#ifndef MCLUTILS_H
#define MCLUTILS_H

#include <minos.h>
#include <iostream>
#include <cstring>
#include <unistd.h>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <sstream>
#include <math.h>


#define BILLION 1000000000ULL
#define tdiff(end,start) BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec
#define K_DOUBLE_OPTION     1
#define AMD_DOUBLE_OPTION   2

namespace mcl {
    /** Finds the largest possible work group size **/
    size_t findMaxWorkGroupSize(void);

    /** Converts a string of openCL code to a file usable by MCL **/
    char* createKernelFile(const char* src);
    /** Removes kernel files, frees string memory and clears vector **/
    void clearKernelFiles(void);

    /** Runs an MCL kernel and waits for its completion, 
    /*  returns execution time ins seconds */
    double time_exec(mcl_handle* hdl, uint64_t* pes, uint64_t* lsize, uint64_t flags);
};



// decide which timer type we are supposed to use
#if defined(_WIN32)
#    define TIMEINFO _timeb
#else
#    define TIMEINFO timespec
#endif


// ****************************************************************************
//  Class:  Timer
//
//  Purpose:
//    Encapsulated a set of hierarchical timers.  Starting a timer
//    returns a handle to a timer.  Pass this handle, and a description,
//    into the timer Stop routine.  Timers can nest and output will
//    be displayed in a tree format.
//
//    Externally, Timer represents time in units of seconds.
//
//  Programmer:  Jeremy Meredith
//  Creation:    August  6, 2004
//
// ****************************************************************************
class Timer
{
  public:
    static Timer *Instance();

    static int    Start();

    // Returns time since start of corresponding timer (determined by handle),
    // in seconds.
    static double Stop(int handle, const std::string &descr);
    static void   Insert(const std::string &descr, double value);

    static void   Dump(std::ostream&);

  private:

    int    real_Start();
    double real_Stop(int, const std::string &);
    void   real_Insert(const std::string &descr, double value);
    void   real_Dump(std::ostream&);

    Timer();
    ~Timer();

    static Timer *instance;

    std::vector<TIMEINFO>    startTimes;
    std::vector<double>      timeLengths;
    std::vector<std::string> descriptions;
    int                      currentActiveTimers;
};

std::string HumanReadable(long long value, long long *rounding=0);
std::vector<std::string> SplitValues(const std::string &buff, char delim);


// ********************************************************
// Function: MCLErrorString
// Author: Alok Kamatar
// Description: 
// Converts an MCL error code to a human readable string
// ********************************************************
inline const char *MCLErrorString(int err)
{
    switch (-1*err)
    {
      case MCL_ERR_INVARG:      return "MCL Invalid Argument";
      case MCL_ERR_MEMALLOC:    return "MCL Unable to allocate memory";
      case MCL_ERR_INVREQ:      return "MCL Invalid Request";
      case MCL_ERR_INVPES:      return "MCL Invalid PES";
      case MCL_ERR_INVKER:      return "MCL Invalid kernel";
      case MCL_ERR_INVDEV:      return "MCL Invalid Device";
      case MCL_ERR_SRVCOMM:     return "MCL Server Communication Error";
      case MCL_ERR_INVTSK:      return "MCL Invalid Task";
      case MCL_ERR_MEMCOPY:     return "MCL Unable to copy memory";
      case MCL_ERR_EXEC:        return "MCL Exec Error";
      case MCL_ERR_INVPRG:      return "MCL Invalid Program";
      default:                  return "UNKNOWN";
  }
}

// ********************************************************
// Function: MCL_CHECK_ERROR
// Author: Alok Kamatar
// Description: 
// Checks if err is an MCL error code. If it is, prints the
// error string, removes the current kernel file, closes the
// connection to the MCL scheduler and exits
// ********************************************************
#define MCL_CHECK_ERROR(err) \
    {                       \
        if (err != 0)                  \
        { \
            std::cerr << __FILE__ << ':' << __LINE__ << ": " << MCLErrorString(err) << std::endl; \
            mcl::clearKernelFiles(); \
            mcl_finit(); \
            exit(1); \
        } \
    }\


#ifdef _WIN32

// On Windows, srand48 and drand48 don't exist.
// Create convenience routines that use srand/rand
// and let developers continue to use the -48 versions.

inline void srand48(unsigned int seed)
{
    srand(seed);
}

inline double drand48()
{
    return double(rand()) / RAND_MAX;
}

#endif // _WIN32


#endif /* !MCLUTILS_H */
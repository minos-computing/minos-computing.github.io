#include "utils.h"
#include <iostream>
#include <cstring>
#include <unistd.h>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <sstream>
#include <math.h>

using std::cerr;
using std::endl;
using std::max;

std::vector<char*> mcl_cur_kernel_paths;

// ****************************************************************************
// Function: findMaxWorkGroupSize
//
// Purpose:
//   This function queries MCL about the largest local work size possible
//
// Arguments:
//
// Returns:  nothing
//
// Programmer: Alok Kamatar
// Creation: Fri May 29 2020
// ****************************************************************************
size_t mcl::findMaxWorkGroupSize(void) {
    int deviceCount = mcl_get_ndev();
    size_t max_wgsize = 0;
    for (int device = 0; device < deviceCount; ++device)
    {
        mcl_device_info deviceProp;
        mcl_get_dev(device, &deviceProp);
        if(deviceProp.wgsize > max_wgsize)
            max_wgsize = deviceProp.wgsize;
    }
    return max_wgsize;
}

// ****************************************************************************
// Function: time_exec
//
// Purpose:
//   This function executes the task set up in an MCL handle, then waits for
//   the result.
//
// Arguments: Same as mcl_exec
//
// Returns:  The runtime of the task in seconds
//
// Programmer: Alok Kamatar
// Creation: Fri May 29 2020
// ****************************************************************************
double  mcl::time_exec(mcl_handle* hdl, uint64_t* pes, uint64_t* lsize, uint64_t flags) {
    struct timespec start, end;
    int err;

    clock_gettime(CLOCK_MONOTONIC,&start);
    err = mcl_exec(hdl, pes, lsize, flags);
    MCL_CHECK_ERROR(err);

    // Wait for the kernel to finish
    if((err = mcl_wait(hdl))) {
        MCL_CHECK_ERROR(err);
    }
    clock_gettime(CLOCK_MONOTONIC,&end);
    
    return ((double)tdiff(end, start))/BILLION;
}

// ****************************************************************************
// Function: createKernelFile
//
// Purpose:
//   Creates a temporary file with source code for use with an mcl kernel
//
// Arguments:
//    src - the source of the OpenCl/MCL kernel
//
// Returns:  the name of the temporary file
//
// Programmer: Alok Kamatar
// Creation: Fri May 29 2020
// ****************************************************************************
char* mcl::createKernelFile(const char* src) {
    char* kernel_path = strdup("tempKernelXXXXXX.cl");
    int fd;
    if((fd = mkstemps(kernel_path, 3)) == -1) {
        perror("Failed to create temporary kernel file");
        exit(1);
    }

    size_t written = 0;
    while(written < strlen(src)){
        ssize_t ret = write(fd, src, strlen(src) + 1);
        if(ret < 0){
            perror("Unable to write to kernel file");
            exit(1);
        }
        written += ret;
    }
    close(fd);
    mcl_cur_kernel_paths.push_back(kernel_path);
    return kernel_path;
}

// ****************************************************************************
// Function: clearKernelFiles
//
// Purpose:
//   Deletes the files created by createKernelFiles and frees associated memory
//
// Arguments:
//
// Returns:  nothing
//
// Programmer: Alok Kamatar
// Creation: 6/16/2020
// ****************************************************************************
void mcl::clearKernelFiles(void) {
    for(const auto& file : mcl_cur_kernel_paths) {
        remove(file);
        free(file);
    }
    mcl_cur_kernel_paths.clear();
}

// ----------------------------------------------------------------------------

Timer *Timer::instance = NULL;

// ----------------------------------------------------------------------------
static double
DiffTime(const struct TIMEINFO &startTime, const struct TIMEINFO &endTime)
{
#if defined(_WIN32)
    //
    // Figure out how many milliseconds between start and end times
    //
    int ms = (int) difftime(endTime.time, startTime.time);
    if (ms == 0)
    {
        ms = endTime.millitm - startTime.millitm;
    }
    else
    {
        ms =  ((ms - 1) * 1000);
        ms += (1000 - startTime.millitm) + endTime.millitm;
    }

    double seconds = (ms/1000.);
#else
    double seconds = double(endTime.tv_sec - startTime.tv_sec) +
                    double(endTime.tv_nsec - startTime.tv_nsec) / 1.0e9;
#endif
    return seconds;
}

static void
GetCurrentTimeInfo(struct TIMEINFO &timeInfo)
{
#if defined(_WIN32)
    _ftime(&timeInfo);
#else
    clock_gettime( CLOCK_REALTIME, &timeInfo );
#endif
}



// ****************************************************************************
//  Constructor:  Timer::Timer
//
//  Programmer:  Jeremy Meredith
//  Creation:    August  9, 2004
//
// ****************************************************************************
Timer::Timer()
{
    // Initialize some timer methods and reserve some space.
    startTimes.reserve(1000);
    timeLengths.reserve(1000);
    descriptions.reserve(1000);
    currentActiveTimers = 0;
}

// ****************************************************************************
//  Destructor:
//
//  Programmer:  Jeremy Meredith
//  Creation:    August  9, 2004
//
// ****************************************************************************
Timer::~Timer()
{
    // nothing to do
}

// ****************************************************************************
//  Method:  Timer::Instance
//
//  Purpose:
//    Return the timer singleton.
//
//  Arguments:
//
//
//  Programmer:  Jeremy Meredith
//  Creation:    August  9, 2004
//
// ****************************************************************************
Timer *Timer::Instance()
{
    if (!instance)
    {
        instance = new Timer;
    }
    return instance;
}

// ****************************************************************************
//  Method:  Timer::Start
//
//  Purpose:
//    Start a timer, and return a handle.
//
//  Arguments:
//    none
//
//  Programmer:  Jeremy Meredith
//  Creation:    August  9, 2004
//
// ****************************************************************************
int Timer::Start()
{
    return Instance()->real_Start();
}

// ****************************************************************************
//  Method:  Timer::Stop
//
//  Purpose:
//    Stop a timer and add its length to our list.
//
//  Arguments:
//    handle       a timer handle returned by Timer::Start
//    desription   a description for the event timed
//
//  Programmer:  Jeremy Meredith
//  Creation:    August  9, 2004
//
// ****************************************************************************
double Timer::Stop(int handle, const std::string &description)
{
    return Instance()->real_Stop(handle, description);
}

// ****************************************************************************
//  Method:  Timer::Insert
//
//  Purpose:
//    Add a user-generated (e.g. calculated) timing to the list
//
//  Arguments:
//    desription   a description for the event timed
//    value        the runtime to insert
//
//  Programmer:  Jeremy Meredith
//  Creation:    October 22, 2007
//
// ****************************************************************************
void Timer::Insert(const std::string &description, double value)
{
    Instance()->real_Insert(description, value);
}

// ****************************************************************************
//  Method:  Timer::Dump
//
//  Purpose:
//    Add timings to on ostream.
//
//  Arguments:
//    out        the stream to print to.
//
//  Programmer:  Jeremy Meredith
//  Creation:    August  9, 2004
//
// ****************************************************************************
void Timer::Dump(std::ostream &out)
{
    return Instance()->real_Dump(out);
}

// ****************************************************************************
//  Method:  Timer::real_Start
//
//  Purpose:
//    the true start routine
//
//  Arguments:
//    none
//
//  Programmer:  Jeremy Meredith
//  Creation:    August  9, 2004
//
// ****************************************************************************
int Timer::real_Start()
{
    int handle = startTimes.size();
    currentActiveTimers++;

    struct TIMEINFO t;
    GetCurrentTimeInfo(t);
    startTimes.push_back(t);

    return handle;
}

// ****************************************************************************
//  Method:  Timer::real_Stop
//
//  Purpose:
//    the true stop routine
//
//  Arguments:
//    handle       a timer handle returned by Timer::Start
//    desription   a description for the event timed
//
//  Programmer:  Jeremy Meredith
//  Creation:    August  9, 2004
//
// ****************************************************************************
double Timer::real_Stop(int handle, const std::string &description)
{
    if ((unsigned int)handle > startTimes.size())
    {
        cerr << "Invalid timer handle '"<<handle<<"'\n";
        exit(1);
    }

    struct TIMEINFO t;
    GetCurrentTimeInfo(t);
    double length = DiffTime(startTimes[handle], t);
    timeLengths.push_back(length);

    char str[2048];
    sprintf(str, "%*s%s", currentActiveTimers*3, " ", description.c_str());
    descriptions.push_back(str);

    currentActiveTimers--;
    return length;
}

// ****************************************************************************
//  Method:  Timer::real_Insert
//
//  Purpose:
//    the true insert routine
//
//  Arguments:
//    desription   a description for the event timed
//    value        the run time to insert
//
//  Programmer:  Jeremy Meredith
//  Creation:    October 22, 2007
//
// ****************************************************************************
void Timer::real_Insert(const std::string &description, double value)
{
#if 0 // can disable inserting just to make sure it isn't broken
    cerr << description << " " << value << endl;
#else
    timeLengths.push_back(value);

    char str[2048];
    sprintf(str, "%*s[%s]",
            (currentActiveTimers+1)*3, " ", description.c_str());
    descriptions.push_back(str);
#endif
}

// ****************************************************************************
//  Method:  Timer::real_Dump
//
//  Purpose:
//    the true dump routine
//
//  Arguments:
//    out        the stream to print to.
//
//  Programmer:  Jeremy Meredith
//  Creation:    August  9, 2004
//
// ****************************************************************************
void Timer::real_Dump(std::ostream &out)
{
    size_t maxlen = 0;
    for (unsigned int i=0; i<descriptions.size(); i++)
        maxlen = max(maxlen, descriptions[i].length());

    out << "\nTimings\n-------\n";
    for (unsigned int i=0; i<descriptions.size(); i++)
    {
        char desc[10000];
        sprintf(desc, "%-*s", (int)maxlen, descriptions[i].c_str());
        out << desc << " took " << timeLengths[i] << endl;
    }
}

// ****************************************************************************
// File:  Utility.h
//
// Purpose:
//   Various generic utility routines having to do with string and number
//   manipulation.
//
// Programmer:  Jeremy Meredith
// Creation:    September 18, 2009
// Modified:    Jan 2010, rothpc
//    Jeremy Meredith, Tue Oct  9 17:25:25 EDT 2012
//    Round is c99, not Windows-friendly.  Assuming we are using
//    positive values, replaced it with an equivalent of int(x+.5).
//
// ****************************************************************************

std::string HumanReadable(long long value, long long *rounding)
{
    std::ostringstream vstr;
    long long pVal;
    if (value>10ll*1024*1024*1024)
    {
        pVal = (long long)(0.5 + value/(1024.0*1024*1024));
        if (rounding)
            *rounding = pVal*1024*1024*1024 - value;
        vstr << pVal << 'G';
    }
    else if (value>10ll*1024*1024)
    {
        pVal = (long long)(0.5 + value/(1024.0*1024));
        if (rounding)
            *rounding = pVal*1024*1024 - value;
        vstr << pVal << 'M';
    }
    else if (value>10ll*1024)
    {
        pVal = (long long)(0.5 + value/(1024.0));
        if (rounding)
            *rounding = pVal*1024 - value;
        vstr << pVal << 'k';
    }
    else
    {
        if (rounding)
            *rounding = 0;
        vstr << value;
    }
    return vstr.str();
}

std::vector<std::string> SplitValues(const std::string &buff, char delim)
{
    std::vector<std::string> output;
    std::string tmp="";
    for (size_t i=0; i<buff.length(); i++)
    {
       if (buff[i] == delim)
       {
          if (!tmp.empty())
             output.push_back(tmp);
          tmp = "";
       }
       else
       {
          tmp += buff[i];
       }
    }
    if (!tmp.empty())
       output.push_back(tmp);

    return output;
}

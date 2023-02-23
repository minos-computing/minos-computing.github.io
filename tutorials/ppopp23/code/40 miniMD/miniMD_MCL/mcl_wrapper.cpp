/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National 
   Laboratories ( http://www.mantevo.org ). The primary 
   authors of miniMD are Steve Plimpton (sjplimp@sandia.gov) , Paul Crozier 
   (pscrozi@sandia.gov) and Christian Trott (crtrott@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you 
   can redistribute it and/or modify it under the terms of the GNU Lesser 
   General Public License as published by the Free Software Foundation; 
   either version 3 of the License, or (at your option) any later 
   version.
  
   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
   Lesser General Public License for more details.
    
   You should have received a copy of the GNU Lesser General Public 
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA.  See also: http://www.gnu.org/licenses/lgpl.txt .

   For questions, contact Paul S. Crozier (pscrozi@sandia.gov) or
   Christian Trott (crtrott@sandia.gov). 

   Please read the accompanying README and LICENSE files.
---------------------------------------------------------------------- */

#include "mcl_wrapper.h"
#include <cstring>
#include <cmath>
#include <cstdlib>

#ifndef MDPREC
#define MDPREC_STR "2"
#elif MDPREC == 2
#define MDPREC_STR "2"
#else
#define MDPREC_STR "1"
#endif

MCLWrapper::MCLWrapper()
{
	buffer = NULL;
	buffersize = 0;
	buffer_flags = 0;
	blockdim = 192;
}

MCLWrapper::~MCLWrapper()
{
	mcl_finit();
}

int MCLWrapper::Init(int argc, char** argv, int workers)
{
	mcl_init(workers, 0x0);
	return 0;
}

void* MCLWrapper::BufferResize(uint64_t newsize)
{
	if(buffer) free(buffer);
	buffer = malloc(newsize);
	buffersize = newsize;
	return buffer;
}

void* MCLWrapper::BufferGrow(uint64_t newsize)
{
	if(newsize<=buffersize) return buffer;

	if(buffer) free(buffer);
	buffer = malloc(newsize);
	buffersize = newsize;
	return buffer;
}

mcl_handle* MCLWrapper::LaunchKernel(const char* kernel_src, const char* kernel_name, int glob_threads, int nwait, mcl_handle** waitlist, int nargs, ...)
{
	va_list args;
	va_start(args,nargs);
	//fprintf(stderr, "Creating task.\n");
	int ret;
	mcl_handle* hdl = mcl_task_create();
	//fprintf(stderr, "Setting kernel.\n");
	ret = mcl_task_set_kernel(hdl, (char*)kernel_src, (char*)kernel_name, nargs, "-DMDPREC=" MDPREC_STR " -cl-mad-enable -DIAMONDEVICE", 0);
	for(int i=0; i<nargs; i++)
	{
		void* arg = va_arg(args,void*);
		unsigned int size = va_arg(args,unsigned int);
		uint64_t flags = va_arg(args, uint64_t);
		//fprintf(stderr, "Setting argument %d: size: %u, flags: %lu.\n", i, size, flags);
		mcl_task_set_arg(hdl, i, arg, size, flags);
	}
	va_end(args);

	size_t grid[3];
	grid[0] = ((glob_threads+blockdim-1)/blockdim)*blockdim;
	grid[1] = 1;
	grid[2] = 1;

	size_t block[3];
	block[0] = blockdim;
	block[1] = 1;
	block[2] = 1;

	//fprintf(stderr, "Executing task, block dim: %ld .\n", blockdim);
	ret = mcl_exec_with_dependencies(hdl, grid, block, MCL_TASK_GPU, nwait, waitlist);
	return hdl;
}

mcl_handle* MCLWrapper::LaunchKernelShared(const char* kernel_src, const char* kernel_name, int glob_threads, int nwait, mcl_handle** waitlist, int nargs, ...)
{
	va_list args;
	va_start(args,nargs);
	//fprintf(stderr, "Creating task.\n");
	int ret;
	mcl_handle* hdl = mcl_task_create_with_props(MCL_HDL_SHARED);
	//fprintf(stderr, "Setting kernel.\n");
	ret = mcl_task_set_kernel(hdl, (char*)kernel_src, (char*)kernel_name, nargs, "-DMDPREC=" MDPREC_STR " -cl-mad-enable -DIAMONDEVICE", 0);
	for(int i=0; i<nargs; i++)
	{
		void* arg = va_arg(args,void*);
		unsigned int size = va_arg(args,unsigned int);
		uint64_t flags = va_arg(args, uint64_t);
		//fprintf(stderr, "Setting argument %d: size: %u, flags: %lu.\n", i, size, flags);
		mcl_task_set_arg(hdl, i, arg, size, flags);
	}
	va_end(args);

	size_t grid[3];
	grid[0] = ((glob_threads+blockdim-1)/blockdim)*blockdim;
	grid[1] = 1;
	grid[2] = 1;

	size_t block[3];
	block[0] = blockdim;
	block[1] = 1;
	block[2] = 1;

	//fprintf(stderr, "Executing task, block dim: %ld .\n", blockdim);
	ret = mcl_exec_with_dependencies(hdl, grid, block, MCL_TASK_GPU, nwait, waitlist);
	return hdl;
}

mcl_handle* MCLWrapper::SetupKernel(const char* kernel_src, const char* kernel_name, uint64_t props, int nargs, ...)
{
	va_list args;
	va_start(args,nargs);

	mcl_handle* hdl = mcl_task_create_with_props(props);
	mcl_task_set_kernel(hdl, (char*)kernel_src, (char*)kernel_name, (nargs/3), NULL, 0);
	for(int i=0; i<nargs; i++)
	{
		void* arg = va_arg(args,void*);
		unsigned int size = va_arg(args,unsigned int);
		uint64_t flags = va_arg(args, uint64_t);
		mcl_task_set_arg(hdl, i, arg, size, flags);
	}
	va_end(args);

	return hdl;
}

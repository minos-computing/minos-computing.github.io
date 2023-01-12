/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator 

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov 

   See the README file in the top-level LAMMPS directory. 

   ----------------------------------------------------------------------- 

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/ 

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany 

   See the README file in the USER-CUDA directory. 

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef _MCL_DATA_H_
#define _MCL_DATA_H_


enum copy_mode {x, xx, xy, yx, xyz, xzy}; // yxz, yzx, zxy, zyx not yet implemented since they were not needed yet
//xx==x in atom_vec x is a member therefore copymode x produces compile errors
#include "mcl_wrapper.h"
#include <ctime>

#include <cstdio>
#include <typeinfo>

template <typename host_type, copy_mode mode>
class cMCLData
{
	protected:
	MCLWrapper* wrapper;
	unsigned int dim[3];
	uint64_t flags;
	host_type* host_data;
	host_type* temp_data;
	unsigned nbytes;
	bool is_continues;
	bool owns_data;

	public:
	cMCLData(MCLWrapper* mcl_wrapper, uint64_t flags, unsigned dim_x, unsigned dim_y=0, unsigned dim_z=0);
	cMCLData(MCLWrapper* mcl_wrapper, host_type* host_data, uint64_t flags, unsigned dim_x, unsigned dim_y=0, unsigned dim_z=0);
	~cMCLData();
	void setHostData(host_type* host_data);
	host_type* hostData() { return host_data;};

	void upload();
	void download();

	unsigned int* getDim() {return dim;};
	unsigned int devSize() {return nbytes;}
	uint64_t mclFlags() {return flags;}
	host_type* devData() {return temp_data ? temp_data : host_data;}
};



template <typename host_type, copy_mode mode>
cMCLData<host_type, mode>
::cMCLData(MCLWrapper* mcl_wrapper, uint64_t mcl_flags, unsigned dim_x, unsigned dim_y, unsigned dim_z)
{
	wrapper = mcl_wrapper;
	is_continues = true;
	owns_data = true;
	flags = mcl_flags;

	unsigned ndev;
	if((mode == x)||(mode==xx))
	{
		ndev = dim_x;
		dim[0] = dim_x;
		dim[1] = 0;
		dim[2] = 0;
		is_continues = true;
	}
	else if(mode == xy || mode == yx )
	{
		ndev = dim_x * dim_y;
		dim[0] = dim_x;
		dim[1] = dim_y;
		dim[2] = 0;
	}
	else
	{
		ndev = dim_x * dim_y * dim_z;
		dim[0] = dim_x;
		dim[1] = dim_y;
		dim[2] = dim_z;
	}

	nbytes = ndev * sizeof(host_type);
	if(nbytes<=0)
	{
		this->host_data=NULL;
		temp_data = NULL;
		return;
	}

	host_type* host_tmp = new host_type[ndev];
	if((mode==x)||(mode==xx))
		host_data = host_tmp;
	if((mode==xy)||(mode==yx))
	{
		host_type** host_tmpx = new host_type*[dim[0]];
		for(unsigned int i=0;i<dim[0];i++)
			host_tmpx[i] = &host_tmp[i*dim[1]];
		host_data = (host_type*) host_tmpx;
	}
	if((mode==xyz)||(mode==xzy))
	{
		host_type*** host_tmpx = new host_type**[dim[0]];
		for(unsigned int i=0;i<dim[0];i++)
		{
			host_tmpx[i] = new host_type*[dim[1]];
			for(unsigned int j=0;j<dim[1];j++)
				host_tmpx[i][j] = &host_tmp[(i*dim[1]+j)*dim[2]];
		}
		host_data = (host_type*) host_tmpx;
	}

	if((mode!=x)&&(mode!=xx)&&(!((mode==xy)&&is_continues))&&(!((mode==xyz)&&is_continues)))
	{
		temp_data = new host_type[ndev];
		mcl_register_buffer(temp_data, nbytes, mcl_flags);
	} else {
		temp_data = NULL;
		mcl_register_buffer(host_data, nbytes, mcl_flags);
	}
}

template <typename host_type, copy_mode mode>
cMCLData<host_type, mode>
::cMCLData(MCLWrapper* mcl_wrapper, host_type* host_data, uint64_t mcl_flags, unsigned dim_x, unsigned dim_y, unsigned dim_z)
{
	wrapper = mcl_wrapper;
	is_continues = false;
	owns_data = false;
	temp_data = NULL;
	flags = mcl_flags;

	this->host_data = host_data;
	unsigned ndev;
	if((mode == x)||(mode==xx))
	{
		ndev = dim_x;
		dim[0] = dim_x;
		dim[1] = 0;
		dim[2] = 0;
		is_continues = true;
	}
	else if(mode == xy || mode == yx )
	{
		ndev = dim_x * dim_y;
		dim[0] = dim_x;
		dim[1] = dim_y;
		dim[2] = 0;
	}
	else
	{
		ndev = dim_x * dim_y * dim_z;
		dim[0] = dim_x;
		dim[1] = dim_y;
		dim[2] = dim_z;
	}
	
	nbytes = ndev * sizeof(host_type);
	if(nbytes<=0)
	{
		this->host_data=NULL;
		temp_data=NULL;
		return;
	}
	
	if((mode!=x)&&(mode!=xx))
	{
		temp_data = new host_type[ndev];
		mcl_register_buffer(temp_data, nbytes, mcl_flags);
	} else {
		temp_data = NULL;
		mcl_register_buffer(host_data, nbytes, mcl_flags);
	}

	//dev_image = wrapper->AllocDevDataImageFloat4(0,imagesize);
}

template <typename host_type, copy_mode mode>
cMCLData<host_type, mode>
::~cMCLData()
{
	if(owns_data)
	{
		host_type* host_tmp;
		if((mode==x)||(mode==xx)) host_tmp=host_data;
		if((mode==xy)||(mode==yx))
		{
			host_tmp=&((host_type**)host_data)[0][0];
			delete [] (host_type**)host_data;
		}
		if((mode==xyz)||(mode==xzy))
		{
			host_tmp=&((host_type***)host_data)[0][0][0];
			for(unsigned int i=0;i<dim[0];i++)
			delete [] ((host_type***)host_data)[i];
			delete [] (host_type***)host_data;
		}
		delete [] host_tmp;
	}
	if(temp_data)
	{
		delete [] temp_data;
	}
}

template <typename host_type, copy_mode mode>
void cMCLData<host_type, mode>
::setHostData(host_type* host_data)
{
	this->host_data = host_data;
}

template <typename host_type, copy_mode mode>
void cMCLData<host_type, mode>
::upload()
{
	switch(mode)
	{
		case x:
		{	
			//temp_data = host_data;
			break;
		}
		
		case xx:
		{
			//temp_data = host_data;
			break;
		}

		case xy:
		{
			for(unsigned i=0; i<dim[0]; ++i)
			{
				host_type* temp = &temp_data[i * dim[1]];
				for(unsigned j=0; j< dim[1]; ++j)
				{
					temp[j] = reinterpret_cast<host_type**>(host_data)[i][j];
				}
			}
			break;
		}
		
		case yx:
		{
			for(unsigned j=0; j< dim[1]; ++j)
			{
				host_type* temp = &temp_data[j * dim[0]];
				for(unsigned i=0; i< dim[0]; ++i)
				{
					temp[i] = reinterpret_cast<host_type**>(host_data)[i][j];
				}
			}
			break;
		}	
		case xyz:
		{
			for(unsigned i=0; i < dim[0]; ++i)
			for(unsigned j=0; j < dim[1]; ++j)
			{
				host_type* temp = &temp_data[(i * dim[1] + j) * dim[2]];
				for(unsigned k=0; k < dim[2]; ++k)
				{
					temp[k] = reinterpret_cast<host_type***>(host_data)[i][j][k];
				}
			}
			break;
		}	

		case xzy:
		{
			for(unsigned i=0; i< dim[0]; ++i)
			for(unsigned k=0; k< dim[2]; ++k)
			{
				host_type* temp = &temp_data[(i* dim[2]+k)* dim[1]];
				for(unsigned j=0; j< dim[1]; ++j)
				{
					temp[j] = reinterpret_cast<host_type***>(host_data)[i][j][k];
				}
			}
			break;
		}	
	}
}

template <typename host_type, copy_mode mode>
void cMCLData<host_type, mode>
::download()
{
	switch(mode)
	{
		case x:
		case xx:
			break;

		case xy:
		{
			for(unsigned i=0; i< dim[0]; ++i)
			{
				host_type* temp = &temp_data[i *  dim[1]];
				for(unsigned j=0; j< dim[1]; ++j)
				{
					reinterpret_cast<host_type**>(host_data)[i][j] = temp[j];
				}
			}
			break;
		}
		
		case yx:
		{
			for(unsigned j=0; j< dim[1]; ++j)
			{
				host_type* temp = &temp_data[j* dim[0]];
				for(unsigned i=0; i< dim[0]; ++i)
				{
					reinterpret_cast<host_type**>(host_data)[i][j] = temp[i];
				}
			}
			break;
		}

		case xyz:
		{
			for(unsigned i=0; i< dim[0]; ++i)
			for(unsigned j=0; j< dim[1]; ++j)
			{
				host_type* temp = &temp_data[(i *  dim[1]+j)* dim[2]];
				for(unsigned k=0; k< dim[2]; ++k)
				{
					reinterpret_cast<host_type***>(host_data)[i][j][k] = temp[k];
				}
			}
			break;
		}

		case xzy:
		{
			for(unsigned i=0; i< dim[0]; ++i)
			for(unsigned k=0; k< dim[2]; ++k)
			{
				host_type* temp = &temp_data[(i *  dim[2]+k)* dim[1]];
				for(unsigned j=0; j< dim[1]; ++j)
				{
					reinterpret_cast<host_type***>(host_data)[i][j][k] = temp[j];
				}
			}
			break;
		}
	}
}
#endif // _MCL_DATA_H_

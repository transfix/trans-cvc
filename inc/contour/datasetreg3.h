/*
  Copyright 2011 The University of Texas at Austin

	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of MolSurf.

  MolSurf is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.


  MolSurf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MolSurf; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
// Datasetreg3 - representation for a 3D time-varying regular grid

#ifndef DATASET_REG3_H
#define DATASET_REG3_H

//#include <Utility/utility.h>
#include <contour/Dataset.h>
#include <contour/datareg3.h>

extern int verbose;

// Datasetreg3 - a 3D scalar time-varying regular grid of data
class Datasetreg3 : public Dataset
{
	private:
		Datareg3**	reg3;
	public:
		Datasetreg3(Data::DataType t, int ndata, int ntime, char* files[]);
		Datasetreg3(Data::DataType t, int ndata, int ntime, int* dim, u_char* data);
		~Datasetreg3()
		{
			if(reg3)
			{
				for(int i=0; i<ntime; i++)
				{
					delete reg3[i];
				}
				free(reg3);
			}
			free(min);
			free(max);
		};

		// member access methods

		int getNData(void)
		{
			return(reg3[0]->getNData());
		}

		// min, max for "0" variable at time step "t"
		float getMin(int t) const
		{
			return(reg3[t]->getMin());
		}
		float getMax(int t) const
		{
			return(reg3[t]->getMax());
		}

		// min, max for "0" variable, over all times
		float getMin()      const
		{
			return(min[0]);
		}
		float getMax()      const
		{
			return(max[0]);
		}

		// min, max for "j" variable, over all times
		float getMinFun(int j)      const
		{
			return(min[j]);
		}
		float getMaxFun(int j)      const
		{
			return(max[j]);
		}

		// "i" time step
		Data*	getData(int i)
		{
			return(reg3[i]);
		}
		Datareg3*	getMesh(int i)
		{
			return(reg3[i]);
		}
};

// Datasetreg3() - usual constructor, reads data from one or more files
inline Datasetreg3::Datasetreg3(Data::DataType t, int nd, int nt, char* fn[])
	: Dataset(t, nd, nt, fn)
{
	int i, j;
	meshtype = 5;
	min = (float*)malloc(sizeof(float)*nd);
	max = (float*)malloc(sizeof(float)*nd);
	for(i=0; i<nd; i++)
	{
		min[i] = 1e10;
		max[i] = -1e10;
	}
	reg3 = (Datareg3**)malloc(sizeof(Datareg3*)*nt);
	ncells = 0;
	maxcellindex=0;
	for(i = 0; i < nt; i++)				// timestep loop
	{
		// min[i] = 1e10;
		// max[i] = -1e10;
		if(verbose)
		{
			printf("loading file: %s\n", fn[i]);
		}
		reg3[i] = new Datareg3(t, nd, fn[i]);
		for(j=0; j<nd; j++)				// per variable loop
		{
			if(reg3[i]->getMin(j) < min[j])
			{
				min[j] = reg3[i]->getMin(j);
			}
			if(reg3[i]->getMax(j) > max[j])
			{
				max[j] = reg3[i]->getMax(j);
			}
		}
		if(reg3[i]->getNCells() > ncells)
		{
			ncells = reg3[i]->getNCells();
		}
		if(reg3[i]->maxCellIndex() > maxcellindex)
		{
			maxcellindex = reg3[i]->maxCellIndex();
		}
	}
	if(verbose)
		for(i = 0; i < nd; i++)
		{
			printf("variable[%d]: min=%f, max=%f\n",i, min[i],max[i]);
		}
}

// Datasetreg3() - alternative constructor for the libcontour library
inline Datasetreg3::Datasetreg3(Data::DataType t, int ndata, int ntime, int* dim, u_char* data)
	: Dataset(t, ndata, ntime, data)
{
	int i, j;				// timestep and variable indices
	int	size = 0;			// size of single timestep of data
	meshtype = 5;
	min = (float*)malloc(sizeof(float)*ndata);
	max = (float*)malloc(sizeof(float)*ndata);
	for(i=0; i<ndata; i++)
	{
		min[i] = 1e10;
		max[i] = -1e10;
	}
	reg3 = (Datareg3**)malloc(sizeof(Datareg3*)*ntime);
	ncells = 0;
	maxcellindex=0;
	switch(t)
	{
		case Data::UCHAR :
			size = dim[0]*dim[1]*dim[2] * ndata * sizeof(u_char);
			break;
		case Data::USHORT :
			size = dim[0]*dim[1]*dim[2] * ndata * sizeof(u_short);
			break;
		case Data::FLOAT :
			size = dim[0]*dim[1]*dim[2] * ndata * sizeof(float);
			break;
	}
	for(i = 0; i < ntime; i++)				// timestep loop
	{
		// min[i] = 1e10;
		// max[i] = -1e10;
		reg3[i] = new Datareg3(t, ndata, dim, data + i*size);
		for(j=0; j<ndata; j++)				// per variable loop
		{
			if(reg3[i]->getMin(j) < min[j])
			{
				min[j] = reg3[i]->getMin(j);
			}
			if(reg3[i]->getMax(j) > max[j])
			{
				max[j] = reg3[i]->getMax(j);
			}
		}
		if(reg3[i]->getNCells() > ncells)
		{
			ncells = reg3[i]->getNCells();
		}
		if(reg3[i]->maxCellIndex() > maxcellindex)
		{
			maxcellindex = reg3[i]->maxCellIndex();
		}
	}
	if(verbose)
		for(i = 0; i < ndata; i++)
		{
			printf("variable[%d]: min=%f, max=%f\n",i, min[i],max[i]);
		}
}

#endif

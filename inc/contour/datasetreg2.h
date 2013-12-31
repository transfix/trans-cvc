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
// Datasetreg2 - representation for a 2D scalar time-varying regular grid

#ifndef DATASET_REG2_H
#define DATASET_REG2_H

#include <contour/Dataset.h>
#include <contour/datareg2.h>

extern int verbose;

class Datasetreg2 : public Dataset
{
	private:				// data member

		Datareg2**	reg2;

	public:				// constructors and destructors

		Datasetreg2(Data::DataType t, int ndata, int ntime, char* files[]);
		Datasetreg2(Data::DataType t, int ndata, int ntime, int* dim,
					u_char* data);
		~Datasetreg2()
		{
			free(min);
			free(max);
		}

		// member access methods

		float getMin(int t) const
		{
			return(reg2[t]->getMin());
		}
		float getMax(int t) const
		{
			return(reg2[t]->getMax());
		}
		float getMin()      const
		{
			return(min[0]);
		}
		float getMax()      const
		{
			return(max[0]);
		}

		float getMinFun(int j)      const
		{
			return(min[j]);
		}
		float getMaxFun(int j)      const
		{
			return(max[j]);
		}


		Data*	getData(int i)
		{
			return(reg2[i]);
		}
		Datareg2*	getMesh(int i)
		{
			return(reg2[i]);
		}
};

// Datasetreg2() - usual constructor, reads data from one or more files
inline Datasetreg2::Datasetreg2(Data::DataType t, int nd, int nt, char* fn[])
	: Dataset(t, nd, nt, fn)
{
	int i,j;
	meshtype = 4;
	reg2 = (Datareg2**)malloc(sizeof(Datareg2*)*nt);
	//Joe R. -- 05/11/2010 -- somehow missed allocation of min/max
	min = (float*)malloc(sizeof(float)*nd);
	max = (float*)malloc(sizeof(float)*nd);
	for(j = 0; j < nd; j++)
	{
		min[j] = 1e10;
		max[j] = -1e10;
	}
	ncells = 0;
	maxcellindex = 0;
	for(i=0; i<nt; i++)
	{
		if(verbose)
		{
			printf("loading file: %s\n", fn[i]);
		}
		reg2[i] = new Datareg2(t, nd, fn[i]);
		for(j = 0; j < nd; j++)
		{
			if(reg2[i]->getMin() < min[j])
			{
				min[j] = reg2[i]->getMin();
			}
			if(reg2[i]->getMax() > max[j])
			{
				max[j] = reg2[i]->getMax();
			}
		}
		if(reg2[i]->getNCells() > ncells)
		{
			ncells = reg2[i]->getNCells();
		}
		if(reg2[i]->maxCellIndex() > maxcellindex)
		{
			maxcellindex = reg2[i]->maxCellIndex();
		}
	}
}

// Datasetreg2() - alternative constructor for the libcontour library
inline Datasetreg2::Datasetreg2(Data::DataType t, int ndata, int ntime,	int* dim, u_char* data)
	: Dataset(t, ndata, ntime, data)
{
	int size = 0;		// size of single timestep of data
	meshtype = 4;
	reg2 = (Datareg2**)malloc(sizeof(Datareg2*)*ntime);
	//Joe R. -- 05/11/2010 -- somehow missed allocation of min/max
	min = (float*)malloc(sizeof(float)*ndata);
	max = (float*)malloc(sizeof(float)*ndata);
	for(int j = 0; j < ndata; j++)
	{
		min[j] = 1e10;
		max[j] = -1e10;
	}
	ncells = 0;
	maxcellindex = 0;
	switch(t)
	{
		case Data::UCHAR :
			size = dim[0] * dim[1] * ndata * sizeof(u_char);
			break;
		case Data::USHORT :
			size = dim[0] * dim[1] * ndata * sizeof(u_short);
			break;
		case Data::FLOAT :
			size = dim[0] * dim[1] * ndata * sizeof(float);
			break;
	}
	for(int i=0; i<ntime; i++)
	{
		reg2[i] = new Datareg2(t, ndata, dim, data + i*size);
		for(int j = 0; j < ndata; j++)
		{
			if(reg2[i]->getMin() < min[j])
			{
				min[j] = reg2[i]->getMin();
			}
			if(reg2[i]->getMax() > max[j])
			{
				max[j] = reg2[i]->getMax();
			}
		}
		if(reg2[i]->getNCells() > ncells)
		{
			ncells = reg2[i]->getNCells();
		}
		if(reg2[i]->maxCellIndex() > maxcellindex)
		{
			maxcellindex = reg2[i]->maxCellIndex();
		}
	}
}

#endif

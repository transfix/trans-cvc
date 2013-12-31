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
#include <contour/datasetvol.h>

// Datasetvol() - usual constructor, reads data from one or more files
Datasetvol::Datasetvol(Data::DataType t, int nd, int nt, char* fn[])
	: Dataset(t, nd, nt, fn)
{
	int i, j;
	meshtype = 3;
	vol = (Datavol**)malloc(sizeof(Datavol*)*nt);
	for(j = 0; j < nd; j++)
	{
		min[j] = 1e10;
		max[j] = -1e10;
	}
	ncells = 0;
	for(i=0; i<nt; i++)
	{
		vol[i] = new Datavol(t, nd, fn[i]);
		for(j = 0; j < nd; j++)
		{
			if(vol[i]->getMin() < min[j])
			{
				min[j] = vol[i]->getMin();
			}
			if(vol[i]->getMax() > max[j])
			{
				max[j] = vol[i]->getMax();
			}
		}
		if(vol[i]->getNCells() > ncells)
		{
			ncells = vol[i]->getNCells();
		}
	}
	maxcellindex=ncells;
}

// Datasetvol() - called by the constructors to initialize the data
Datasetvol::Datasetvol(Data::DataType t, int ndata, int ntime, u_int nverts, u_int ncells, double* verts, u_int* cells, int* celladj, u_char* data)
	: Dataset(t, ndata, ntime, data)
{
	int		i;		// timestep index variable
	int		j;		// a variable index
	int		size = 0;	// size of single timestep of data
	meshtype = 3;
	vol = (Datavol**)malloc(sizeof(Datavol*)*ntime);
	for(j = 0; j < ndata; j++)
	{
		min[j] = 1e10;
		max[j] = -1e10;
	}
	Datasetvol::ncells = ncells;
	switch(t)
	{
		case Data::UCHAR :
			size = nverts * ndata * sizeof(u_char);
			break;
		case Data::USHORT :
			size = nverts * ndata * sizeof(u_short);
			break;
		case Data::FLOAT :
			size = nverts * ndata * sizeof(float);
			break;
	}
	for(i=0; i<ntime; i++)
	{
		vol[i] = new Datavol(t, ndata, nverts, ncells, verts, cells, celladj,
							 data + i*size);
		for(j = 0; j < ndata; j++)
		{
			if(vol[i]->getMin() < min[j])
			{
				min[j] = vol[i]->getMin();
			}
			if(vol[i]->getMax() > max[j])
			{
				max[j] = vol[i]->getMax();
			}
		}
		if(vol[i]->getNCells() > ncells)
		{
			ncells = vol[i]->getNCells();
		}
	}
	maxcellindex=ncells;
}

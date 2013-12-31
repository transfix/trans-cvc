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
// Dataset - representation for a scalar time-varying dataset

#ifndef DATASET_H
#define DATASET_H

//#include <Utility/utility.h>
#include <contour/data.h>

class Dataset
{
	private:
		Data::DataType	type;		// data type: uchar, ushort, float
		int		ndata;
		char**		filenames;	// data filenames

	protected:
		int		ntime;		// number of timesteps
		u_int	ncells;     		// number of cells
		int	meshtype;		// 2d unstr, reg2, 3d unstr, reg3
		int	maxcellindex;		// maximum number of cells
		float*	min, *max;		// min/max values for each variable

	public:	Dataset(Data::DataType t, int ndata, int ntime, char* files[]);
		Dataset(Data::DataType t, int ndata, int ntime, u_char* data);

		virtual ~Dataset()
		{
		}

		// member access methods
		Data::DataType	dataType(void) const
		{
			return(type);
		}
		int		meshType(void) const
		{
			return(meshtype);
		}
		int		nTime(void) const
		{
			return(ntime);
		}
		int		nData(void) const
		{
			return(ndata);
		}
		char**		fileNames(void) const
		{
			return(filenames);
		}

		virtual float getMin(int t) const = 0;	// min, max for "0" variable
		virtual float getMax(int t) const = 0;	// at time step "t"
		virtual float getMin() const = 0;	// min, max for "0" variable
		virtual float getMax() const = 0;	// over all times

		virtual float getMinFun(int f) const = 0; // min, max for "j" variable
		virtual float getMaxFun(int f) const = 0;// over all times (reg3 only)

		virtual Data* getData(int i) = 0;

		u_int         getNCells(void)
		{
			return(ncells);
		}
		int           maxCellIndex(void)
		{
			return(maxcellindex);
		}
};

// Dataset() - the usual constructor, initializes some data
inline Dataset::Dataset(Data::DataType t, int nd, int nt, char* fn[])
{
	type      = t;			// base type of the data
	ntime     = nt;			// number of time step
	ndata     = nd;			// number of data: add by fan
	filenames = fn;			// filenames containing the data
}

// Dataset() - alternative constructor for the libcontour library
inline Dataset::Dataset(Data::DataType t, int nd, int nt, u_char* data)
{
	type      = t;
	ntime     = nt;
	ndata     = nd;
	filenames = NULL;
}

#endif

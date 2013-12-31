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
// Data - representation for a single timestep of data
//	- data elements can be scalar or vector
//	- ndata/nvars gives dimension of each element (1 for scalar)

#ifndef DATA_H
#define DATA_H

//#include <Utility/utility.h>
#include <contour/endian_io.h>
#include <contour/basic.h>

#include <cstdlib>
#include <cstring>

const unsigned int BAD_INDEX = 0xFFFFFFFF;
extern int verbose;
extern void (*errorHandler)(char*, int);

// Data - a scalar dataset
class Data
{
	public:
		// DataType - supported types of raw data
		typedef enum
		{
			UCHAR=0,
			USHORT,
			FLOAT
		} DataType;

		// datatypes - a union to hold all of the above types
		typedef union
		{
			u_char* ucdata;
			u_short* usdata;
			float* fdata;
		} datatypes;

		// constructors and destructors
		Data(DataType t, int ndata, char* rawfile)
		{
			commonConstructor(t, ndata, rawfile);
		}
		Data(DataType t, int _ndata, u_char* data);
		virtual ~Data()
		{
			unsigned int f;
			for(f=0; f<ndata; f++)
			{
				if(data[f].ucdata)
				{
					std::free(data[f].ucdata);
					data[f].ucdata=NULL;
				}
			}
			if(data)
			{
				std::free(data);
				data=NULL;
			}
			if(min)
			{
				std::free(min);
				min=NULL;
			}
			if(max)
			{
				std::free(max);
				max = NULL;
			}
		}

		// member access methods

		u_int getNVerts(void)
		{
			return(nverts);
		}
		u_int getNCells(void)
		{
			return(ncells);
		}
		u_int getNData(void)
		{
			return(ndata);
		}

#ifdef VARIABLE
		static int funcontour;
		static int funcolor;
		static void setContourFun(int f)
		{
			funcontour = f;
		}
		static void setColorFun(int f)
		{
			funcolor = f;
		}
		static int getContourFun()
		{
			return funcontour;
		}
		static int getColorFun()
		{
			return funcolor;
		}

		int	variable;
		int	colorvar;
		void	setVariable(int	v) const
		{
			variable = v;
		}
		void	setColorVariable(int v) const
		{
			colorvar = v;
		}
		int	getVariable()
		{
			return variable;
		}
		int	getColorVariable()
		{
			return colorvar;
		}

		int	getVarMin(int var) const
		{
			return min[var];
		}
		int	getVarMax(int var) const
		{
			return max[var];
		}
		void	setVarMin(int var, float m)
		{
			min[var] = m;
		}
		void	setVarMax(int var, float m)
		{
			max[var] = m;
		}
#else
		int	funcontour;
		int	funcolor;

		void	setContourFun(int f)
		{
			funcontour = f;
		}
		void	setColorFun(int f)
		{
			funcolor = f;
		}
		int	getContourFun()
		{
			return funcontour;
		}
		int	getColorFun()
		{
			return funcolor;
		}
#endif /* of VARIABLE */

		// for topology
		static int funtopol1;
		static int funtopol2;
		static void setFunTopol1(int f)
		{
			funtopol1 = f;
		}
		static void setFunTopol2(int f)
		{
			funtopol2 = f;
		}
		static int getFunTopol1()
		{
			return funtopol1;
		}
		static int getFunTopol2()
		{
			return funtopol2;
		}

		float getMin() const
		{
			return(min[funcontour]);
		}
		float getMax() const
		{
			return(max[funcontour]);
		}
		void setMin(float m)
		{
			min[funcontour] = m;
		}
		void setMax(float m)
		{
			max[funcontour] = m;
		}

		// get and set the min/max values, modify by fan: f=0 not set
		float getMin(int f) const
		{
			return(min[f]);
		}
		float getMax(int f) const
		{
			return(max[f]);
		}
		void setMin(float m, int f)
		{
			min[f] = m;
		}
		void setMax(float m, int f)
		{
			max[f] = m;
		}

		virtual int maxCellIndex(void)
		{
			return(ncells);
		}

		virtual u_int getNCellVerts(void) = 0;
		virtual u_int getNCellFaces(void) = 0;
		virtual void getCellRange(int c, float& min, float& max) = 0;
		virtual void getFaceRange(u_int c, u_int f, float& min, float& max) = 0;
		virtual int getCellAdj(int c, int f) = 0;
		virtual u_int getAdjIndex(int c1, int c2)
		{
			u_int f;
			for(f=0; f<getNCellFaces(); f++)
				if(getCellAdj(c1,f) == c2)
				{
					return(f);
				}
			return(BAD_INDEX);
		}
		virtual u_int getCellVert(int c, int v) = 0;

		// signature functions
		virtual int getNFunctions(void) = 0; // # of signature functions
		virtual float* compFunction(int, u_int&, float**) = 0;
		virtual char* fName(int) = 0; // # signature function name

		// overloading
		virtual float* compFunction(int, u_int&, float***, float***, float***) = 0;

		// get and set the min/max values, modify by fan: f=0 not set
		float getValue(int i, int f)
		{
			if(type==UCHAR)
			{
				return(data[f].ucdata[i]);
			}
			else if(type==USHORT)
			{
				return(data[f].usdata[i]);
			}
			else if(type==FLOAT)
			{
				return(data[f].fdata[i]);
			}
			return(0.0);
		}
		void* getValues(int f)
		{
			if(type==UCHAR)
			{
				return(data[f].ucdata);
			}
			else if(type==USHORT)
			{
				return(data[f].usdata);
			}
			else if(type==FLOAT)
			{
				return(data[f].fdata);
			}
			return(NULL);
		}
		float getValue(int i)
		{
			if(type==UCHAR)
			{
				return(data[funcontour].ucdata[i]);
			}
			else if(type==USHORT)
			{
				return(data[funcontour].usdata[i]);
			}
			else if(type==FLOAT)
			{
				return(data[funcontour].fdata[i]);
			}
			return(0.0);
		}
		void* getValues()
		{
			if(type==UCHAR)
			{
				return(data[funcontour].ucdata);
			}
			else if(type==USHORT)
			{
				return(data[funcontour].usdata);
			}
			else if(type==FLOAT)
			{
				return(data[funcontour].fdata);
			}
			return(NULL);
		}

		int getDataSize(void)
		{
			if(type==UCHAR)
			{
				return(sizeof(u_char));
			}
			else if(type==USHORT)
			{
				return(sizeof(u_short));
			}
			else if(type==FLOAT)
			{
				return(sizeof(float));
			}
			return(0);
		}

		void getExtent(float min[3], float max[3])
		{
			std::memcpy(min, minext, sizeof(float[3]));
			std::memcpy(max, maxext, sizeof(float[3]));
		}

	protected:
		// called by the constructor functions
		inline void commonConstructor(DataType, int ndata, char*);
		inline void readData(void);
		inline void preprocessData(u_char* data);

		FILE*	fp;		// file descriptor for file containing data

	protected :

		u_int	nverts;		// number of vertices for unstruct'd data (?)
		u_int	ncells;		// number of cells for unstructured data
		u_int	ndata;		// number of variables in a grid element
		DataType	type;		// data of each variable in a grid element

		char*	filename;	// name of file containing the data

		float*	min;		// minimum and maximum data values
		float*	max;		// one for each variable (set externally!!)

		float	minext[3];	// dataset extents (bounding box?)
		float	maxext[3];

		datatypes* data;		// array containing data
		// width x height [x depth] x nvars
};


// commonConstructor() - called by the constructors to initialize the data
void Data::commonConstructor(DataType t, int _ndata, char* fn)
{
	type  = t;
	ndata = _ndata;
	filename = fn;
	min = NULL;
	max = NULL;
	if(ndata > 1)			// default setup add by fan
	{
		funcolor = 1;
		funcontour = 0;
		funtopol1 = 0;
		funtopol2 = 1;
	}
	else
	{
		funcontour = 0;
		funcolor = 0;
	}
#ifdef WIN32
	if(filename != NULL && (fp = fopen(filename, "rb")))
#else
	if(filename != NULL && (fp = fopen(filename, "r")))
#endif
	{
		if(verbose)
		{
			std::printf("reading extent\n");
		}
		fread_float(minext, sizeof(float), 3, fp);
		fread_float(maxext, sizeof(float), 3, fp);
		if(verbose)
			std::printf(" min = %f %f %f max = %f %f %f\n", minext[0], minext[1], minext[2],
				 maxext[0], maxext[1], maxext[2]);
		fread_int(&nverts, sizeof(int), 1, fp);
		fread_int(&ncells, sizeof(int), 1, fp);
		if(verbose)
		{
			std::printf("%d verts, %d cells\n", nverts, ncells);
		}
	}
	else
	{
		char	str[256];
		std::sprintf(str, "Data::commonConstructor: couldn't open file: %s",
				filename);
		errorHandler(str, true);
		fp = NULL;
	}
}

// readData() - a function to read and preprocess an array of scalar data
void Data::readData(void)
{
	u_int i;
	u_int f;
	float v;
	static float min_cutoff;
	data = (datatypes*)std::malloc(sizeof(datatypes)*ndata);
	if(verbose)
	{
		printf("reading data values\n");
	}
	switch(type)
	{
		case UCHAR:
			for(f=0; f<ndata; f++)
			{
				data[f].ucdata = (u_char*)std::malloc(sizeof(u_char)*nverts);
			}
			break;
		case USHORT:
			for(f=0; f<ndata; f++)
			{
				data[f].usdata = (u_short*)std::malloc(sizeof(u_short)*nverts);
			}
			break;
		case FLOAT:
			for(f=0; f<ndata; f++)
			{
				data[f].fdata = (float*)std::malloc(sizeof(float)*nverts);
			}
			break;
	}
	min = (float*)std::malloc(sizeof(float)*ndata);
	max = (float*)std::malloc(sizeof(float)*ndata);
	min_cutoff = 1e10;
	for(f = 0; f < ndata; f++)
	{
		if(verbose)
		{
			printf("reading size %d into %p\n", getDataSize(), getValues(f));
		}
		min[f] = 1e10;
		max[f] = -1e10;
		switch(type)
		{
			case UCHAR:
			        //freadSafely(getValues(f), getDataSize(), nverts, fp);
			        fread(getValues(f), getDataSize(), nverts, fp);
				break;
			case USHORT:
				fread_short(getValues(f), getDataSize(), nverts, fp);
				break;
			case FLOAT:
				fread_float(getValues(f), getDataSize(), nverts, fp);
				break;
		}
		for(i=0; i<nverts; i++)
		{
			if((v=getValue(i,f)) < min[f])
			{
				min[f] = v;
				if(min_cutoff > v)
				{
					min_cutoff = v;
					funcontour = f;
					funcolor = f;
				}
			}
			if(v > max[f])
			{
				max[f] = v;
			}
		}
		if(verbose)
		{
			printf("min = %f, max = %f\n", min[f], max[f]);
		}
	}
}

// preprocessData() - a function to preprocess an array of scalar data
void Data::preprocessData(u_char* data)
{
	u_int	i;
	u_int	f;
	float	v;
	static float min_cutoff;
	Data::data = (datatypes*)std::malloc(sizeof(datatypes)*ndata);
	switch(type)
	{
		case UCHAR:
			for(f = 0; f < ndata; f++)
			{
				Data::data[f].ucdata = ((u_char*)data) + f * nverts;
			}
			break;
		case USHORT:
			for(f = 0; f < ndata; f++)
			{
				Data::data[f].usdata = ((u_short*)data) + f * nverts;
			}
			break;
		case FLOAT:
			for(f = 0; f < ndata; f++)
			{
				Data::data[f].fdata = ((float*)data) + f * nverts;
			}
			break;
	}
	min = (float*)std::malloc(sizeof(float)*ndata);
	max = (float*)std::malloc(sizeof(float)*ndata);
	min_cutoff = 1e10;
	for(f = 0; f < ndata; f++)
	{
		if(verbose)
		{
			printf("preprocessing size %d into %p\n", getDataSize(), getValues(f));
		}
		min[f] = 1e10;
		max[f] = -1e10;
		for(i=0; i<nverts; i++)
		{
			if((v=getValue(i,f)) < min[f])
			{
				min[f] = v;
				if(min_cutoff > v)
				{
					min_cutoff = v;
					funcontour = f;
					funcolor = f;
				}
			}
			if(v > max[f])
			{
				max[f] = v;
			}
		}
		if(verbose)
		{
			printf("min = %f, max = %f\n", min[f], max[f]);
		}
	}
}

#endif

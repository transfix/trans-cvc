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
// Datareg2d - class for a 2d regular grid of scalar data
// Datareg2d - a volume of scalar data.

#ifndef DATAREG2D_H
#define DATAREG2D_H

#include <contour/data.h>

class dirSeedsReg2;
class regProp2;
class respProp2;
class seedChkr2;
class Conplotreg2;

class Datareg2 : public Data
{
	private:

		friend class dirSeedsReg2;		// private friends
		friend class regProp2;
		friend class respProp2;
		friend class seedChkr2;

		u_int dim[2];			// data members
		float orig[2];
		float span[2];

		int xbits, ybits;
		int xmask, ymask;
		int yshift;

	public:				// constructors and destructors

		Datareg2(Data::DataType, int ndata, char* rawfile);
		Datareg2(Data::DataType, int ndata, int* dim, u_char* data);
		~Datareg2() {}

		// member access methods

		void	getDim(u_int* v)
		{
			memcpy(v, dim,  2 * sizeof(u_int));
		}
		void	getOrig(float* v)
		{
			memcpy(v, orig, 2 * sizeof(float));
		}
		void	getSpan(float* v)
		{
			memcpy(v, span, 2 * sizeof(float));
		}

		void	setOrig(float* v)
		{
			memcpy(orig, v, 2 * sizeof(float));
		}
		void	setSpan(float* v)
		{
			memcpy(span, v, 2 * sizeof(float));
		}

		float xCoord(int i)
		{
			return(orig[0] + i*span[0]);
		}
		float yCoord(int j)
		{
			return(orig[1] + j*span[1]);
		}

		int maxCellIndex(void)
		{
			return(index2cell(dim[0]-2, dim[1]-2));
		}

		// signature function methods

		int getNFunctions(void)
		{
			return(4);
		}

		float* compFunction(int, u_int&, float**);
		float* compFunction(int, u_int&, float***,
							float***, float***)
		{
			return(NULL);   // add by fan
		}
		char* fName(int);

	protected :				// signature functions

		float* compLength(u_int&, float**);
		float* compArea(u_int&, float**);
		float* compMaxArea(u_int&, float**);
		float* compGradient(u_int&, float**);

	public :	// get data or gradient approximations (by differencing)

		void getCellValues(int c, float* val)
		{
			int i,j;
			cell2index(c,i,j);
			getCellValues(i,j,val);
		}
		void getCellValues(int i, int j, float* val)
		{
			val[0] = getValue(index2vert(i,j));
			val[1] = getValue(index2vert(i+1,j));
			val[2] = getValue(index2vert(i+1,j+1));
			val[3] = getValue(index2vert(i,j+1));
		}

		u_int   getCellVert(int c, int v)
		{
			int i, j;
			cell2index(c,i,j);
			switch(v)
			{
				case 0:
					return(index2vert(i,j));
				case 1:
					return(index2vert(i+1,j));
				case 2:
					return(index2vert(i+1,j+1));
				case 3:
					return(index2vert(i,j+1));
			}
			return(BAD_INDEX);
		}

		u_int getNCellVerts(void)
		{
			return(4);
		}
		u_int getNCellFaces(void)
		{
			return(4);
		}
		int getCellAdj(int c, int f)
		{
			int i, j;
			cell2index(c,i,j);
			switch(f)
			{
				case 0:
					return(j==0? -1 : index2cell(i,j-1));
				case 1:
					return(i==(signed int)dim[0]-2? -1 : index2cell(i+1,j));
				case 2:
					return(j==(signed int)dim[1]-2? -1 : index2cell(i,j+1));
				case 3:
					return(i==0? -1 : index2cell(i-1,j));
			}
			return(-1);
		}

		void getCellRange(int c, float& min, float& max)
		{
			float t;
			u_int i;
			max = min = getValue(getCellVert(c,0));
			for(i=1; i<getNCellVerts(); i++)
				if((t=getValue(getCellVert(c,i))) < min)
				{
					min = t;
				}
				else if(t > max)
				{
					max = t;
				}
		}

		void getFaceRange(u_int c, u_int f, float& min, float& max)
		{
			float t;
			min = max = getValue(getCellVert(c,f));
			if((t=getValue(getCellVert(c,f<3?f+1:0))) < min)
			{
				min = t;
			}
			else if(t > max)
			{
				max = t;
			}
		}

	protected:

		friend class Conplotreg2;

		void cell2index(int c, int& i, int& j)
		{
			int _left;
			i = c&xmask;
			_left = c>>xbits;
			j = _left&ymask;
		}

		int index2cell(int i, int j)
		{
			return((j << yshift) | i);
		}

		int index2vert(int i, int j)
		{
			return(i*dim[1] + j);
		}

};

#endif

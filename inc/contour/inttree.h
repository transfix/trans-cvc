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
// intTree.h - interval tree data structure

#ifndef INT_TREE_H
#define INT_TREE_H

//#include <Utility/utility.h>
#include <contour/CellSearch.h>

// Interval Tree Class
class IntTree : public CellSearch
{
	public:
		IntTree(u_int n = 0, float* v = NULL);
		~IntTree();

		void Init(u_int n, float* v);
		void InsertSeg(u_int cellid, float min, float max);
		void Dump(void);
		void Info(void);
		void Traverse(float, void (*f)(u_int, void*), void*);
		u_int getCells(float, u_int*);
		void Done(void);

	protected:
		u_int addSeed(u_int id, float mn, float mx)
		{
			u_int n = nseed++;
			if(n >= seedsize)
			{
				if(seedsize == 0)
				{
					seedsize=5;
					cellid = (u_int*)malloc(sizeof(u_int) * seedsize);
					min = (float*)malloc(sizeof(float) * seedsize);
					max = (float*)malloc(sizeof(float) * seedsize);
				}
				else
				{
					seedsize*=2;
					cellid = (u_int*)realloc(cellid, sizeof(u_int) * seedsize);
					min = (float*)realloc(min, sizeof(float) * seedsize);
					max = (float*)realloc(max, sizeof(float) * seedsize);
				}
			}
			cellid[n]=id;
			min[n]=mn;
			max[n]=mx;
			return(n);
		}

		u_int seedID(u_int n)
		{
			return(cellid[n]);
		}
		float seedMin(u_int n)
		{
			return(min[n]);
		}
		float seedMax(u_int n)
		{
			return(max[n]);
		}

		static int mincmp(const void*, const void*);
		static int maxcmp(const void*, const void*);
		static void travFun(u_int n, void* data);

	private:
		u_int nseed;
		u_int seedsize;
		u_int* cellid;
		float* min;
		float* max;

		int nleaf;
		float* vals;
		CellBucket* minlist;
		CellBucket* maxlist;

		void (*travCB)(u_int, void*);
		void* travData;
};

#endif

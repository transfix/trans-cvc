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
// conplot.h - class for preprocessing and extraction of isocontours

#ifndef CONPLOT_H
#define CONPLOT_H

//#include <Utility/utility.h>
#include <contour/BucketSearch.h>
#include <contour/CellQueue.h>
#include <contour/contour2d.h>
#include <contour/contour3d.h>
#include <contour/Dataset.h>
#include <contour/edgehash.h>
#include <contour/inttree.h>
#include <contour/range.h>
#include <contour/seedcells.h>
#include <contour/segtree.h>

#define USE_INT_TREE

class Conplot
{
	public:
		Conplot(Dataset* v);
		virtual ~Conplot();

		// preprocess the volume to compute the seed set
		void Preprocess(int t, void (*func)(int, void*) = NULL, void * = NULL);

		// extract an isosurface of the given value
		u_int Extract(float isovalue)
		{
			return(ExtractAll(isovalue));
		}

		// select a timestep
		void  setTime(int t);

		// mark a cell or test whether a cell is marked
		void TouchCell(u_int);
		int CellTouched(u_int);
		void ClearTouched(void);

		void ResetAll(void)
		{
			for(int t=0; t<data->nTime(); t++)
			{
				Reset(t);
			}
		}

		int getCells(float val)
		{
			return(tree[curtime].getCells(val, int_cells));
		}

		SeedCells*	getSeeds()
		{
			return &seeds[curtime];
		}

		Contour2d*	getContour2d()
		{
			return &contour2d[curtime];
		}

		Contour3d*	getContour3d()
		{
			return &contour3d[curtime];
		}

		// routines to support isocontour component output
		void BeginWrite(char* fprefix)
		{
			ncomponents = 0;
			filePrefix = fprefix;
		}

		inline void EndWrite()
		{
			filePrefix = NULL;
		}

	protected:
		// extract an isosurface
		u_int ExtractAll(float isovalue);
		// build the segment tree for the seed set
		void BuildSegTree(int t);
		virtual void Reset(int) = 0;
		virtual int  Size(int)  = 0;
		virtual int  isDone(int)  = 0;
		virtual void Done(int)  = 0;
		// track a contour from a seed cell
		virtual void TrackContour(float, int)  = 0;
		Dataset*	data;
		CellQueue queue;
		SeedCells* seeds;
		Contour2d* contour2d;
		Contour3d* contour3d;
		int	curtime;
		int	ncomponents;		// number of isocontour components
		char*	filePrefix;		// isocontour component file prefix
	private:
#ifdef USE_SEG_TREE
		SegTree* tree;
#elif defined USE_INT_TREE
		IntTree* tree;
#elif defined USE_BUCKETS
		BucketSearch* tree;
#endif
		u_int*	int_cells;
		u_char*	touched;
};

#endif

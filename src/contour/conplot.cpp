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
// conplot.C - preprocess and extract contours from 3d scalar data
// Copyright (c) 1997 Dan Schikore

#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif
#include <memory.h>
#include <string.h>
#ifndef WIN32
#include <unistd.h>
#include <sys/time.h>
#endif

#include <contour/BucketSearch.h>
#include <contour/Conplot.h>
#include <contour/inttree.h>
#include <contour/range.h>
#include <contour/segtree.h>
#include <stdlib.h>

extern int verbose;

// CellTouched() - test if a cell has been visited
int Conplot::CellTouched(u_int id)
{
	int byte;
	int bit;
	byte = id>>3;
	bit  = id &0x7;
	return(touched[byte] & (1<<bit));
}

// TouchCell() - mark a cell as visited
void Conplot::TouchCell(u_int id)
{
	int byte;
	int bit;
	byte = id>>3;
	bit  = id &0x7;
	touched[byte] |= (1<<bit);
}

// ClearTouched() - clear the bit array of 'touched' cells
void Conplot::ClearTouched(void)
{
	memset(touched, 0, sizeof(char)*((data->maxCellIndex()+7)>>3));
}


// Conplot() - create a contour plot for the given volume.
Conplot::Conplot(Dataset* d)
{
	data	= d;
	contour2d = NULL;
	contour3d = NULL;
	filePrefix = NULL;
	if(verbose)
	{
		printf("***** Data Characteristics\n");
		printf("cells: %d\n", data->getNCells());
		printf("*****\n");
	}
	// initialize the bit array of 'touched' (visited) cells
	touched = (u_char*)malloc(sizeof(u_char) * (data->maxCellIndex()+7)>>3);
	int_cells = (u_int*)malloc(sizeof(u_int) * data->maxCellIndex());
	if(verbose)
	{
		printf("initializing %d trees\n", data->nTime());
	}
	tree = NULL;
#ifdef USE_SEG_TREE
	tree = new SegTree[data->nTime()];
#elif defined USE_INT_TREE
	tree = new IntTree[data->nTime()];
#elif defined USE_BUCKETS
	tree = new BucketSearch[data->nTime()];
#endif
	seeds = new SeedCells[data->nTime()];	// initialize seed data array
	curtime = 0;
}

// ~Conplot() - destroy a plot
Conplot::~Conplot()
{
	delete [] tree;
	delete [] seeds;
	if(int_cells)
	{
		free(int_cells);
		int_cells = NULL;
	}
	if(touched)
	{
		free(touched);
		touched = NULL;
	}
}

void Conplot::setTime(int t)
{
	curtime = t;
}

// ExtractAll() - extract an isosurface by propagation in 3d.  Data is
//                assumed to reside in memory.
//            isovalue  = surface value of interest
u_int Conplot::ExtractAll(float isovalue)
{
	int n;
	int cur;
#ifdef TIME_SEARCH
	time_t start, finish;
	int t;
#endif
	if(isDone(curtime))
	{
		return(Size(curtime));
	}
#ifdef TIME_SEARCH
	start = clock();
	for(t=0; t<NQUERY; t++)
	{
		n=tree[curtime].getCells(isovalue, int_cells);
	}
	finish = clock();
	printf("%f seconds for %d queries\n", (finish-start)/(float)CLOCKS_PER_SEC, NQUERY);
	printf("%f seconds/query\n", (finish-start)/((float)(CLOCKS_PER_SEC)*NQUERY));
#endif
	// find the intersected seeds
	n = tree[curtime].getCells(isovalue, int_cells);
	if(verbose)
	{
		printf("%d intersected seeds\n", n);
	}
	// flush the old surface
	Reset(curtime);
	// clear bit array of 'touched' cells
	ClearTouched();
	// loop through the seeds in order
	for(cur = 0; cur < n; cur++)
	{
		if(!CellTouched(int_cells[cur]))
		{
			TouchCell(int_cells[cur]);
			TrackContour(isovalue, int_cells[cur]);
		}
	}
	if(verbose)
		if(contour3d)
		{
			printf("%d triangles\n", contour3d->getNTri());
		}
	Done(curtime);
#ifdef WRITE
	if(contour3d)
	{
		contour3d->write("output.tmesh");
	}
#endif
	return(Size(curtime));
}

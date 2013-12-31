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
// intTree.C - interval Tree manipulation
// Copyright (c) 1997 Dan Schikore

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <contour/inttree.h>

#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif
#ifndef WIN32
#include <unistd.h>
#endif

IntTree* global_tree = NULL;

#define DEBUG_TREENo

extern int verbose;

// IntTree() - construct a segment tree for the given range of values
IntTree::IntTree(u_int n, float* val)
{
	nseed=0;
	seedsize=0;
	cellid = NULL;
	min = NULL;
	max = NULL;
	if(n==0)
	{
		nleaf=0;
		vals=NULL;
		minlist = NULL;
		maxlist = NULL;
		return;
	}
	Init(n, val);
}

// Init() - Initialize the segment tree for the given set of values
void IntTree::Init(u_int n, float* val)
{
	nleaf = n;
	vals = (float*)malloc(sizeof(float)*nleaf);
	memcpy(vals, val, sizeof(float)*nleaf);
	minlist = new CellBucket[nleaf];
	maxlist = new CellBucket[nleaf];
}

// ~IntTree() - free storage for a segment tree
IntTree::~IntTree()
{
	if(verbose)
	{
		printf("IntTree destructor\n");
	}
	free(vals);
	/* should free inside buckets here */
	delete [] minlist;
	delete [] maxlist;
	if(min)
	{
		free(min);
		min = NULL;
	}
	if(max)
	{
		free(max);
		max = NULL;
	}
	if(cellid)
	{
		free(cellid);
		cellid = NULL;
	}
}

int IntTree::mincmp(const void* p1, const void* p2)
{
	u_int s1 = *((u_int*)p1);
	u_int s2 = *((u_int*)p2);
	if(global_tree->seedMin(s1) < global_tree->seedMin(s2))
	{
		return(-1);
	}
	if(global_tree->seedMin(s1) > global_tree->seedMin(s2))
	{
		return(1);
	}
	return(0);
}

int IntTree::maxcmp(const void* p1, const void* p2)
{
	u_int s1 = *((u_int*)p1);
	u_int s2 = *((u_int*)p2);
	if(global_tree->seedMax(s1) > global_tree->seedMax(s2))
	{
		return(-1);
	}
	if(global_tree->seedMax(s1) < global_tree->seedMax(s2))
	{
		return(1);
	}
	return(0);
}

void IntTree::Done(void)
{
	int i;
	global_tree = this;
	for(i=0; i<nleaf; i++)
	{
		qsort(minlist[i].getCells(), maxlist[i].nCells(), sizeof(u_int), mincmp);
		qsort(maxlist[i].getCells(), maxlist[i].nCells(), sizeof(u_int), maxcmp);
	}
}

// InsertSet() - recursively insert a segment into the tree
void IntTree::InsertSeg(u_int cellid, float min, float max)
{
	u_int left, right, root;
	u_int n;
	n = addSeed(cellid, min, max);
	left = 0;
	right = nleaf-1;
	while(left < right)
	{
		root = (left + right) >> 1;
		if(min <= vals[root] && vals[root] <= max)
		{
			minlist[root].insert(n);
			maxlist[root].insert(n);
			return;
		}
		if(min > vals[root])
		{
			left=root+1;
		}
		else   /* max < vals[root] */
		{
			right=root-1;
		}
	}
	// left == right
	minlist[left].insert(n);
	maxlist[left].insert(n);
}

void IntTree::travFun(u_int n, void* data)
{
	IntTree* tree = (IntTree*)data;
	tree->travCB(tree->seedID(n), tree->travData);
}

// Traverse() - Traverse the tree, calling the given function for
//              each stored segment containing the given value
void IntTree::Traverse(float val, void (*f)(u_int, void*), void* data)
{
	int left, right, root;
	left = 0;
	right = nleaf-1;
	travCB = f;
	travData = data;
	while(left < right)
	{
		root = (left + right) >> 1;
		if(vals[root] > val)
		{
			minlist[root].traverseCells(travFun, this);
			right=root-1;
		}
		else
		{
			maxlist[root].traverseCells(travFun, this);
			left=root+1;
		}
	}
}

// getCells() - traverse the tree, storing the cell id's of all
//              segments containing the given value in a list
u_int IntTree::getCells(float val, u_int* cells)
{
	int left, right, root;
	u_int ncells;
	int i;
	left = 0;
	right = nleaf-1;
	ncells=0;
	while(left < right)
	{
		root = (left + right) >> 1;
		if(vals[root] > val)
		{
			// for all cells in minlist, we know max > val
			// search the minlist for all cells with min < val
			for(i=0; i<minlist[root].nCells(); i++)
				if(seedMin(minlist[root].getCell(i)) < val)
				{
					cells[ncells++]  = seedID(minlist[root].getCell(i));
				}
				else
				{
					break;
				}
			right=root-1;
		}
		else
		{
			// for all cells in maxlist, we know min < val
			// search the maxlist for all cells with max > val
			for(i=0; i<maxlist[root].nCells(); i++)
				if(seedMax(maxlist[root].getCell(i)) > val)
				{
					cells[ncells++]  = seedID(maxlist[root].getCell(i));
				}
				else
				{
					break;
				}
			left=root+1;
		}
	}
	return(ncells);
}

// Dump() - dump the tree
void IntTree::Dump(void)
{
	int i, j;
	for(i=0; i<nleaf; i++)
	{
		printf("%d: value %f\n", i, vals[i]);
		minlist[i].dump((char*)"   MIN:");
		maxlist[i].dump((char*)"   MAX:");
		printf("seeds: ");
		for(j=0; j<minlist[i].nCells(); j++)
		{
			printf("(%d %f %f)", seedID(minlist[i].getCell(j)),
				   seedMin(minlist[i].getCell(j)),
				   seedMax(minlist[i].getCell(j)));
		}
		printf("\n");
	}
}

// Info() - print some stats about the tree
void IntTree::Info(void)
{
	int total, max;
	printf("______INTERVAL TREE STATS_____\n");
	printf("%d total segments\n", nseed);
	printf("%d values in segment tree (%d buckets)\n", nleaf, nleaf*2);
	total = max = 0;
	for(int i=0; i<nleaf; i++)
	{
		total += minlist[i].nCells();
		total += maxlist[i].nCells();
		if(minlist[i].nCells() > max)
		{
			max = minlist[i].nCells();
		}
		if(maxlist[i].nCells() > max)
		{
			max = maxlist[i].nCells();
		}
	}
	printf("total labels in tree: %d\n", total);
	printf("maximum labels in one list: %d\n", max);
	printf("______INTERVAL TREE STATS_____\n");
}

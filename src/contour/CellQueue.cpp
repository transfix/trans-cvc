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
// cellQueue   - queue of cell identifiers.  The circular queue dyanmically
//               resizes itself when full.  Elements in the queue are of
//               indicies of i,j,k (specialized for 3d structured grids)
// Copyright (c) 1997 Dan Schikore
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <contour/CellQueue.h>

// reset to empty
void CellQueue::Reset(void)
{
	nel = 0;
}

// check if queue is empty
int CellQueue::Empty(void)
{
	return(nel == 0);
}


// CellQueue() - create a new cell queue with elements of specified size
CellQueue::CellQueue(int size)
{
	nel    = 0;
	start  = 0;
	cellsize = size;
	cells    = (unsigned int*)malloc(sizeof(unsigned int) * cellsize);
}

// ~CellQueue() - free storage
CellQueue::~CellQueue()
{
	if(cells != NULL)
	{
		free(cells);
	}
}

// Add() - add an item to the queue
void CellQueue::Add(unsigned int c)
{
	int n;
	int oldsize;
	int atend;
	n = nel++;
	// resize the queue if needed
	if(nel > cellsize)
	{
		oldsize = cellsize;
		cellsize *= 2;
		cells = (unsigned int*)realloc(cells, sizeof(int) * cellsize);
		// move everything from 'start' to the end
		if(start != 0)
		{
			atend = oldsize - start;
			memmove(&cells[cellsize-atend], &cells[start], sizeof(unsigned int)*atend);
			start = cellsize-atend;
		}
	}
	n += start;
	if(n >= cellsize)
	{
		n-=cellsize;
	}
	cells[n] = c;
}

// Get() - return the top item from the queue
int CellQueue::Get(int& c)
{
	if(Peek(c) == -1)
	{
		return(-1);
	}
	Pop();
	return(1);
}

// Peek() - return the top item, but don't remove it
int CellQueue::Peek(int& c)
{
	if(nel == 0)
	{
		return(-1);
	}
	c = cells[start];
	return(1);
}

// Pop() - delete the top item in the queue
void CellQueue::Pop(void)
{
	start++;
	if(start == cellsize)
	{
		start=0;
	}
	nel--;
}

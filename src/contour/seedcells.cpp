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
// seedCells.C - maintain the list of seed cells

#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif
#include <contour/seedcells.h>

// SeedCells() - initialize the list
SeedCells::SeedCells()
{
	ncells = 0;
	cell_size = 10000;
	cells = (SeedCellP)malloc(sizeof(struct SeedCell) * cell_size);
}

// ~SeedCells() - free storage
SeedCells::~SeedCells()
{
	free(cells);
}

// AddSeed() - add a seed cell, increasing storage as necessary
int SeedCells::AddSeed(u_int id, float min, float max)
{
	int n = ncells++;
	if(n >= cell_size)
	{
		cell_size *= 2;
		cells = (SeedCellP)realloc(cells, sizeof(struct SeedCell) * cell_size);
	}
	cells[n].cell_id = id;
	cells[n].min = min;
	cells[n].max = max;
	return(n);
}

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
// seedCells.h - maintain a list of seed cells

#ifndef SEED_CELLS_H
#define SEED_CELLS_H

//#include <Utility/utility.h>
#include <contour/basic.h>

typedef struct SeedCell
{
	float min, max;
	u_int cell_id;
} *SeedCellP;

class SeedCells
{
	public:
		SeedCells();
		~SeedCells();

		int    getNCells(void)
		{
			return(ncells);
		}
		u_int  getCellID(int i)
		{
			return(cells[i].cell_id);
		}
		float  getMin(int i)
		{
			return(cells[i].min);
		}
		float  getMax(int i)
		{
			return(cells[i].max);
		}
		void   Clear(void)
		{
			ncells = 0;
		}
		SeedCell* getCellPointer()
		{
			return(cells);
		}

		int AddSeed(u_int, float, float);
		void AddToRange(u_int i, float mn, float mx)
		{
			if(mn < cells[i].min)
			{
				cells[i].min = mn;
			}
			if(mx > cells[i].max)
			{
				cells[i].max = mx;
			}
		}

	private:
		int ncells;
		int cell_size;
		SeedCellP cells;
};

#endif

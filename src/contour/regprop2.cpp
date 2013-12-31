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
// regProp2.C - preprocessing of 2d volumes for seed set extraction

#include <stdlib.h>
#include <memory.h>
#include <contour/regprop2.h>
#include <contour/datareg2.h>

#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif
#ifndef WIN32
#include <unistd.h>
#endif

extern int verbose;

// Preprocess() - build a segment tree for O(log n) queries
void regProp2::compSeeds(void)
{
	Datareg2& reg2 = (Datareg2&)data;
	int i, j;
	int xdim, ydim;
	float val[4];
	Range* _prop_x, *prop_x;
	Range prop_y;
	Range propagated;
	Range c_prop;
	Range responsibility, c_respons;
	Range delay;
	Range y_comp;
	float min_x, min_y, max_x, max_y;
	float min_in, max_in, min4, max4;
	int nseed;
	xdim = reg2.dim[0];
	ydim = reg2.dim[1];
	_prop_x = new Range[ydim];
	// proceed through the slices computing seeds
	nseed=0;
	// process the k'th slab
	for(i=0; i<xdim-1; i++)
		for(j=0; j<ydim-1; j++)
		{
			prop_x = &_prop_x[j];
			// load the voxel data
			reg2.getCellValues(i, j, val);
			min_x = MIN2(val[0], val[3]);
			max_x = MAX2(val[0], val[3]);
			min_y = MIN2(val[0], val[1]);
			max_y = MAX2(val[0], val[1]);
			// set the incoming values if on a border
			if(i==0)
			{
				prop_x->Set(min_x, max_x);
			}
			if(j==0)
			{
				prop_y.Set(min_y, max_y);
			}
			// merge incoming information
			y_comp = prop_y.Complement(min_y, max_y);
			propagated = prop_y + ((*prop_x)-y_comp);
			// compute complement of incoming ranges
			min_in = MIN2(min_x, min_y);
			max_in = MAX2(max_x, max_y);
			c_prop.Set(min_in,max_in);
			c_prop -= propagated;
			// compute responsibility ranges
			min4 = MIN2(min_in, val[2]);
			max4 = MAX2(max_in, val[2]);
			responsibility.Set(min4, max4);
			responsibility-=c_prop;
			c_respons = responsibility.Complement(min4, max4);
			// determine range which can be delayed
			delay.MakeEmpty();
			if(i < xdim-2)
				delay+=Range(MIN2(val[1], val[2]),
							 MAX2(val[1], val[2]));
			if(j < ydim-2)
				delay+=Range(MIN2(val[2], val[3]),
							 MAX2(val[2], val[3]));
			// test for propagation of entire responsibility range
			if(responsibility.Empty() || (!delay.Empty() &&
										  delay.MinAll() <= responsibility.MinAll() &&
										  delay.MaxAll() >= responsibility.MaxAll()))
			{
				// propagate first to the next x-slice
				if(i == xdim-2)
				{
					prop_x->MakeEmpty();
				}
				else
				{
					prop_x->Set(MIN2(val[1], val[2]), MAX2(val[1], val[2]));
					*prop_x-=c_respons;
				}
				c_respons += *prop_x;
				// all remaining propagated in y-dir
				if(j == ydim-2)
				{
					prop_y.MakeEmpty();
				}
				else
				{
					prop_y.Set(MIN2(val[2], val[3]), MAX2(val[2], val[3]));
					prop_y-= c_respons;
				}
			}
			else
			{
				// can't propagate all responsiblity, cell must be a seed
				seeds.AddSeed(reg2.index2cell(i,j), responsibility.MinAll(),
							  responsibility.MaxAll());
				nseed++;
				prop_y.MakeEmpty();
				prop_x->MakeEmpty();
			}
		}
	if(verbose)
	{
		printf("computed %d seeds\n", nseed);
	}
}

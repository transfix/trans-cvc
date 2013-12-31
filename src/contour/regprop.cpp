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
// regProp.C - preprocessing of 3d volumes for seed set extraction

#include <stdlib.h>
#include <contour/regprop.h>
#include <contour/datareg3.h>
#include <memory.h>

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
void regProp::compSeeds(void)
{
	Datareg3& reg3 = (Datareg3&)data;
	int i, j, k;
	int xdim, ydim, zdim;
	float val[8];
	Range* _prop_z, *prop_z;
	Range* _prop_y, *prop_y;
	Range prop_x;
	Range propagated;
	Range c_prop;
	Range responsibility, c_respons;
	Range delay;
	Range x_comp;
	float min_x, min_y, min_z, max_x, max_y, max_z;
	float min_in, max_in, min8, max8;
	int nseed;
	xdim = reg3.dim[0];
	ydim = reg3.dim[1];
	zdim = reg3.dim[2];
	_prop_z = new Range[xdim * ydim];
	_prop_y = new Range[xdim];
	// proceed through the slices computing seeds
	nseed=0;
	for(k=0; k<zdim-1; k++)
	{
		if(verbose)
			if(k % 10 == 0)
			{
				printf("slice %d, %d seeds\n", k, nseed);
			}
		// process the k'th slab
		for(j=0; j<ydim-1; j++)
			for(i=0; i<xdim-1; i++)
			{
				prop_y = &_prop_y[i];
				prop_z = &_prop_z[j*(xdim-1)+i];
				// load the voxel data
				reg3.getCellValues(i, j, k, val);
				min_x = MIN4(val[0], val[3], val[4], val[7]);
				max_x = MAX4(val[0], val[3], val[4], val[7]);
				min_y = MIN4(val[0], val[1], val[2], val[3]);
				max_y = MAX4(val[0], val[1], val[2], val[3]);
				min_z = MIN4(val[0], val[1], val[4], val[5]);
				max_z = MAX4(val[0], val[1], val[4], val[5]);
				// set the incoming values if on a border
				if(i==0)
				{
					prop_x.Set(min_x, max_x);
				}
				if(j==0)
				{
					prop_y->Set(min_y, max_y);
				}
				if(k==0)
				{
					prop_z->Set(min_z, max_z);
				}
				// merge incoming information
				x_comp = prop_x.Complement(min_x, max_x);
				propagated = prop_x + ((*prop_y)+(*prop_z)-x_comp);
				// compute complement of incoming ranges
				min_in = MIN3(min_x, min_y, min_z);
				max_in = MAX3(max_x, max_y, max_z);
				c_prop.Set(min_in,max_in);
				c_prop -= propagated;
				// compute responsibility ranges
				min8 = MIN2(min_in, val[6]);
				max8 = MAX2(max_in, val[6]);
				responsibility.Set(min8, max8);
				responsibility-=c_prop;
				c_respons = responsibility.Complement(min8, max8);
				// determine range which can be delayed
				delay.MakeEmpty();
				if(i < xdim-2)
					delay+=Range(MIN4(val[1], val[2], val[5], val[6]),
								 MAX4(val[1], val[2], val[5], val[6]));
				if(j < ydim-2)
					delay+=Range(MIN4(val[4], val[5], val[6], val[7]),
								 MAX4(val[4], val[5], val[6], val[7]));
				if(k < zdim-2)
					delay+=Range(MIN4(val[2], val[3], val[6], val[7]),
								 MAX4(val[2], val[3], val[6], val[7]));
				// test for propagation of entire responsibility range
				if(responsibility.Empty() || (!delay.Empty() &&
											  delay.MinAll() <= responsibility.MinAll() &&
											  delay.MaxAll() >= responsibility.MaxAll()))
				{
					// propagate first to the next z-slice
					if(k == zdim-2)
					{
						prop_z->MakeEmpty();
					}
					else
					{
						prop_z->Set(MIN4(val[2], val[3], val[6], val[7]),
									MAX4(val[2], val[3], val[6], val[7]));
						*prop_z -= c_respons;
					}
					c_respons += *prop_z;
					// propagate in y-direction next
					if(j == ydim-2)
					{
						prop_y->MakeEmpty();
					}
					else
					{
						prop_y->Set(MIN4(val[4], val[5], val[6], val[7]),
									MAX4(val[4], val[5], val[6], val[7]));
						*prop_y-=c_respons;
					}
					c_respons += *prop_y;
					// all remaining propagated in x-dir
					if(i == xdim-2)
					{
						prop_x.MakeEmpty();
					}
					else
					{
						prop_x.Set(MIN4(val[1], val[2], val[5], val[6]),
								   MAX4(val[1], val[2], val[5], val[6]));
						prop_x-= c_respons;
					}
				}
				else
				{
					// can't propagate all responsiblity, cell must be a seed
					seeds.AddSeed(reg3.index2cell(i,j,k), responsibility.MinAll(),
								  responsibility.MaxAll());
					nseed++;
					prop_z->MakeEmpty();
					prop_y->MakeEmpty();
					prop_x.MakeEmpty();
				}
			}
	}
	delete [] _prop_z;
	delete [] _prop_y;
}

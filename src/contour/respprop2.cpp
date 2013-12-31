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
// respProp2.C - preprocessing of 2d volumes for seed set extraction

#include <stdlib.h>
#include <memory.h>
#include <contour/respprop2.h>
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

void respProp2::compSeeds(void)
{
	Datareg2& reg2 = (Datareg2&)data;
	int i, j;
	int xdim, ydim;
	float val[4];
	Range prop, c_prop, done, resp, out;
	float min_x, min_y, max_x, max_y;
	int nseed;
	xdim = reg2.dim[0];
	ydim = reg2.dim[1];
	// proceed through the slices computing seeds
	nseed=0;
	// process the k'th slab
	for(i=0; i<xdim-1; i++)
		for(j=0; j<ydim-1; j++)
		{
			// load the voxel data
			reg2.getCellValues(i, j, val);
			min_x = MIN2(val[0], val[3]);
			max_x = MAX2(val[0], val[3]);
			min_y = MIN2(val[0], val[1]);
			max_y = MAX2(val[0], val[1]);
			// set the incoming values if on a border
			if(j==0)
			{
				prop.Set(min_y, max_y);
				c_prop.MakeEmpty();
			}
			if(i==0)
			{
				done.MakeEmpty();
				resp = Range(min_x, max_x);
			}
			else
			{
				done = Range(min_x, max_x);
				resp.MakeEmpty();
			}
			done += c_prop;
			resp = (prop + Range(MIN2(val[1],val[2]),MAX2(val[1],val[2]))) - done;
			if(j < ydim-2)
			{
				out = Range(MIN2(val[2], val[3]), MAX2(val[2], val[3]));
			}
			else
			{
				out.MakeEmpty();
			}
			// test for propagation of entire responsibility range
			if(resp.Empty() || (!out.Empty() &&
								out.MinAll() <= resp.MinAll() &&
								out.MaxAll() >= resp.MaxAll()))
			{
				prop = out - done;
				c_prop = out - prop;
			}
			else
			{
				// can't propagate all responsiblity, cell must be a seed
				seeds.AddSeed(reg2.index2cell(i,j), resp.MinAll(), resp.MaxAll());
				nseed++;
				prop.MakeEmpty();
				c_prop = out;
			}
		}
	if(verbose)
	{
		printf("computed %d seeds\n", nseed);
	}
}

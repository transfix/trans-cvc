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
// seedChkr2.C - preprocessing of 2d volumes for seed set extraction

#include <memory.h>
#include <stdlib.h>
#include <contour/seedchkr2.h>
#include <contour/datareg2.h>

#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif
#ifndef WIN32
#include <unistd.h>
#endif

void seedChkr2::compSeeds(void)
{
	Datareg2& reg2 = (Datareg2&)data;
	int i, j;
	int xdim, ydim;
	float val[4];
	float min4, max4;
	int nseed = 0;
	xdim = reg2.dim[0];
	ydim = reg2.dim[1];
	// proceed through the slices computing seeds
	// process the k'th slab
	for(i=0; i<xdim-1; i+=2)
		for(j=0; j<ydim-1; j+=2)
		{
			// load the voxel data
			reg2.getCellValues(i, j, val);
			min4 = MIN4(val[0], val[1], val[2], val[3]);
			max4 = MAX4(val[0], val[1], val[2], val[3]);
			seeds.AddSeed(reg2.index2cell(i,j), min4, max4);
			nseed++;
		}
	for(i=1; i<xdim-1; i+=2)
		for(j=1; j<ydim-1; j+=2)
		{
			// load the voxel data
			reg2.getCellValues(i, j, val);
			min4 = MIN4(val[0], val[1], val[2], val[3]);
			max4 = MAX4(val[0], val[1], val[2], val[3]);
			seeds.AddSeed(reg2.index2cell(i,j), min4, max4);
			nseed++;
		}
}

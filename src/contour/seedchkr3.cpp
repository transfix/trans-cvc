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
// seedChkr3.C - preprocessing of 3d volumes for seed set extraction

#include <stdlib.h>
#include <memory.h>
#include <contour/seedchkr3.h>
#include <contour/datareg3.h>

#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif
#ifndef WIN32
#include <unistd.h>
#endif

void seedChkr3::compSeeds(void)
{
	Datareg3& reg3 = (Datareg3&)data;
	int i, j, k;
	int xdim, ydim, zdim;
	float val[8];
	float min8, max8;
	int nseed = 0;
	xdim = reg3.dim[0];
	ydim = reg3.dim[1];
	zdim = reg3.dim[2];
	// proceed through the slices computing seeds
	// process the k'th slab
	for(i=0; i<xdim-1; i+=2)
		for(j=0; j<ydim-1; j+=2)
			for(k=0; k<zdim-1; k+=2)
			{
				// load the voxel data
				reg3.getCellValues(i, j, k, val);
				MIN8(min8, val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
				MAX8(max8, val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
				seeds.AddSeed(reg3.index2cell(i,j,k), min8, max8);
				nseed++;
			}
	for(i=1; i<xdim-1; i+=2)
		for(j=1; j<ydim-1; j+=2)
			for(k=1; k<zdim-1; k+=2)
			{
				// load the voxel data
				reg3.getCellValues(i, j, k, val);
				MIN8(min8, val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
				MAX8(max8, val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
				seeds.AddSeed(reg3.index2cell(i,j,k), min8, max8);
				nseed++;
			}
}

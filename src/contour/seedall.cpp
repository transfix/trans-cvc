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
// seedAll.C - preprocessing of 2d volumes for seed set extraction

#include <stdlib.h>
#include <memory.h>
#include <contour/seedall.h>
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

void seedAll::compSeeds(void)
{
	float min, max;
	// proceed through the slices computing seeds
	for(u_int c=0; c<data.getNCells(); c++)
	{
		// load the voxel data
		data.getCellRange(c, min, max);
		seeds.AddSeed(c, min, max);
	}
}

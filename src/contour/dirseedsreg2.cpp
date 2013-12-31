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
#include <stdio.h>
#include <contour/dirseedsreg2.h>
#include <contour/datareg2.h>
#include <contour/basic.h>

#define DEBUGNo
#define DEBUGSLEEPNo
#define sgn(x) ((x)>0 ? 1 : ((x)<0?-1:0))

extern int verbose;

void dirSeedsReg2::dirSweep(Datareg2& reg2)
{
	u_int i, j;
	Range resp;
	float min, max, t;
	float gradx;
	float grad1ya, grad1yb;
	float grad2ya, grad2yb;
	int keepflat;
	int prev;
	for(i=0; i<reg2.dim[0]-1; i++)
	{
		keepflat = 1;
		prev = -1;
		for(j=0; j<reg2.dim[1]-1; j++)
		{
			resp.MakeEmpty();
			// test responsiblity for each face
			// left
			if(i == 0)
			{
				min = max = reg2.getValue(reg2.index2vert(i,j));
				if((t=reg2.getValue(reg2.index2vert(i,j+1))) < min)
				{
					min = t;
				}
				if(t > max)
				{
					max = t;
				}
				if(min != max)
				{
					resp += Range(min,max);
				}
			}
			else
			{
				// never do anything other than on boundary
			}
			// right
			if(i == reg2.dim[0]-2)
			{
				// never keep a right boundary seed
			}
			else
			{
				// never keep a right boundary seed
			}
			// general case: bottom edge in middle
			// cell (i,j) and (i,j-1) share this x-grad
			gradx = reg2.getValue(reg2.index2vert(i+1,j)) -
					reg2.getValue(reg2.index2vert(i,j));
			// compute y-grad at (i,j) and (i+1,j)
			grad1ya = reg2.getValue(reg2.index2vert(i,j+1)) -
					  reg2.getValue(reg2.index2vert(i,j));
			grad1yb = reg2.getValue(reg2.index2vert(i+1,j+1)) -
					  reg2.getValue(reg2.index2vert(i+1,j));
			if(keepflat)
			{
				// check to see if gradient has 'turned'
				// only a seed if gradx & grady disagree in sign
				// note that 0 gradient is not considered opposite
				if(sgn(grad1ya) == 0 && sgn(grad1yb) == 0)
				{
					// flat cell (in y dim) - continue
				}
				else if(sgn(gradx) == -sgn(grad1ya) || sgn(gradx) == -sgn(grad1yb))
				{
					// extreme occurs if y components oppose each other
					// note that 0 gradient is not considered opposite
					min = max = reg2.getValue(reg2.index2vert(i,j));
					if((t=reg2.getValue(reg2.index2vert(i+1,j))) < min)
					{
						min = t;
					}
					if(t > max)
					{
						max = t;
					}
					resp += Range(min,max);
					keepflat = 0;
				}
			}
			// top
			if(j == reg2.dim[1]-2)
			{
				if(keepflat)
				{
					min = max = reg2.getValue(reg2.index2vert(i,j+1));
					if((t=reg2.getValue(reg2.index2vert(i+1,j+1))) < min)
					{
						min = t;
					}
					if(t > max)
					{
						max = t;
					}
					resp += Range(min,max);
				}
			}
			else
			{
				// only consider the top at the boundary
				if(!keepflat)
				{
					gradx = reg2.getValue(reg2.index2vert(i+1,j+1)) -
							reg2.getValue(reg2.index2vert(i,j+1));
					grad2ya = reg2.getValue(reg2.index2vert(i,j+1)) -
							  reg2.getValue(reg2.index2vert(i,j));
					grad2yb = reg2.getValue(reg2.index2vert(i+1,j+1)) -
							  reg2.getValue(reg2.index2vert(i+1,j));
					if(sgn(gradx) != 0 && (sgn(gradx) == sgn(grad2ya) || sgn(gradx) == sgn(grad2yb)))
					{
						keepflat=1;
					}
				}
				else
				{
					gradx = reg2.getValue(reg2.index2vert(i+1,j+1)) -
							reg2.getValue(reg2.index2vert(i,j+1));
					grad2ya = reg2.getValue(reg2.index2vert(i,j+1)) -
							  reg2.getValue(reg2.index2vert(i,j));
					grad2yb = reg2.getValue(reg2.index2vert(i+1,j+1)) -
							  reg2.getValue(reg2.index2vert(i+1,j));
					if(sgn(gradx) == -sgn(grad2ya) || sgn(gradx) == -sgn(grad2yb))
					{
						keepflat=0;
					}
				}
			}
			if(!resp.Empty())
			{
				if(prev == -1)
				{
					if(i!=0)
						prev = seeds.AddSeed(reg2.index2cell(i,j), resp.MinAll(),
											 resp.MaxAll());
					else
						seeds.AddSeed(reg2.index2cell(i,j), resp.MinAll(),
									  resp.MaxAll());
				}
				else
				{
					seeds.AddToRange(prev, resp.MinAll(), resp.MaxAll());
					prev = -1;
				}
			}
			else
			{
				prev = -1;
			}
		}
	}
}

void
dirSeedsReg2::compSeeds(void)
{
	if(verbose)
	{
		printf("------- computing seeds\n");
	}
	// clear the array of mark bits
	seeds.Clear();
	dirSweep((Datareg2&)data);
	if(verbose)
	{
		printf("computed %d seeds\n", seeds.getNCells());
	}
}

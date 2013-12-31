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
#include <contour/seeddirreg3.h>
#include <contour/datareg3.h>
#include <contour/basic.h>

#define sgn(x) ((x)>0 ? 1 : ((x)<0?-1:0))

void seedDirReg3::dirSweep(Datareg3& reg3)
{
	u_int i, j, k;
	Range resp;
	float min, max, t;
	float gradz;
	float grad1xa, grad1xb;
	float grad1ya, grad1yb;
	float grad2xa, grad2xb;
	float grad2ya, grad2yb;
	int keepflat, *keepflat_y;
	keepflat_y = (int*)malloc(sizeof(int)*reg3.dim[0]);
	for(k=0; k<reg3.dim[2]-1; k++)
	{
		for(i=0; i<reg3.dim[0]-1; i++)
		{
			keepflat_y[i] = 1;
		}
		for(j=0; j<reg3.dim[1]-1; j++)
		{
			keepflat = 1;
			for(i=0; i<reg3.dim[0]-1; i++)
			{
				resp.MakeEmpty();
				// test responsiblity for each face
				// minimum z
				if(k == 0)
				{
					min = max = reg3.getValue(i,j,k);
					if((t=reg3.getValue(i,j+1,k)) < min)
					{
						min = t;
					}
					if(t > max)
					{
						max = t;
					}
					if((t=reg3.getValue(i+1,j+1,k)) < min)
					{
						min = t;
					}
					if(t > max)
					{
						max = t;
					}
					if((t=reg3.getValue(i+1,j,k)) < min)
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
				// maximum z
				if(k == reg3.dim[2]-2)
				{
					// never keep a top boundary seed
				}
				else
				{
					// never keep a top boundary seed
				}
				// general case: bottom edge in middle
				// cell (i,j) and (i,j-1) share this x-grad
				gradz = reg3.getValue(i,j,k+1) - reg3.getValue(i,j,k);
				// compute grad at (i,j,k) and (i,j,k+1)
				grad1xa = reg3.getValue(i+1,j,k) -
						  reg3.getValue(i,j,k);
				grad1xb = reg3.getValue(i+1,j,k+1) -
						  reg3.getValue(i,j,k+1);
				grad1ya = reg3.getValue(i,j+1,k) -
						  reg3.getValue(i,j,k);
				grad1yb = reg3.getValue(i,j+1,k+1) -
						  reg3.getValue(i,j,k+1);
				if(keepflat && keepflat_y[i])
				{
					// check to see if gradient has 'turned'
					// only a seed if gradx & grady disagree in sign
					// note that 0 gradient is not considered opposite
					if(sgn(grad1xa) == 0 && sgn(grad1xb) == 0)
					{
						// flat cell (in x dim) - continue
					}
					else if(sgn(grad1ya) == 0 && sgn(grad1yb) == 0)
					{
						// flat cell (in y dim) - continue
					}
					else if((sgn(gradz) == -sgn(grad1xa) && sgn(gradz) == -sgn(grad1ya)) ||
							(sgn(gradz) == -sgn(grad1xb) && sgn(gradz) == -sgn(grad1yb)))
					{
						// extreme occurs if y components oppose each other
						// note that 0 gradient is not considered opposite
						min = max = reg3.getValue(i,j,k);
						if((t=reg3.getValue(i,j,k+1)) < min)
						{
							min = t;
						}
						if(t > max)
						{
							max = t;
						}
						resp += Range(min,max);
						keepflat = 0;
						keepflat_y[i] = 0;
					}
				}
				else
				{
				}
				// top
				if(i == reg3.dim[0]-2)
				{
					if(keepflat && keepflat_y[i])
					{
						// reached end at a flat.. add the edge values
						min = max = reg3.getValue(i+1,j,k);
						if((t=reg3.getValue(i+1,j,k+1)) < min)
						{
							min = t;
						}
						if(t > max)
						{
							max = t;
						}
						resp += Range(min,max);
					}
					if(j == reg3.dim[1]-2)
					{
						if(keepflat && keepflat_y[i])
						{
							// reached end at a flat.. add the edge values
							min = max = reg3.getValue(i+1,j+1,k);
							if((t=reg3.getValue(i+1,j+1,k+1)) < min)
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
						// do we need to set keepflat_y[i]?
						gradz = reg3.getValue(i,j+1,k+1) -
								reg3.getValue(i,j+1,k);
						grad2xa = reg3.getValue(i+1,j+1,k+1) - reg3.getValue(i,j+1,k+1);
						grad2xb = reg3.getValue(i+1,j+1,k)   - reg3.getValue(i,j+1,k);
						grad2ya = reg3.getValue(i,j+1,k+1) - reg3.getValue(i,j,k+1);
						grad2yb = reg3.getValue(i,j+1,k)   - reg3.getValue(i,j,k);
						keepflat_y[i] = (sgn(gradz) != 0 &&
										 ((sgn(gradz) == -sgn(grad2xa) && sgn(gradz) == sgn(grad2ya))
										  ||
										  (sgn(gradz) == -sgn(grad2xb) && sgn(gradz) == sgn(grad2yb))));
					}
				}
				else
				{
					if(j == reg3.dim[1]-2)
					{
						if(keepflat && keepflat_y[i])
						{
							// reached end at a flat.. add the edge values
							min = max = reg3.getValue(i+1,j+1,k);
							if((t=reg3.getValue(i+1,j+1,k+1)) < min)
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
						// do we need to set keepflat_y[i]?
						gradz = reg3.getValue(i,j+1,k+1) -
								reg3.getValue(i,j+1,k);
						grad2xa = reg3.getValue(i+1,j+1,k+1) - reg3.getValue(i,j+1,k+1);
						grad2xb = reg3.getValue(i+1,j+1,k)   - reg3.getValue(i,j+1,k);
						grad2ya = reg3.getValue(i,j+1,k+1) - reg3.getValue(i,j,k+1);
						grad2yb = reg3.getValue(i,j+1,k)   - reg3.getValue(i,j,k);
						keepflat_y[i] = (sgn(gradz) != 0 &&
										 ((sgn(gradz) == -sgn(grad2xa) && sgn(gradz) == sgn(grad2ya))
										  ||
										  (sgn(gradz) == -sgn(grad2xb) && sgn(gradz) == sgn(grad2yb))));
					}
					// do we need to set keepflat?
					gradz = reg3.getValue(i+1,j,k+1) -
							reg3.getValue(i+1,j,k);
					grad2xa = reg3.getValue(i+1,j,k+1)   - reg3.getValue(i,j,k+1);
					grad2xb = reg3.getValue(i+1,j,k)     - reg3.getValue(i,j,k);
					grad2ya = reg3.getValue(i+1,j+1,k+1) - reg3.getValue(i+1,j,k+1);
					grad2yb = reg3.getValue(i+1,j+1,k)   - reg3.getValue(i+1,j,k);
					keepflat = (sgn(gradz) != 0 &&
								((sgn(gradz) == sgn(grad2xa) && sgn(gradz) == -sgn(grad2ya))
								 ||
								 (sgn(gradz) == sgn(grad2xb) && sgn(gradz) == -sgn(grad2yb))));
				}
				if(!resp.Empty())
				{
					seeds.AddSeed(reg3.index2cell(i,j,k), resp.MinAll(), resp.MaxAll());
				}
			}
		}
	}
}

void seedDirReg3::compSeeds(void)
{
	// clear the array of mark bits
	seeds.Clear();
	dirSweep((Datareg3&)data);
}

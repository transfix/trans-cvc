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
#include <contour/dirseeds.h>
#include <contour/datavol.h>
#include <contour/dataslc.h>

#include <cstdio>
#include <cstring>
#include <cmath>

#define DEBUGNo
#define DEBUGSLEEPNo
#define sgn(x) ((x)>0 ? 1 : ((x)<0?-1:0))

extern int verbose;

void dirSeeds::dirSweep(Dataslc& slc)
{
	u_int c, f;
	Range resp;
	float g1[3], g2[3];
	float norm[2];
	int adjc;
	float min, max;
	for(c=0; c<slc.getNCells(); c++)
	{
		resp.MakeEmpty();
		slc.getCellGrad(c, g1);
		for(f=0; f<slc.getNCellFaces(); f++)
		{
			adjc = slc.getCellAdj(c,f);
			if(adjc != -1)
			{
				slc.normalToFace(c,f,norm);
				if(norm[1] >= 0.0)
				{
					slc.getCellGrad(adjc, g2);
					if(sgn(g1[0])==sgn(g1[1]) && g1[1] * g2[1] < 0.0)
					{
						slc.getFaceRange(c,f,min,max);
						resp += Range(min,max);
					}
				}
			}
			else
			{
				// boundary case... do something special?
				slc.normalToFace(c,f,norm);
				// first condition:  all left boundary cells are selected
				// second: top/bottom sides *may* be selected
				//         right hand cells should never be selected (sgn==0.0)
				if((std::fabs(norm[1]) < 0.0000001 && norm[0] < 0.0) ||
						//                (sgn(norm[1])*g1[0] > 0.0)) {
						(sgn(norm[1]) *(sgn(g1[0]*g1[1])) > 0.0))
				{
					slc.getFaceRange(c,f,min,max);
					resp += Range(min,max);
				}
			}
		}
		if(!resp.Empty())
		{
			seeds.AddSeed(c, resp.MinAll(), resp.MaxAll());
		}
	}
}

void dirSeeds::dirSweep(Datavol& vol)
{
	u_int c;
	for(c=0; c<vol.getNCells(); c++)
	{
	}
}

void dirSeeds::compSeeds(void)
{
	if(verbose)
	{
		printf("------- computing seeds\n");
	}
	// clear the array of mark bits
	seeds.Clear();
	dirSweep((Dataslc&)data);
	if(verbose)
	{
		printf("computed %d seeds\n", seeds.getNCells());
	}
}

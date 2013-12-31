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
#ifndef DIR_SEEDS_REG2_H
#define DIR_SEEDS_REG2_H

//#include <Utility/utility.h>
#include <contour/Conplot.h>
#include <contour/data.h>
#include <contour/range.h>
#include <contour/seedcells.h>

class Datareg2;

class dirSeedsReg2
{
	public:
		dirSeedsReg2(Data& d, SeedCells& s, Conplot& p) : data(d), seeds(s), plot(p) {}
		~dirSeedsReg2() {}
		void compSeeds(void);
	private:
		void dirSweep(Datareg2& reg);
		Data& data;
		SeedCells& seeds;
		Conplot&   plot;
};

#endif

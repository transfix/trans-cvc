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
#ifndef REG_PROP2_H
#define REG_PROP2_H

//#include <Utility/utility.h>
#include <contour/Conplot.h>
#include <contour/data.h>
#include <contour/range.h>
#include <contour/seedcells.h>

class Datavol;
class Dataslc;

class regProp2
{
	public:
		regProp2(Data& d, SeedCells& s, Conplot& p) : data(d), seeds(s), plot(p) {}
		~regProp2() {}
		void compSeeds(void);
	private:
		Data& data;
		SeedCells& seeds;
		Conplot&   plot;
};

#endif

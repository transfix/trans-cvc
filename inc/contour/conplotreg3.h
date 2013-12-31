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
// conPlot2d.h - class for preprocessing and extraction of surfaces from 3d data

#ifndef CONPLOT_RECT_3D_H
#define CONPLOT_RECT_3D_H

//#include <Utility/utility.h>
#include <contour/CellQueue.h>
#include <contour/Conplot.h>
#include <contour/contour3d.h>
#include <contour/Dataset.h>
#include <contour/datasetreg3.h>
#include <contour/edgehash.h>
#include <contour/range.h>
#include <contour/seedcells.h>
#include <contour/segtree.h>

class Conplotreg3 : public Conplot
{
	public:
		Conplotreg3(Datasetreg3* d);
		~Conplotreg3();
		Contour3d* con3, *curcon;

	protected:
		// extract in 3d (from memory) or slice-by-slice (swap from disk)
		u_int ExtractAll(float isovalue);
		void interpRect3Dpts_x(int, int, int, float*, float*, float [3][8], int, int, float, float*, float*, float*);
		void interpRect3Dpts_y(int, int, int, float*, float*, float [3][8], int, int, float, float*, float*, float*);
		void interpRect3Dpts_z(int, int, int, float*, float*, float [3][8], int, int, float, float*, float*, float*);
		int InterpEdge(float*, float*, float [3][8], float, int, int, int, int);
		// track a contour from a seed cell
		void TrackContour(float, int);
		// enqueue faces for propagation of surface
		inline void EnqueueFaces(int, u_int, u_int, u_int, CellQueue&);
		void Reset(int t)
		{
			con3[t].Reset();
		}
		int  Size(int t)
		{
			return(con3[t].getSize());
		}
		int  isDone(int t)
		{
			return(con3[t].isDone());
		}
		void Done(int t)
		{
			con3[t].Done();
		}

	private:

		Datasetreg3*	reg3;
		Datareg3*	curreg3;
};

#endif

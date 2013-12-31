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
// contour2d.h - Class for a 2d isocontour polyline
// contour2d is a class for representing a 2d isocontour polyline.

#ifndef _CONTOUR_2D_H
#define _CONTOUR_2D_H

//#include <Utility/utility.h>
#include <contour/basic.h>

#include <cstring>

class Contour2d
{
	public:
		// arrays of vertices, and edges
		float(*vert)[2];			// polyline vertex array
		unsigned int (*edge)[2];		// array of polyline edges

		// constructor
		Contour2d();

		// destructor
		~Contour2d();

		// reset (delete all vertices and triangles)
		void Reset(void);
		void Done(void);
		int  isDone(void)
		{
			return(done);
		}

		// add a vertex with the given position and normal
		int AddVert(float p[2])
		{
			return(AddVert(p[0], p[1]));
		}
		int AddVert(float, float);

		// add an edge indexed by the given 2 vertices
		int AddEdge(u_int v[2])
		{
			return(AddEdge(v[0], v[1]));
		}
		int AddEdge(u_int, u_int);

		// get the number of vertices or edges
		int getSize(void)
		{
			return(nedge);
		}
		int getNVert(void)
		{
			return(nvert);
		}
		int getNEdge(void)
		{
			return(nedge);
		}

		// write vertices and triangles to a file
		int write(char* filename);

		void setExtent(float min[3], float max[3])
		{
			std::memcpy(minext, min, sizeof(float[3]));
			std::memcpy(maxext, max, sizeof(float[3]));
		}

	protected :

		int done; // done with isocontour ?
		// the size of the vertex and edge arrays
		int vsize, tsize;
		// the number of vertices and edges
		int nvert, nedge;
		float minext[3], maxext[3];

};

#endif

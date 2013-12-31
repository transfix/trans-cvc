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
// contour3d is a class for representing an isocontour surface (or any 3d triangular mesh).

#ifndef _CONTOUR_3D_H
#define _CONTOUR_3D_H

#ifdef USEDICT
#include <contour/dict.h>
#endif

//#include <Utility/utility.h>
#include <contour/basic.h>

#include <cstring>

class Contour3d
{

	public:
		Contour3d(int fn=0);
		~Contour3d();

		// to color isosurface using another variable (for multivariate data)
		void colorByFun(int c)
		{
			color =  c;
		}
		void minmaxFun(float mn, float mx)
		{
			fmin=mn;
			fmax=mx;
		}

		// reset (delete all vertices and triangles)
		void Reset(void);
		void Done(void);
		int  isDone(void)
		{
			return(done);
		}

		// add a vertex with the given position and normal
		int AddVert(float p[3], float n[3], float f=0.0)
		{
			return(AddVert(p[0], p[1], p[2], n[0], n[1], n[2], f));
		}
		int AddVert(float, float, float, float, float, float, float f=0.0);

		int AddVertUnique(float p[3], float n[3], float f=0.0)
		{
			return(AddVertUnique(p[0], p[1], p[2], n[0], n[1], n[2], f));
		}
		int AddVertUnique(float, float, float, float, float, float, float f=0.0);

		// add a triangle indexed by the given 3 vertices
		int AddTri(u_int v[3])
		{
			return(AddTri(v[0], v[1], v[2]));
		}
		int AddTri(u_int, u_int, u_int);

		// get the number of vertices or triangles
		int getSize(void)
		{
			return(ntri);
		}
		int getNVert(void)
		{
			return(nvert);
		}
		int getNTri(void)
		{
			return(ntri);
		}

		// write vertices and triangles to a file
		int write(char* filename);

		void setExtent(float min[3], float max[3])
		{
			std::memcpy(minext, min, sizeof(float[3]));
			std::memcpy(maxext, max, sizeof(float[3]));
		}

	protected:

		int	done;	// done with isocontour?

		// the size of the vertex and triangle arrays
		int	vsize, tsize;

		// the number of vertices and triangles
		int	nvert, ntri;

		float	minext[3], maxext[3];

		/* Added by Joe R - 12/21/2005, used for fast searching of duplicate vertices */
#ifdef USEDICT
		dict_t    vertex_dict;
#endif

	public :				// made public by Emilio

		// true if colored by function on contour
		int	color;
		int	vf;				// variable used for coloring
		float	fmin, fmax;			// min and max color values

		// arrays of vertices, vertex normals, and triangles
		float(*vert)[3];			// isosurface vertex array
		float(*vnorm)[3];			// array of vertex normals
		float(*vfun);			// color values at vertices

		unsigned int (*tri)[3];			// triangle mesh array
};

#endif

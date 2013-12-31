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
//      Contour2d - Class for a 2d contour curve
//      Copyright (c) 1998 Emilio Camahort, Dan Schikore

#include <stdio.h>
#include <contour/contour2d.h>

#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif


// Contour2d() - basic constructor
Contour2d::Contour2d()
{
	done = 0;
	nvert = 0;
	nedge = 0;
	vsize =  500;
	tsize = 1000;
	vert  = (float(*)[2])malloc(sizeof(float[2]) * vsize);
	edge  = (u_int(*)[2])malloc(sizeof(u_int[2]) * tsize);
}

// ~Contour2d() - free allocated memory
Contour2d::~Contour2d()
{
	free(vert);
	free(edge);
}

// AddVert() - add a vertex with the given (unit) normal
int Contour2d::AddVert(float x, float y)
{
	int n = nvert++;
	if(nvert > vsize)
	{
		vsize<<=1;
		vert  = (float(*)[2])realloc(vert, sizeof(float[2]) * vsize);
	}
	vert[n][0] = x;
	vert[n][1] = y;
	return(n);
}

// AddEdge() - add an edge indexed by its 2 vertices
int Contour2d::AddEdge(u_int v1, u_int v2)
{
	int n = nedge++;
	if(nedge > tsize)
	{
		tsize<<=1;
		edge = (u_int(*)[2])realloc(edge, sizeof(u_int[2]) * tsize);
	}
	edge[n][0] = v1;
	edge[n][1] = v2;
	return(n);
}

// Reset() - clear vertex and edge info
void Contour2d::Reset(void)
{
	nvert = 0;
	nedge  = 0;
	done = 0;
}

void Contour2d::Done(void)
{
	done = 1;
}

// write() - write vertex and triangles to a file

int Contour2d::write(char* filename)
{
	FILE* fp;
	int v, t;
	fp = fopen(filename, "w");
	// silent failure --> changed by Emilio: return 1 = ERROR
	if(fp == NULL)
	{
		return 1;
	}
	fprintf(fp, "%d %d 0 0 0 0 0\n0 0 0\n", nvert, nedge);
	// Emilio: following construct gives an out-of-bounds warning, I'll just
	//	      write a zero instead for the time being, since writing nothing
	//	      would screw up the file format
#ifdef EMILIO
	for(v=0; v<nvert; v++)
	{
		fprintf(fp, "%g %g %g\n", vert[v][0], vert[v][1], vert[v][2]);
	}
#endif /* of EMILIO */
	for(v=0; v<nvert; v++)
	{
		fprintf(fp, "%g %g %g\n", vert[v][0], vert[v][1], 0.0);
	}
	fprintf(fp, "0 0\n");
	for(t=0; t<nedge; t++)
	{
		fprintf(fp, "%d %d\n", edge[t][0], edge[t][1]);
	}
	fclose(fp);
	return 0;
}

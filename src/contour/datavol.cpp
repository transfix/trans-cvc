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
// volume.C - class for a regular volume of scalar data
// Copyright (c) 1997 Dan Schikore - updated by Emilio Camahort, 1999

#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif

#include <stdio.h>
#include <string.h>
#include <contour/datavol.h>
#include <contour/compute.h>
#include <contour/endian_io.h>

#include <cmath>

#define FSAMPLES 256

extern int verbose;

// Datavol() - alternative constructor for the libcontour library
Datavol::Datavol(Data::DataType t, u_int ndata, u_int nverts, u_int ncells,
				 double* _verts, u_int* _cells, int* _celladj, u_char* data)
	: Data(t, ndata, data)
{
	u_int i;
	Datavol::nverts = nverts;			// initializations
	Datavol::ncells = ncells;
	verts = (float(*)[3])_verts;
	cells = (u_int(*)[4])_cells;
	celladj = (int (*)[4])_celladj;
	if(verbose)
	{
		printf("computing extent\n");    // compute data extents
	}
	minext[0] = minext[1] = minext[2] = 1e10;
	maxext[0] = maxext[1] = maxext[2] = -1e10;
	for(i = 0; i < nverts; i++)
	{
		if(verts[i][0] < minext[0])
		{
			minext[0] = verts[i][0];
		}
		if(verts[i][0] > maxext[0])
		{
			maxext[0] = verts[i][0];
		}
		if(verts[i][1] < minext[1])
		{
			minext[1] = verts[i][1];
		}
		if(verts[i][1] > maxext[1])
		{
			maxext[1] = verts[i][1];
		}
		if(verts[i][2] < minext[2])
		{
			minext[2] = verts[i][2];
		}
		if(verts[i][2] > maxext[2])
		{
			maxext[2] = verts[i][2];
		}
	}
	if(verbose)
		printf("  min = %f %f %f  max = %f %f %f\n",
			   minext[0], minext[1], minext[2],
			   maxext[0], maxext[1], maxext[2]);
	if(verbose)
	{
		printf("%d verts, %d cells\n", Datavol::nverts, Datavol::ncells);
	}
	// compute gradients
	grad = (float(*)[3])malloc(sizeof(float[3])*getNVerts());
	for(i=0; i<getNCells(); i++)
	{
		if(cells[i][0] == 100 ||
				cells[i][1] == 100 ||
				cells[i][2] == 100 ||
				cells[i][3] == 100)
		{
			if(verbose)
			{
				printf("%d %d %d %d\n", cells[i][0], cells[i][1], cells[i][2], cells[i][3]);
			}
		}
		if(cells[i][0] == 101 ||
				cells[i][1] == 101 ||
				cells[i][2] == 101 ||
				cells[i][3] == 101)
		{
			if(verbose)
			{
				printf("%d %d %d %d\n", cells[i][0], cells[i][1], cells[i][2], cells[i][3]);
			}
		}
		if(verbose > 1)
			printf("cell %d: %d %d %d %d (%d %d %d %d)\n", i,
				   cells[i][0],
				   cells[i][1],
				   cells[i][2],
				   cells[i][3],
				   celladj[i][0],
				   celladj[i][1],
				   celladj[i][2],
				   celladj[i][3]);
	}
	for(i=0; i<getNCells(); i++)
	{
		for(u_int j=0; j<getNCellFaces(); j++)
		{
			int adj = celladj[i][j];
			int same = 0;
			if(adj != -1)
			{
				for(int k=0; k<4; k++)
					for(int l=0; l<4; l++)
						if(cells[i][k] == cells[adj][l])
						{
							same++;
						}
				if(verbose)
					if(same != 3)
						printf("cell %d (%d %d %d %d) not adj to %d (%d %d %d %d)\n",
							   i, cells[i][0], cells[i][1], cells[i][2], cells[i][3],
							   adj, cells[adj][0], cells[adj][1], cells[adj][2], cells[adj][3]);
			}
		}
	}
	preprocessData(data);
	compGrad();
}

//  compLength() -
float* Datavol::compLength(u_int& len, float** funx)
{
	float* val = (float*)malloc(sizeof(float)*FSAMPLES);
	float* fx = (float*)malloc(sizeof(float)*FSAMPLES);
	u_int c;
	u_int* v;
	len = FSAMPLES;
	memset(val, 0, sizeof(float)*len);
	*funx = fx;
	for(c=0; c<len; c++)
	{
		fx[c] = getMin() + (c/(len-1.0f)) * (getMax()-getMin());
	}
	for(c=0; c<getNCells(); c++)
	{
		v = getCellVerts(c);
		tetSurfIntegral(getVert(v[0]), getVert(v[1]), getVert(v[2]), getVert(v[3]),
						getValue(v[0]), getValue(v[1]), getValue(v[2]), getValue(v[3]),
						fx, val, len, getMin(), getMax(), 1.0);
	}
	return(val);
}

//  compGradient() -
float* Datavol::compGradient(u_int& len, float** funx)
{
	float* val = (float*)malloc(sizeof(float)*FSAMPLES);
	float* fx = (float*)malloc(sizeof(float)*FSAMPLES);
	float cellgrad[4], scaling;
	u_int c;
	u_int* v;
	len = FSAMPLES;
	memset(val, 0, sizeof(float)*len);
	*funx = fx;
	for(c=0; c<len; c++)
	{
		fx[c] = getMin() + (c/(len-1.0f)) * (getMax()-getMin());
	}
	for(c=0; c<getNCells(); c++)
	{
		v = getCellVerts(c);
		getCellGrad4(c, cellgrad);
		scaling = (sqr(cellgrad[0]) + sqr(cellgrad[1]) + sqr(cellgrad[2])) /
				  sqr(cellgrad[3]);
		tetSurfIntegral(getVert(v[0]), getVert(v[1]), getVert(v[2]), getVert(v[3]),
						getValue(v[0]), getValue(v[1]), getValue(v[2]), getValue(v[3]),
						fx, val, len, getMin(), getMax(), (float)std::fabs(scaling));
	}
	return(val);
}

//  compArea() -
float* Datavol::compArea(u_int& len, float** funx)
{
	float* val = (float*)malloc(sizeof(float)*FSAMPLES);
	float* cum = (float*)malloc(sizeof(float)*FSAMPLES);
	float* fx = (float*)malloc(sizeof(float)*FSAMPLES);
	float sum;
	u_int c;
	u_int* v;
	len = FSAMPLES;
	memset(val, 0, sizeof(float)*len);
	memset(cum, 0, sizeof(float)*len);
	*funx = fx;
	for(c=0; c<len; c++)
	{
		fx[c] = getMin() + (c/(len-1.0f)) * (getMax()-getMin());
	}
	for(c=0; c<getNCells(); c++)
	{
		v = getCellVerts(c);
		tetVolIntegral(getVert(v[0]), getVert(v[1]), getVert(v[2]), getVert(v[3]),
					   getValue(v[0]), getValue(v[1]), getValue(v[2]), getValue(v[3]),
					   fx, val, cum, len, getMin(), getMax(), 1.0);
	}
	// sum the results to add all
	sum=0;
	for(c=0; c<len; c++)
	{
		val[c] += sum;
		sum+=cum[c];
	}
	return(val);
}

//  compMaxArea() -
float* Datavol::compMaxArea(u_int& len, float** funx)
{
	float* val;
	float max;
	u_int i;
	val = compArea(len, funx);
	max = val[len-1];
	for(i=0; i<len; i++)
	{
		val[i] = max-val[i];
	}
	return(val);
}

//  compFunction(), fName() -
float* Datavol::compFunction(int n, u_int& len, float** fx)
{
	switch(n)
	{
		case 0:
			return(compLength(len, fx));
		case 1:
			return(compGradient(len, fx));
		case 2:
			return(compArea(len, fx));
		case 3:
			return(compMaxArea(len, fx));
	}
	return(NULL);
}

char* Datavol::fName(int n)
{
	switch(n)
	{
		case 0:
			return((char*)"Surface Area");
		case 1:
			return((char*)"Gradient");
		case 2:
			return((char*)"Min Volume");
		case 3:
			return((char*)"Max Volume");
	}
	return(NULL);
}

// compGrad() - compute gradients
void Datavol::compGrad(void)
{
	float u[4], v[4], w[4], x[4], y[4], z[4], g[4];
	int v0, v1, v2, v3;
	float len;
	float weight;
	memset(grad, 0, sizeof(float[3])*getNVerts());
	for(int i=0; i<getNCells(); i++)
	{
		if(verbose > 1)
		{
			printf("grad for cell %d\n", i);
		}
		v0 = cells[i][0];
		v1 = cells[i][1];
		v2 = cells[i][2];
		v3 = cells[i][3];
		u[0] = verts[v1][0] - verts[v0][0];
		u[1] = verts[v1][1] - verts[v0][1];
		u[2] = verts[v1][2] - verts[v0][2];
		u[3] = getValue(v1) - getValue(v0);
		v[0] = verts[v2][0] - verts[v0][0];
		v[1] = verts[v2][1] - verts[v0][1];
		v[2] = verts[v2][2] - verts[v0][2];
		v[3] = getValue(v2) - getValue(v0);
		w[0] = verts[v3][0] - verts[v0][0];
		w[1] = verts[v3][1] - verts[v0][1];
		w[2] = verts[v3][2] - verts[v0][2];
		w[3] = getValue(v3) - getValue(v0);
		x[0] = verts[v3][0] - verts[v1][0];
		x[1] = verts[v3][1] - verts[v1][1];
		x[2] = verts[v3][2] - verts[v1][2];
		y[0] = verts[v3][0] - verts[v2][0];
		y[1] = verts[v3][1] - verts[v2][1];
		y[2] = verts[v3][2] - verts[v2][2];
		z[0] = verts[v2][0] - verts[v1][0];
		z[1] = verts[v2][1] - verts[v1][1];
		z[2] = verts[v2][2] - verts[v1][2];
		g[0] =   u[1] * (v[2]*w[3] - v[3]*w[2])
				 + u[2] * (v[3]*w[1] - v[1]*w[3])
				 + u[3] * (v[1]*w[2] - v[2]*w[1]);
		g[1] =   u[0] * (v[2]*w[3] - v[3]*w[2])
				 + u[2] * (v[3]*w[0] - v[0]*w[3])
				 + u[3] * (v[0]*w[2] - v[2]*w[0]);
		g[2] =   u[0] * (v[1]*w[3] - v[3]*w[1])
				 + u[1] * (v[3]*w[0] - v[0]*w[3])
				 + u[3] * (v[0]*w[1] - v[1]*w[0]);
		if(verbose > 1)
		{
			printf(" grad %f %f %f\n", g[0], g[1], g[2]);
		}
		if(verbose)
			if(v0 == 101 || v1 == 101 || v2 == 101 || v3 == 101)
			{
				printf("v100: %f %f %f\n", g[0], g[1], g[2]);
			}
		weight =   u[0]*(v[1]*w[2] - v[2]*w[1])
				   - u[1]*(v[2]*w[0] - v[0]*w[2])
				   + u[2]*(v[0]*w[1] - v[1]*w[0]);
		weight=1;
		grad[v0][0] += g[0] * weight;
		grad[v0][1] += g[1] * weight;
		grad[v0][2] += g[2] * weight;
		weight = (-u[0])*(x[1]*z[2] - x[2]*z[1])
				 - (-u[1])*(x[2]*z[0] - x[0]*z[2])
				 + (-u[2])*(x[0]*z[1] - x[1]*z[0]);
		weight=1;
		grad[v1][0] += g[0] * weight;
		grad[v1][1] += g[1] * weight;
		grad[v1][2] += g[2] * weight;
		weight = (-v[0])*((-z[1])*y[2] - (-z[2])*y[1])
				 - (-v[1])*((-z[2])*y[0] - (-z[0])*y[2])
				 + (-v[2])*((-z[0])*y[1] - (-z[1])*y[0]);
		weight=1;
		grad[v2][0] += g[0] * weight;
		grad[v2][1] += g[1] * weight;
		grad[v2][2] += g[2] * weight;
		weight = (-w[0])*((-y[1])*(-x[2]) - (-y[2])*(-x[1]))
				 - (-w[1])*((-y[2])*(-x[0]) - (-y[0])*(-x[2]))
				 + (-w[2])*((-y[0])*(-x[1]) - (-y[1])*(-x[0]));
		weight=1;
		grad[v3][0] += g[0] * weight;
		grad[v3][1] += g[1] * weight;
		grad[v3][2] += g[2] * weight;
	}
	for(int i=0; i<getNVerts(); i++)
	{
		if(verbose > 1)
		{
			printf("scaling vgrad %d\n", i);
		}
		len = (float)(sqrt(grad[i][0]*grad[i][0] + grad[i][1]*grad[i][1] +
						   grad[i][2]*grad[i][2]));
		if(len != 0.0)
		{
			grad[i][0] /= len;
			grad[i][1] /= len;
			grad[i][2] /= len;
		}
	}
	if(verbose)
	{
		printf("grad101 = %f %f %f\n", grad[101][0], grad[101][1], grad[101][2]);
	}
}

// ~Datavol() - destroy a volume
Datavol::~Datavol()
{
	if(filename)
	{
		free(verts);
		free(cells);
		free(celladj);
	}
}

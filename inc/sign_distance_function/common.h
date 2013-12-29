/*
  Copyright 2004-2005 The University of Texas at Austin

        Authors: Lalit Karlapalem <ckl@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of sign_distance_function.

  sign_distance_function is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  sign_distance_function is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef CCV_SDF_COMMON_H
#define CCV_SDF_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <map>

#include <sign_distance_function/head.h>

namespace SDFLibrary {


bool initSDF();

void readGeom(int nverts, float* verts, int ntris, int* tris);

void adjustData();

void compute();



int isEqual (double one, double two);

int isZero(double num);

int isNegative(double num);

int isBetween(double one, double two, double num);

int isZero(SDFLibrary::myPoint one);

int isSame(SDFLibrary::myPoint one, SDFLibrary::myPoint two);

void init_all_vars();



void propagate_left(int i, int j, int k);

void propagate_bottom(int i, int j, int k);

void propagate_inside(int i, int j, int k);

void propagate_right(int i, int j, int k);

void propagate_top(int i, int j, int k);

void propagate_outside(int i, int j, int k);

void apply_distance_transform(int vi, int vj, int vk);

void insert_bound_vert(int vert);




int index2vert(int i, int j, int k);

void _vert2index(int c, int &i, int &j, int &k);

int index2cell(int i, int j, int k);

void _cell2index(int c, int &i, int &j, int &k);

double xCoord(int i);

double yCoord(int i);

double zCoord(int i);

void object2octree(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int &ci, int &cj, int &ck);

double getTime();

	extern double MAX_DIST;
	extern int size;

	extern triangle* surface;
	extern myVert* vertices;
	extern myPoint* normals;
	extern double* distances;
	extern cell*** sdf;
	extern voxel* values;
	extern int total_points, total_triangles, all_verts_touched;
	extern double minx, miny, minz, maxx, maxy, maxz;

	extern double TOLERANCE;

	extern int octree_depth;
	extern int flipNormals;

	extern bool *bverts;
	extern int *queues;

	extern double minext[3];
	extern double maxext[3];
	extern double span[3];
	
}; //namespace SDFLibrary

#endif

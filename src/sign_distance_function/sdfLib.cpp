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

#include <algorithm>
#include <boost/scoped_array.hpp>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sign_distance_function/common.h>
#include <sign_distance_function/sdfLib.h>

using namespace SDFLibrary;

void free_memory()
{
  int i, j, k;
  SDFLibrary::listnode* temp;
  SDFLibrary::listnode* currNode;

  printf("starting memory de-allocation\n");

  //1. Octree
  for (i = 0; i < SDFLibrary::size; i++)
    {
      for (j = 0; j < SDFLibrary::size; j++)
        {
          for (k = 0; k < SDFLibrary::size; k++)
            {
              currNode = SDFLibrary::sdf[i][j][k].tindex;

              while(currNode != NULL)
                {
                  temp = currNode;
                  currNode = currNode->next;
                  free(temp);
                }
            }
          free(SDFLibrary::sdf[i][j]);
        }
      free(SDFLibrary::sdf[i]);
    }   
  free(SDFLibrary::sdf);

  free(SDFLibrary::values);

  if (SDFLibrary::vertices != NULL)
    free(SDFLibrary::vertices);

  if (SDFLibrary::surface != NULL)
    free(SDFLibrary::surface);

  if (SDFLibrary::normals != NULL)
    free(SDFLibrary::normals);

  if (SDFLibrary::distances != NULL)
    free(SDFLibrary::distances);

  if (SDFLibrary::queues != NULL)
    free(SDFLibrary::queues);

  if (SDFLibrary::bverts != NULL)
    free(SDFLibrary::bverts);

  printf("Memory de-allocated successfully! \n");
}

void SDFLibrary::setParameters(int Size, int isNormalFlip, float* mins, float* maxs)
{
  //First the default values.
  SDFLibrary::init_all_vars();
        
  //Then, assign the actual input values.
  SDFLibrary::size = Size;
  SDFLibrary::flipNormals = isNormalFlip;

  SDFLibrary::minext[0] = mins[0];      SDFLibrary::minext[1] = mins[1];        SDFLibrary::minext[2] = mins[2];
  SDFLibrary::maxext[0] = maxs[0];      SDFLibrary::maxext[1] = maxs[1];        SDFLibrary::maxext[2] = maxs[2];
  SDFLibrary::span[0] = (maxs[0]-mins[0])/(SDFLibrary::size);
  SDFLibrary::span[1] = (maxs[1]-mins[1])/(SDFLibrary::size);
  SDFLibrary::span[2] = (maxs[2]-mins[2])/(SDFLibrary::size);

  if ((Size!=16) && (Size!=32) &&(Size!=64) && (Size!=128) && (Size!=256) &&(Size!=512) &&(Size!=1024))
    {
      printf("size is incorrect\n");
      exit(1);
    }
}

float* SDFLibrary::computeSDF(int nverts, float* verts, int ntris, int* tris)
{
  int i, numb;
  float* sdfValues =NULL;
  float isoval;

  //Set up the volume grid
  if( !initSDF() ) return 0;

  //Read in the Geometry
  readGeom(nverts, verts, ntris, tris);

  //Setup the Octree
  adjustData();

  //Compute the SDF
  compute();

  //Return the SDF
  numb = (SDFLibrary::size+1)*(SDFLibrary::size+1)*(SDFLibrary::size+1);
  sdfValues = (float*)(malloc(sizeof(float)*(numb)));
  isoval = 100.0f;

  for (i=0; i<numb; i++)
    sdfValues[i] = SDFLibrary::values[i].value * SDFLibrary::values[i].signe;

  free_memory();

  return (sdfValues);
}       

RAWIV_header* SDFLibrary::getVolumeInfo()
{
  int i;

  RAWIV_header* ret = (RAWIV_header*)(malloc(sizeof(RAWIV_header)*1));

  for (i=0; i<3; i++)
    {
      ret->minext[i] = SDFLibrary::minext[i];
      ret->maxext[i] = SDFLibrary::maxext[i];
      ret->span[i] = SDFLibrary::span[i];
      ret->origin[i] = 0.0f;
      ret->dim[i] = SDFLibrary::size+1;
    }

  ret->ngridpts = (SDFLibrary::size+1)*(SDFLibrary::size+1)*(SDFLibrary::size+1);
  ret->ncells = (SDFLibrary::size)*(SDFLibrary::size)*(SDFLibrary::size);
  ret->size = SDFLibrary::size;

  return ret;
}


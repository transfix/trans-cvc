/*
  Copyright 2008 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeRover.

  VolumeRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <cvc/algorithm.h>
#include <cvc/utility.h>

#ifdef CVC_USING_MULTI_SDF
#include <multi_sdf/multi_sdf.h>
#endif

#ifdef CVC_USING_SDFLIB
#include <sign_distance_function/sdfLib.h>
#endif

#include <iostream>
#include <cstring>
#include <algorithm>

namespace
{
  CVC_DEF_EXCEPTION(sign_distance_function_error);

#ifdef CVC_USING_SDFLIB
  CVC_NAMESPACE::volume sdf_library(const CVC_NAMESPACE::geometry& geom,
				    const CVC_NAMESPACE::dimension& dim,
				    const CVC_NAMESPACE::bounding_box& bbox)
  {
    using namespace CVC_NAMESPACE;

    int size;
    int flipNormals = 0;

    //This library is limited to volumes with dimx == dimy == dimz and dimx == 2^n up to 1024...
    //but lets fake support for arbitrary dimensions by using VolMagick::Volume::resize()
  
    unsigned int maxdim = std::max(dim[0],std::max(dim[1],dim[2]));
    maxdim = upToPowerOfTwo(maxdim);
    size = std::min((unsigned int)1024,maxdim);

    float mins[3], maxs[3];
    point_t geom_min = geom.min_point();
    point_t geom_max = geom.max_point();

    if(bbox.isNull())
      {
	for(int i = 0; i < 3; i++) mins[i] = geom_min[i];
	for(int i = 0; i < 3; i++) maxs[i] = geom_max[i];
      }
    else
      {
	mins[0] = bbox[0];
	mins[1] = bbox[1];
	mins[2] = bbox[2];
	maxs[0] = bbox[3];
	maxs[1] = bbox[4];
	maxs[2] = bbox[5];
      }

    SDFLibrary::setParameters(size,flipNormals,mins,maxs);

    //TODO: avoid this copy
    int nverts = geom.num_points();
    boost::scoped_array<float> verts(new float[nverts*3]);
    for(int i = 0; i < nverts; i++)
      {
	verts[i*3+0] = geom.const_points()[i][0];
	verts[i*3+1] = geom.const_points()[i][1];
	verts[i*3+2] = geom.const_points()[i][2];
      }

    int ntris = geom.num_tris();
    boost::scoped_array<int> tris(new int[ntris*3]);
    for(int i = 0; i < ntris; i++)
      {
	tris[i*3+0] = geom.const_tris()[i][0];
	tris[i*3+1] = geom.const_tris()[i][1];
	tris[i*3+2] = geom.const_tris()[i][2];
      }

    float* values = SDFLibrary::computeSDF(geom.num_points(),
					   verts.get(),
					   ntris,
					   tris.get());
    if(!values)
      throw sign_distance_function_error("SDFLibrary::computeSDF() failed.\n");

    boost::scoped_array<float> choppedValues(new float[size*size*size]);
    {
      int i, j, k;
      int c=0;
      for( i=0; i<=size; i++ )
	{
	  for( j=0; j<=size; j++ )
	    {
	      for( k=0; k<=size; k++ )
		{
		  if( i!=size && j!=size && k!=size )
		    choppedValues[c++] = values[i*(size+1)*(size+1) + j*(size+1) + k];
		}
	    }
	}
    }
    delete [] values;

    volume vol(dimension(size,size,size),
	       Float,
	       bounding_box(mins[0],mins[1],mins[2],
			    maxs[0],maxs[1],maxs[2]));
    memcpy(*vol,choppedValues.get(),size*size*size*sizeof(float));
    vol.resize(dim);
    vol.desc("sign_distance_function");
    return vol;
  }
#endif
}

namespace CVC_NAMESPACE
{
  volume sdf(const geometry& geom,
	     const dimension& dim,
	     const bounding_box& bbox,
	     sdf_method method)
  {
    volume vol;

    switch(method)
      {
      case MULTI_SDF:
#ifdef CVC_USING_MULTI_SDF
        vol = multi_sdf::signedDistanceFunction(boost::shared_ptr<Geometry>(new Geometry(convert(geom))),dim,bbox);
        vol.desc("Signed Distance Function - multi_sdf");
#else
        throw unsupported_exception("multi_sdf unsupported");
#endif
        break;
      case SDFLIB:
#ifdef CVC_USING_SDFLIB
        vol = sdf_library(geom,dim,bbox);
        vol.desc("Signed Distance Function - SDFLibrary");
#else
        throw unsupported_exception("SDFLibrary unsupported");
#endif
        break;
      }

    return vol;
  }
}
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
#include <cvc/app.h>

#ifdef CVC_USING_MULTI_SDF
#include <multi_sdf/mesh_io.h>
#include <multi_sdf/sdf.h>
#include <multi_sdf/kdtree.h>
#include <multi_sdf/matrix.h>
#include <multi_sdf/dt.h>
#endif

#ifdef CVC_USING_SDFLIB
#include <sign_distance_function/sdfLib.h>
#endif

//fast_contouring is broken right now... :(
#include <fast_contouring/FastContouring.h>

//libcontour headers
#include <contour/contour.h>
#include <contour/datasetreg3.h>

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
    thread_info ti(BOOST_CURRENT_FUNCTION);

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

#ifdef CVC_USING_MULTI_SDF
  CVC_NAMESPACE::volume multi_sdf_library(const CVC_NAMESPACE::geometry& geom,
					  const CVC_NAMESPACE::dimension& dim,
					  const CVC_NAMESPACE::bounding_box& bbox)
  {
    using namespace CVC_NAMESPACE;
    using namespace multi_sdf;
    thread_info ti(BOOST_CURRENT_FUNCTION);

    int dimx = dim[0], dimy = dim[1], dimz = dim[2];

    Mesh mesh;
    cvcapp.log(3,"Reading input mesh");

    //TODO: avoid this copy
    int nverts = geom.num_points();
    boost::scoped_array<float> verts(new float[nverts*3]);
    for(int i = 0; i < nverts; i++)
      {
	verts[i*3+0] = geom.const_points()[i][0];
	verts[i*3+1] = geom.const_points()[i][1];
	verts[i*3+2] = geom.const_points()[i][2];
      }

    boost::scoped_array<float> colors(new float[nverts*3]);
    if(geom.const_colors().empty())
      for(int i = 0; i < nverts*3; i++)
	colors[i] = 1.0f;
    else
      for(int i = 0; i < nverts; i++)
	{
	  colors[i*3+0] = geom.const_colors()[i][0];
	  colors[i*3+1] = geom.const_colors()[i][1];
	  colors[i*3+2] = geom.const_colors()[i][2];
	}

    int ntris = geom.num_tris();
    boost::scoped_array<int> tris(new int[ntris*3]);
    for(int i = 0; i < ntris; i++)
      {
	tris[i*3+0] = geom.const_tris()[i][0];
	tris[i*3+1] = geom.const_tris()[i][1];
	tris[i*3+2] = geom.const_tris()[i][2];
      }

    read_labeled_mesh(mesh,
		      nverts, verts.get(), colors.get(), ntris, tris.get());
    cvcapp.log(3,"done.");

    // build a bounding box around the input and store the
    // origin, span etc.
    //  vector<double> bbox;
    //  construct_bbox(mesh, bbox);
    bounding_box box(bbox);
    if(box.isNull())
      {
	float mins[3], maxs[3];
	point_t geom_min = geom.min_point();
	point_t geom_max = geom.max_point();

	box[0] = geom_min[0];
	box[1] = geom_min[1];
	box[2] = geom_min[2];
	box[3] = geom_max[0];
	box[4] = geom_max[1];
	box[5] = geom_max[2];
      }

    // construct a kd-tree of all the non-isolated mesh_vertices.
    vector<VECTOR3> points;
    vector<Point> pts;
    for(int i = 0; i < mesh.get_nv(); i ++)
      {
	if( mesh.vert_list[i].iso() ) continue;
	Point p = mesh.vert_list[i].point();
	pts.push_back(p);
	points.push_back(VECTOR3(CGAL::to_double(p.x()),
				 CGAL::to_double(p.y()),
				 CGAL::to_double(p.z())));
      }
    KdTree kd_tree(points, 20);
    kd_tree.setNOfNeighbours(1);

    // Now perform a reconstruction to build a tetrahedralized solid
    // with in-out marked.
    Triangulation triang;
    recon(pts, triang);

    // assign weight to each triangle.
    vector<double> weights;
    // assign_sdf_weight(mesh, weights); // comment out for uniform weight.

    volume vol;

    cvcapp.log(3,"SDF");
    try
      {
	vol.voxel_dimensions(dimension(dimx,dimy,dimz));
	vol.voxelType(Float);
	vol.boundingBox(box);

	for(unsigned int k=0; k<vol.ZDim(); k++)
	  {
	    for(unsigned int j=0; j<vol.YDim(); j++)
	      {
		for(unsigned int i=0; i<vol.XDim(); i++)
		  {
		    double x = vol.XMin() + i*vol.XSpan();
		    double y = vol.YMin() + j*vol.YSpan();
		    double z = vol.ZMin() + k*vol.ZSpan();
		    double fn_val = sdf(Point(x,y,z), mesh, weights, kd_tree, triang);
		    vol(i,j,k, fn_val);
		  }
	      }
	    cvcapp.threadProgress(float(k)/float(vol.ZDim()-1));
	  }

	vol.desc("multi_sdf");
      }
    catch(std::exception &e)
      {
	cerr << e.what() << endl;
      }
    cvcapp.log(3,"done.");
    
    return vol;
  }  
#endif

  CVC_NAMESPACE::geometry convert(const FastContouring::TriSurf& geo)
  {
    using namespace std;
    CVC_NAMESPACE::geometry ret_geom;
    ret_geom.points().resize(geo.verts.size()/3);
    memcpy(&(ret_geom.points()[0]),
           &(geo.verts[0]),
           geo.verts.size()*sizeof(double));
    ret_geom.normals().resize(geo.normals.size()/3);
    memcpy(&(ret_geom.normals()[0]),
           &(geo.normals[0]),
           geo.normals.size()*sizeof(double));
    ret_geom.colors().resize(geo.colors.size()/3);
    memcpy(&(ret_geom.colors()[0]),
           &(geo.colors[0]),
           geo.colors.size()*sizeof(double));
    ret_geom.tris().resize(geo.tris.size()/3);
    memcpy(&(ret_geom.tris()[0]),
           &(geo.tris[0]),
           geo.tris.size()*sizeof(unsigned int));
    return ret_geom;
  }
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
        vol = multi_sdf_library(geom,dim,bbox);
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

  // ---
  // iso
  // ---
  // Purpose: 
  //   Returns geometry representing an isosurface of the specified volume.
  // ---- Change History ----
  // 12/29/2013 -- Joe R. -- Creation.
  geometry iso(const volume& vol, double isovalue, double r, double g, double b)
  {
    thread_info ti(BOOST_CURRENT_FUNCTION);

#if 0
    FastContouring::ContourExtractor contourExtractor;
    FastContouring::Volume v;

    double mapped_isoval = 255.0*(isovalue - vol.min())/(vol.max() - vol.min());

    volume localvol = vol;
    localvol.map(0.0,255.0);
    localvol.voxelType(CVC_NAMESPACE::UChar);
    v.data = *localvol;
    v.xdim = vol.XDim(); v.ydim = vol.YDim(); v.zdim = vol.ZDim();
    v.xmin = vol.XMin(); v.ymin = vol.YMin(); v.zmin = vol.ZMin();
    v.xmax = vol.XMax(); v.ymax = vol.YMax(); v.zmax = vol.ZMax();
    contourExtractor.setVolume(v);
    return convert(contourExtractor.extractContour(mapped_isoval,r,g,b));
#endif

    ConDataset* the_data;
    Contour3dData *contour3d;
    int dim[3];
    float span[3], orig[3];

    dim[0] = vol.XDim(); dim[1] = vol.YDim(); dim[2] = vol.ZDim();
    span[0] = vol.XSpan(); span[1] = vol.YSpan(); span[2] = vol.ZSpan();
    orig[0] = vol.XMin(); orig[1] = vol.YMin(); orig[2] = vol.ZMin();
    
    //convert it to a supported libcontour type and load it into libisocontour
    volume localvol = vol;
    switch(localvol.voxelType())
      {
      case UInt:
      case Double:
      case UInt64:
      case Float:
	localvol.voxelType(Float);
	the_data = newDatasetReg(CONTOUR_FLOAT, CONTOUR_REG_3D, 1, 1, dim, *localvol);
	break;
      case UChar:
	the_data = newDatasetReg(CONTOUR_UCHAR, CONTOUR_REG_3D, 1, 1, dim, *localvol);
	break;
      case UShort:
	the_data = newDatasetReg(CONTOUR_USHORT, CONTOUR_REG_3D, 1, 1, dim, *localvol);
	break;
      }

    ((Datareg3 *)the_data->data->getData(0))->setOrig(orig);
    ((Datareg3 *)the_data->data->getData(0))->setSpan(span);

    //actually extract the contour
    contour3d = getContour3d(the_data,0,0,isovalue,NO_COLOR_VARIABLE);

    //convert the result to geometry object
    geometry result;
    points_t& points = result.points(); points.resize(contour3d->nvert);
    colors_t& colors = result.colors(); colors.resize(contour3d->nvert);
    normals_t& normals = result.normals(); normals.resize(contour3d->nvert);
    tris_t& tris = result.tris(); tris.resize(contour3d->ntri);
    for(int i = 0; i < contour3d->nvert; i++)
      {
	for(int j=0; j<3; j++)
	  points[i][j] = contour3d->vert[i][j];
	for(int j=0; j<3; j++)
	  normals[i][j] = contour3d->vnorm[i][j];

	colors[i][0] = r;
	colors[i][1] = g;
	colors[i][2] = b;
      }
    for(int i = 0; i < contour3d->ntri; i++)
      for(int j=0; j<3; j++)
	tris[i][j] = contour3d->tri[i][j];

    delete the_data;
    delete contour3d;

    return result;
  }
}

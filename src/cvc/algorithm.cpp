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

#include <iostream>
#include <cstring>
#include <algorithm>

namespace
{
  CVC_DEF_EXCEPTION(sign_distance_function_error);

  CVC_NAMESPACE::volume sdf_library(const CVC_NAMESPACE::geometry& geom,
				    const CVC_NAMESPACE::dimension& dim,
				    const CVC_NAMESPACE::bounding_box& bbox)
  {
    return CVC_NAMESPACE::volume();
  }

  CVC_NAMESPACE::geometry cvc_mesher(const CVC_NAMESPACE::volume& vol, double isovalue)
  {
    return CVC_NAMESPACE::geometry();
  }
}

namespace CVC_NAMESPACE
{
  // ---
  // sdf
  // ---
  // Purpose: 
  //   Returns a volume representing the signed distance function of the input geometry.
  // ---- Change History ----
  // 12/29/2013 -- Joe R. -- Creation.
  volume sdf(const geometry& geom,
	     const dimension& dim,
	     const bounding_box& bbox)
  {
    thread_info ti(BOOST_CURRENT_FUNCTION);
    volume vol = sdf_library(geom,dim,bbox);
    vol.desc("Signed Distance Function - SDFLibrary");
    return vol;
  }

  // ---
  // iso
  // ---
  // Purpose: 
  //   Returns geometry representing an isosurface of the specified volume.
  // ---- Change History ----
  // 12/29/2013 -- Joe R. -- Creation.
  // 01/08/2014 -- Joe R. -- Removing color args and preparing for cvc-mesher.
  geometry iso(const volume& vol, double isovalue)
  {
    thread_info ti(BOOST_CURRENT_FUNCTION);
    geometry geo = cvc_mesher(vol,isovalue);
    return geo;
  }
}

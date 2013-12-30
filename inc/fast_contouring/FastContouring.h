/*
  Copyright 2002-2003,2008-2011 The University of Texas at Austin

        Authors: Anthony Thane <thanea@ices.utexas.edu>
                 John Wiggins <prok@ices.utexas.edu>
                 Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of FastContouring.

  FastContouring is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  FastContouring is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef __LBIE__FASTCONTOURING_H__
#define __LBIE__FASTCONTOURING_H__

#include <vector>

#include <fast_contouring/Matrix.h>
#include <fast_contouring/MarchingCubesBuffers.h>

#include <boost/cstdint.hpp>

namespace FastContouring
{
  //dead simple contour surface description
  struct TriSurf
  {
    std::vector<double> verts;
    std::vector<double> normals;
    std::vector<double> colors;
    std::vector<unsigned int> tris;
  };

  struct Volume
  {
    boost::uint64_t xdim, ydim, zdim;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    unsigned char* data;
  };

  class ContourExtractor
  {
  public:
    ContourExtractor()
      {
	setVolume(m_Data); //make sure buffers are initialized
      }
    ContourExtractor(const ContourExtractor& copy)
      : m_Data(copy.m_Data),
        m_Buffers(copy.m_Buffers),
        m_SaveMatrix(copy.m_SaveMatrix)
      {}
    ~ContourExtractor() {}
    
    ContourExtractor& operator=(const ContourExtractor& copy)
    {
      m_Data = copy.m_Data;
      m_Buffers = copy.m_Buffers;
      m_SaveMatrix = copy.m_SaveMatrix;
      return *this;
    }

    void setVolume(const Volume& vol);
    const Volume& getVolume() const { return m_Data; }
    
    TriSurf extractContour(double isovalue,
			   double R = 1.0, double G = 1.0, double B = 1.0);

  private:
    void classifyVertices(unsigned int k, unsigned int* cacheMemory, float isovalue) const;
    void getNormal(unsigned int i, unsigned int j, unsigned int k,
		   float& nx, float& ny, float& nz) const;
    unsigned int determineCase(/*unsigned char* data, */
			       unsigned int* offsetTable, unsigned int index
			       /*, GLfloat isovalue*/) const;

    Volume m_Data;
    // buffers used to speed up marching cubes
    MarchingCubesBuffers m_Buffers;
    Matrix m_SaveMatrix;
  };
			 
};

#endif

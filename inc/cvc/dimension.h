/*
  Copyright 2007-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of libCVC.

  libCVC is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  libCVC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: Dimension.h 4741 2011-10-21 21:22:06Z transfix $ */

#ifndef __CVC_DIMENSION_H__
#define __CVC_DIMENSION_H__

#include <cvc/types.h>

// If your compiler complains the "The class "cvc::Dimension" has no member "xdim"."
// Add your architecture Q_OS_XXXX flag (see qglobal.h) in this list.
//#if defined (Q_OS_IRIX) || defined (Q_OS_AIX) || defined (Q_OS_HPUX)
# define UNION_NOT_SUPPORTED
//#endif

namespace CVC_NAMESPACE
{
  class dimension
  {
  public:
    /* The internal data representation is public. */
#if defined (DOXYGEN) || defined (UNION_NOT_SUPPORTED)
    uint64 xdim, ydim, zdim;
#else
    union
    {
      struct { uint64 xdim, ydim, zdim; };
      uint64 dim_[3];
    };
#endif
    
    /* Default constructor */
    dimension() : xdim(0), ydim(0), zdim(0) {}
      
    /* Standard constructor */
    dimension(uint64 x, uint64 y, uint64 z) :
      xdim(x), ydim(y), zdim(z) {}

    /*
      Universal explicit converter from any class to dimension (as long as that class implements
      operator[]).
    */
    template <class C> explicit dimension(const C& m) : 
      xdim(m[0]), ydim(m[1]), zdim(m[2]) {}

    dimension& operator=(const dimension& d)
    {
      xdim = d.xdim; ydim = d.ydim; zdim = d.zdim;
      return *this;
    }

    bool operator==(const dimension& d) const
    {
      return (xdim == d.xdim) && (ydim == d.ydim) && (zdim == d.zdim);
    }

    bool operator!=(const dimension& d) const
    {
      return !((*this)==d);
    }

    bool operator<(const dimension& d) const
    {
      return (*this <= d) && (*this != d);
    }

    bool operator>(const dimension& d) const
    {
      return (*this >= d) && (*this != d);
    }

    bool operator<=(const dimension& d) const
    {
      return (xdim <= d.xdim) && (ydim <= d.ydim) && (zdim <= d.zdim);
    }

    bool operator>=(const dimension& d) const
    {
      return (xdim >= d.xdim) && (ydim >= d.ydim) && (zdim >= d.zdim);
    }

    void setDim(uint64 x, uint64 y, uint64 z) { xdim = x; ydim = y; zdim = z; }

    /* Bracket operator with a constant return value. */
    uint64 operator[](int i) const 
      {
#ifdef UNION_NOT_SUPPORTED
	return (&xdim)[i];
#else
	return dim_[i]; 
#endif      
      }

    /* Bracket operator returning an l-value. */
    uint64& operator[](int i) 
      {
#ifdef UNION_NOT_SUPPORTED
	return (&xdim)[i];
#else 
	return dim_[i];
#endif 
      }

    bool isNull() const { return xdim == 0 && ydim == 0 && zdim == 0; }

    //returns the number of voxels for this dimension
    uint64 size() const { return xdim*ydim*zdim; }

    uint64 XDim() const { return xdim; }
    uint64 YDim() const { return ydim; }
    uint64 ZDim() const { return zdim; }
  };
};

#endif

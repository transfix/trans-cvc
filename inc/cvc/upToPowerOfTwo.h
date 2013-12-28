/*
  Copyright 2005-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of libcvc.

  libcvc is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  libcvc is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef __CVC_UPTOPOWEROFTWO__
#define __CVC_UPTOPOWEROFTWO__

#include <cvc/namespace.h>

namespace CVC_NAMESPACE
{
  static inline unsigned int upToPowerOfTwo(unsigned int value)
    {
      unsigned int c = 0;
      unsigned int v = value;

      // round down to nearest power of two 
      while (v>1) {
        v = v>>1;
        c++;
      }

      // if that isn't exactly the original value 
      if ((v<<c)!=value) {
        // return the next power of two 
        return (v<<(c+1));
      }
      else {
        // return this power of two 
        return (v<<c);
      }
    }
}

#endif

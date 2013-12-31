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
// basic definitions

#ifndef __BASIC_H
#define __BASIC_H

//#include <Utility/utility.h>

typedef unsigned char  u_char;
typedef unsigned short u_short;
typedef unsigned int   u_int;

// Constant used to designate NULL when indexing into an array.
const int INDEX_NULL = -1;


// Return the square of #x#.
template<class T>
T sqr(T x)
{
	return x*x;
}

/// Return the cube of #x#.
template<class T>
T cub(T x)
{
	return x*x*x;
}

#ifdef NO_STL
/// Return the minimum between #x# and #y#.
template<class T>
T min(T x, T y)
{
	return x <= y ? x : y;
}

/// Return the maximum between #x# and #y#.
template<class T>
T max(T x, T y)
{
	return x >= y ? x : y;
#endif

/// Return the sign of #x#.
template<class T>
double sign(T x)
{
	return ((x >= 0.0) ? 1.0 : -1.0);
}

#endif

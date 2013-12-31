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
// endian_io.h - routines to read and (possibly) perform little to big
//		 endian conversions on input data

#ifndef _ENDIAN_IO_H_
#define _ENDIAN_IO_H_

#include <cstdio>

// convert_short() - convert a short int to big endian format
short convert_short(short i);

// convert_long() - convert a long int to big endian format
long convert_long(long i);

// convert_float() - convert a single precision float to big endian format
float convert_float(float i);

// convert_double() - convert double precision real to big endian format
double convert_double(double i);

// fread_short() - read (and possibly convert) short integer data
size_t fread_short(void* ptr, size_t size, size_t nitems, std::FILE* stream);

// fread_int() - read (and possibly convert) long integer data
size_t fread_int(void* ptr, size_t size, size_t nitems, std::FILE* stream);

// fread_float() - read (and possibly convert) single precision data
size_t fread_float(void* ptr, size_t size, size_t nitems, std::FILE* stream);

// fread_double() - read (and possibly convert) double precision data
size_t fread_double(void* ptr, size_t size, size_t nitems, std::FILE* stream);

#endif /* of _ENDIAN_IO_H_ */

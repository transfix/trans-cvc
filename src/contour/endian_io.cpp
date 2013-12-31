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
#include <contour/endian_io.h>

// convert_short() - convert a short int to big endian format
short convert_short(short i)
{
	return (((i & 0xff) << 8) | ((i & 0xff00) >> 8));
}

// convert_long() - convert a long int to big endian format
long convert_long(long i)
{
	return (((i & 0xff) << 24) | ((i & 0xff00) << 8) |
			((i & 0xff0000) >> 8) | ((i & 0xff000000) >> 24));
}

// convert_float() - convert a single precision float to big endian format
float convert_float(float i)
{
	float f;
	unsigned char* iptr = (unsigned char*)&i;
	unsigned char* optr = (unsigned char*)&f;
	*(optr+3) = *iptr;
	*(optr+2) = *(iptr+1);
	*(optr+1) = *(iptr+2);
	*optr = *(iptr+3);
	return f;
}

// convert_double() - convert double precision real to big endian format
double convert_double(double i)
{
	double d;
	unsigned char* iptr = (unsigned char*)&i;
	unsigned char* optr = (unsigned char*)&d;
	*(optr+7) = *iptr;
	*(optr+6) = *(iptr+1);
	*(optr+5) = *(iptr+2);
	*(optr+4) = *(iptr+3);
	*(optr+3) = *(iptr+4);
	*(optr+2) = *(iptr+5);
	*(optr+1) = *(iptr+6);
	*optr = *(iptr+7);
	return d;
}

// fread_short() - read (and possibly convert) short integer data
size_t fread_short(void* ptr, size_t size, size_t nitems, FILE* stream)
{
#ifdef _LITTLE_ENDIAN
	unsigned int i;				// an item index variable
	short*	item;				// an item pointer index
	size_t	ret_size;			// number of items read
	ret_size = fread(ptr, size, nitems, stream);
	item = (short*)ptr;
	for(i = 0; i < ret_size; i++)
	{
		*item = convert_short(*item);
		item++;
	}
	return ret_size;
#else
	return fread(ptr, size, nitems, stream);
#endif
}

// fread_int() - read (and possibly convert) long integer data
size_t fread_int(void* ptr, size_t size, size_t nitems, FILE* stream)
{
#ifdef _LITTLE_ENDIAN
	unsigned int i;				// an item index variable
	long*	item;				// an item pointer index
	size_t	ret_size;			// number of items read
	ret_size = fread(ptr, size, nitems, stream);
	item = (long*)ptr;
	for(i = 0; i < ret_size; i++)
	{
		*item = convert_long(*item);
		item++;
	}
	return ret_size;
#else
	return fread(ptr, size, nitems, stream);
#endif
}

// fread_float() - read (and possibly convert) single precision data
size_t fread_float(void* ptr, size_t size, size_t nitems, FILE* stream)
{
#ifdef _LITTLE_ENDIAN
	unsigned int i;				// an item index variable
	float*	item;				// an item pointer index
	size_t	ret_size;			// number of items read
	ret_size = fread(ptr, size, nitems, stream);
	item = (float*)ptr;
	for(i = 0; i < ret_size; i++)
	{
		*item = convert_float(*item);
		item++;
	}
	return ret_size;
#else
	return fread(ptr, size, nitems, stream);
#endif
}

// fread_double() - read (and possibly convert) double precision data
size_t fread_double(void* ptr, size_t size, size_t nitems, FILE* stream)
{
#ifdef _LITTLE_ENDIAN
	unsigned int i;				// an item index variable
	double*	item;				// an item pointer index
	size_t	ret_size;			// number of items read
	ret_size = fread(ptr, size, nitems, stream);
	item = (double*)ptr;
	for(i = 0; i < ret_size; i++)
	{
		*item = convert_double(*item);
		item++;
	}
	return ret_size;
#else
	return fread(ptr, size, nitems, stream);
#endif
}

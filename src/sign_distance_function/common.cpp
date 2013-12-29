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

#include <math.h>

#include <sign_distance_function/common.h>

#include <time.h>
#ifdef _WIN32
         #include <sys/types.h>
         #include <sys/timeb.h>
#else
         #include <sys/time.h>
#endif
 

int SDFLibrary::isEqual (double one, double two)
{
	if ( (-1*SDFLibrary::TOLERANCE <= (one-two)) && ((one-two) <= SDFLibrary::TOLERANCE) )
		return true;
	return false;
}

int SDFLibrary::isZero(double num)
{
	if ( (-1*SDFLibrary::TOLERANCE <= num) && (num <= SDFLibrary::TOLERANCE) )
		return true;
	return false;
}

int SDFLibrary::isNegative(double num)
{
	if (num <0)	return true;
	return false;
}

int SDFLibrary::isBetween(double one, double two, double num)
{
	if ( ((one<=num) && (num<=two)) || ((isEqual(num, one)) || (isEqual(num, two))) )
		return true;
	return false;
}

int SDFLibrary::isZero(SDFLibrary::myPoint one)
{
	double val = sqrt(one.x*one.x + one.y*one.y + one.z*one.z);

	if (isZero(val))
		return true;
	return false;
}

int SDFLibrary::isSame(SDFLibrary::myPoint one, SDFLibrary::myPoint two)
{
	double val = sqrt( (one.x-two.x)*(one.x-two.x) + (one.y-two.y)*(one.y-two.y) + (one.z-two.z)*(one.z-two.z) );

	if (isZero(val))
		return true;
	return false;
}



void SDFLibrary::_vert2index(int c, int &i, int &j, int &k)
{
	int _left;

	i = c%(SDFLibrary::size+1);

	_left = c/(SDFLibrary::size+1);
	j = _left%(SDFLibrary::size+1);

	_left = _left/(SDFLibrary::size+1);
	k = _left;

	if (i<0) i=0;		
	if (j<0) j=0;			
	if (k<0) k=0;
	if (i>SDFLibrary::size+1) i=SDFLibrary::size+1;	
	if (j>SDFLibrary::size+1) j=SDFLibrary::size+1;		
	if (k>SDFLibrary::size+1) k=SDFLibrary::size+1;
}

int SDFLibrary::index2vert(int i, int j, int k)
{
	return(k*(SDFLibrary::size+1)*(SDFLibrary::size+1) + j*(SDFLibrary::size+1) + i);
}

void SDFLibrary::_cell2index(int c, int &i, int &j, int &k)
{
	int _left;

	i = c%(SDFLibrary::size);

	_left = c/(SDFLibrary::size);
	j = _left%(SDFLibrary::size);

	_left = _left/(SDFLibrary::size);
	k = _left;

	if (i<0) i=0;		
	if (j<0) j=0;			
	if (k<0) k=0;
	if (i>SDFLibrary::size) i=SDFLibrary::size;	
	if (j>SDFLibrary::size) j=SDFLibrary::size;		
	if (k>SDFLibrary::size) k=SDFLibrary::size;
}

int SDFLibrary::index2cell(int i, int j, int k)
{
	return(k*(SDFLibrary::size)*(SDFLibrary::size) + j*(SDFLibrary::size) + i);
}

double SDFLibrary::xCoord(int i)
{
	if ((0<=i) && (i<= SDFLibrary::size+1))
		return ((double)(SDFLibrary::minext[0] + i*SDFLibrary::span[0]));

	printf("X Index out of bounds. now exiting %d\n", i);
	return -1;
}

double SDFLibrary::yCoord(int j)
{
	if ((0<=j) && (j<= SDFLibrary::size+1))
		return ((double)(SDFLibrary::minext[1] + j*SDFLibrary::span[1]));

	printf("Y Index out of bounds. now exiting %d\n", j);
	return -1;
} 

double SDFLibrary::zCoord(int k)
{
	if ((0<=k) && (k<= SDFLibrary::size+1))
		return ((double)(SDFLibrary::minext[2] + k*SDFLibrary::span[2]));

	printf("Z Index out of bounds. now exiting %d\n", k);
	return -1;
}

void SDFLibrary::object2octree(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int &ci, int &cj, int &ck)
{
	//Get the indices of the Closest vertex of the Octree cell by interpolating 'tween in the Object space
	int i, j, k;

	ci = (int)((xmin-SDFLibrary::minext[0])/(SDFLibrary::span[0]));
	i = (int)((xmax-SDFLibrary::minext[0])/(SDFLibrary::span[0]));
	ci = (i+ci)/2;

	cj = (int)((ymin-SDFLibrary::minext[1])/(SDFLibrary::span[1]));
	j = (int)((ymax-SDFLibrary::minext[1])/(SDFLibrary::span[1]));
	cj = (j+cj)/2;

	ck = (int)((zmin-SDFLibrary::minext[2])/(SDFLibrary::span[2]));
	k = (int)((zmax-SDFLibrary::minext[2])/(SDFLibrary::span[2]));
	ck = (k+ck)/2;
	
	if ( (i!=ci+1) || (j!=cj+1) || (k!=ck+1) )
		printf("cannot make a good Octree\n");
}

//Get the current time in seconds as a double value 
double SDFLibrary::getTime()
{
	#ifdef _WIN32
		 time_t ltime;
		 _timeb tstruct;
		 time( &ltime );
		 _ftime( &tstruct );
		 return (double) (ltime + 1e-3*(tstruct.millitm));
	#else
		 struct timeval t;
		 gettimeofday( &t, NULL );
		 return (double)(t.tv_sec + 1e-6*t.tv_usec);
	#endif
}

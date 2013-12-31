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
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <contour/range.h>

#ifndef WIN32
#include <unistd.h>
#endif

#define MIN2(x,y) ((x)<(y)?(x):(y))
#define MAX2(x,y) ((x)>(y)?(x):(y))

//  Check - check for consistency in range
void Range::Check(void)
{
	for(int i=0; i<nrange; i++)
	{
		if(min[i] > max[i])
		{
			printf("invalid range!\n");
#ifndef WIN32
			sleep(3);
#endif
		}
		if(i<nrange-1 && max[i] > min[i+1])
		{
			printf("invalid range(s)!\n");
#ifndef WIN32
			sleep(3);
#endif
		}
	}
}

//  Compress - compress a range
void Range::Compress(void)
{
	int j;
	for(int i=0; i<nrange; i++)
	{
		for(j=i+1; j<nrange && min[j] <= max[i]; j++)
		{
			/* scan over overlapping regions */
			if(max[j] > max[i])
			{
				max[i] = max[j];
			}
		}
		if(j != i+1)
		{
			memcpy(&min[i+1], &min[j], sizeof(float)*(j-(i+1)));
			memcpy(&max[i+1], &max[j], sizeof(float)*(j-(i+1)));
			nrange -= (j-(i+1));
		}
	}
}

//  Print - print a range
void Range::Print()
{
	if(Empty())
	{
		printf("empty\n");
		return;
	}
	for(int i=0; i<nrange; i++)
	{
		printf("%f->%f%s", min[i], max[i], (i==(nrange-1))?"\n":", ");
	}
}

//  AddRange - add a range to a Range
void Range::AddRange(float min_val, float max_val)
{
	int i;
	for(i=nrange-1; i>=0 && min[i] > min_val; i--)
	{
		min[i+1] = min[i];
		max[i+1] = max[i];
	}
	min[i+1] = min_val;
	max[i+1] = max_val;
	nrange++;
}

//  += - add the given range
Range& Range::operator+=(const Range& r)
{
	int i;
	for(i=0; i<r.nrange; i++)
	{
		AddRange(r.min[i], r.max[i]);
	}
	Compress();
#ifdef CHECK
	Check();
#endif
	return(*this);
}

// Difference --  subtract argument from 'this'
Range& Range::operator-=(const Range& r2)
{
	static Range r;
	Range* result = &r;
	int i, j, n;
	float curmin;
	i = 0;
	j = 0;
	result->nrange = 0;
	if(nrange != 0)
	{
		curmin = min[0];
		while(i < nrange && j < r2.nrange)
		{
			if(curmin <= r2.min[j])
			{
				if(max[i] < r2.min[j])
				{
					/* the remainder of this can be taken */
					n = result->nrange++;
					result->min[n] = curmin;
					result->max[n] = max[i];
					i++;
					if(i < nrange)
					{
						curmin = min[i];
					}
				}
				else
				{
					if(max[i] < r2.max[j])
					{
						if(curmin != r2.min[j])
						{
							n = result->nrange++;
							result->min[n] = curmin;
							result->max[n] = r2.min[j];
						}
						i++;
						if(i < nrange)
						{
							curmin = min[i];
						}
					}
					else
					{
						if(curmin != r2.min[j])
						{
							n = result->nrange++;
							result->min[n] = curmin;
							result->max[n] = r2.min[j];
						}
						curmin = r2.max[j];
						j++;
						if(curmin == max[i])
						{
							i++;
							if(i < nrange)
							{
								curmin = min[i];
							}
						}
					}
				}
			}
			else
			{
				if(max[i] <= r2.max[j])
				{
					i++;
					if(i < nrange)
					{
						curmin = min[i];
					}
				}
				else if(r2.max[j] <= curmin)
				{
					j++;
				}
				else
				{
					curmin = r2.max[j];
					j++;
					if(curmin == max[i])
					{
						i++;
						if(i < nrange)
						{
							curmin = min[i];
						}
					}
				}
			}
		}
		/* anything left in this remains in result */
		while(i < nrange)
		{
			n = result->nrange++;
			result->min[n] = curmin;
			result->max[n] = max[i];
			if(i+1 < nrange)
			{
				curmin = min[i+1];
			}
			i++;
		}
#ifdef CHECK
		result->Check();
#endif
		*this = *result;
	}
	return(*this);
}

// Intersect --  ***** WARNING - this is broken
Range& Range::operator^(const Range& r2)
{
	static Range r;
	Range* result = &r;
	int i, j, n;
	i = 0;
	j = 0;
	result->nrange = 0;
	while(i < nrange && j < r2.nrange)
	{
		if(min[i] <= r2.min[j])
		{
			if(max[i] >= r2.min[j])
			{
				if(max[i] < r2.max[j])
				{
					n = result->nrange++;
					result->min[n] = r2.min[j];
					result->max[n] = max[i];
					i++;
				}
				else
				{
					n = result->nrange++;
					result->min[n] = r2.min[j];
					result->max[n] = r2.max[j];
					j++;
				}
			}
			else
			{
				i++;
			}
		}
		else
		{
			if(r2.max[j] >= min[i])
			{
				if(r2.max[j] < max[i])
				{
					n = result->nrange++;
					result->min[n] = min[i];
					result->max[n] = r2.max[j];
					j++;
				}
				else
				{
					n = result->nrange++;
					result->min[n] = min[i];
					result->max[n] = max[i];
				}
			}
			else
			{
				j++;
			}
		}
	}
#ifdef CHECK
	RangeCheck(result);
#endif
	return(*result);
}

// Equal - test for equivalence of two ranges
int Range::operator==(const Range& r2) const
{
	if(nrange != r2.nrange)
	{
		return(0);
	}
	for(int i=0; i<nrange; i++)
	{
		if(min[i] != r2.min[i] || max[i] != r2.max[i])
		{
			return(0);
		}
	}
	return(1);
}

// Complement - take the complement with respect to a range
Range Range::Complement(float min_val, float max_val)
{
	static Range r;
	Range* result=&r;
	int i, n;
	/* handle complement of empty set */
	if(Empty())
	{
		result->nrange = 1;
		result->min[0] = min_val;
		result->max[0] = max_val;
		return(*result);
	}
	result->nrange = 0;
	if(min_val < min[0])
	{
		n = result->nrange++;
		result->min[n] = min_val;
		result->max[n] = min[0];
	}
	for(i=0; i<nrange-1; i++)
	{
		n = result->nrange++;
		result->min[n] = max[i];
		result->max[n] = min[i+1];
	}
	if(max_val > max[nrange-1])
	{
		n = result->nrange++;
		result->min[n] = max[nrange-1];
		result->max[n] = max_val;
	}
#ifdef CHECK
	result->Check();
#endif
	return(r);
}

// Disjoint - test if two ranges are disjoint
int Range::Disjoint(const Range& r2) const
{
	static Range diff;
	diff = (*this) - r2;
	return((diff == *this));

}

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
// edgeHash.C - hash to lookup vertices which are already computed
// This is a *very* basic hash class to aid in finding vertices which
// have already been computed, avoiding recomputation and duplication
// Copyright (c) 1997 Dan Schikore

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <contour/edgehash.h>

#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif

// use a small prime # of buckets
#define NBUCKETS 3001

extern int verbose;

// EdgeHash - construct a new hash
EdgeHash::EdgeHash()
{
	nbuckets = NBUCKETS;
	nitems  = (int*)malloc(sizeof(int) * nbuckets);
	buckets  = (EdgeHashBucket*)malloc(sizeof(EdgeHashBucket) * nbuckets);
	// initialize each bucket
	for(int b=0; b<nbuckets; b++)
	{
		nitems[b] = 0;
		buckets[b].elsize = 5;
		buckets[b].items = (EdgeHashEl*)malloc(sizeof(EdgeHashEl)* buckets[b].elsize);
	}
}

// LookupBucket - search a given bucket for a given key
int EdgeHash::LookupBucket(int* nitems, EdgeHashBucket* b, int key)
{
	int vnum;
	// loop through the items
	for(int i=0; i<(*nitems); i++)
	{
		if(b->items[i].key == key)
		{
			// found the requested key
			vnum = b->items[i].vnum;
			if(++(b->items[i].nref) == 4)
			{
				// edges referenced 4 times will not be used again
				if((*nitems) > 1)
				{
					b->items[i] = b->items[(*nitems)-1];
				}
				(*nitems)--;
			}
			return(vnum);
		}
	}
	return(-1);
}

// InsertBucket - insert an item in the given bucket
void EdgeHash::InsertBucket(int* nitems, EdgeHashBucket* b, int key, int vnum)
{
	int n = (*nitems)++;
	if(n >= b->elsize)
	{
		b->elsize*=2;
		b->items = (EdgeHashEl*)realloc(b->items, sizeof(EdgeHashEl)* b->elsize);
		if(verbose > 1)
		{
			printf("hash size: %d\n", b->elsize);
		}
	}
	b->items[n].key  = key;
	b->items[n].vnum = vnum;
	b->items[n].nref = 1;
}

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
// edgeHash.h - hash to lookup vertices which are already computed
//
// This is a *very* basic hash class to aid in finding vertices which
// have already been computed, avoiding recomputation and duplication

#ifndef EDGE_HASH_H
#define EDGE_HASH_H

#include <cstdlib>

// a bucket entry containing a key, vertex number, and the number of times referenced so far.
// We assume that a vertex referenced four times can be discarded.
typedef struct EdgeHashElT
{
	int key;
	int vnum;
	int nref;
} EdgeHashEl;


// a hash bucket.  Each bucket contains a dynamic list of vertices which  map to the given bucket
typedef struct EdgeHashBucketT
{
	int elsize;
	EdgeHashEl* items;
} EdgeHashBucket;


// EdgeHash
typedef struct EdgeHash
{
	public:
		EdgeHash();
		~EdgeHash()
		{
			std::free(nitems);
			std::free(buckets);
		}

		// insert an item (key, vertex number) into the hash
		void Insert(int key, int vnum)
		{
			int bucket=CompBucket(key);
			InsertBucket(&nitems[bucket], &buckets[bucket], key, vnum);
		}

		// lookup the vertex with the given key
		int  Lookup(int key)
		{
			int bucket = CompBucket(key);
			return(LookupBucket(&nitems[bucket], &buckets[bucket], key));
		}

		// reset the hash
		void Reset(void)
		{
			memset(nitems, 0, sizeof(int)*nbuckets);
		}


	protected:
		// compute the bucket mapping for a given key
		int  CompBucket(int key)
		{
			return (key % nbuckets);
		}

		// search a given bucket for a given key
		int  LookupBucket(int*, EdgeHashBucket*, int);

		// insert a key in the given bucket
		void InsertBucket(int*, EdgeHashBucket*, int, int);

	private:
		// number of buckets
		int nbuckets;

		// number of items for each bucket
		int* nitems;

		// buckets of storage
		EdgeHashBucket* buckets;
} EdgeHash;

#endif


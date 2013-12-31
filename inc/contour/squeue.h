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
#ifndef __SQUEUE_H
#define __SQUEUE_H

#include <contour/queue.h>

//@ManMemo: A templated FIFO queue
/*@Doc: Implementation of a FIFO queue based on a circular array.
  Memory is allocated in blocks, given as a parameter when the
  object is created, or if no parameter is given, using the default
  block size.  Note: Currently, when a new block of memory is
  allocated, all of the current SQueue is recopied (because it is
  circular you cannot just increase the array size), so use caution
  when picking the block size. */
template <class T>
class SQueue : public Queue<T>
{
	public:
		//@ManDoc: Constructor with user define block size.
		SQueue(int blocksize=0);

		//@ManDoc: Destructor.
		virtual ~SQueue();

		int find(T&);

	private:
};

template <class T>
SQueue< T >::SQueue(int blocksize) : Queue<T>(blocksize)
{
}

template <class T>
SQueue< T >::~SQueue()
{
}

template <class T>
int SQueue< T >::find(T& e)
{
	int i, j;
	for(j = this->head, i = 0; i < this->length; i++)
	{
		if(this->q[j] == e)
		{
			return(j);
		}
		j++;
		if(j == this->room)
		{
			j = 0;
		}
	}
	return(-1);
}

#endif

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
#ifndef __QUEUE_H
#define __QUEUE_H

//#include <Utility/utility.h>
#include <contour/utilities.h>

//@ManMemo: A templated FIFO queue
/*@Doc: Implementation of a FIFO queue based on a circular array.
  Memory is allocated in blocks, given as a parameter when the
  object is created, or if no parameter is given, using the default
  block size.  Note: Currently, when a new block of memory is
  allocated, all of the current Queue is recopied (because it is
  circular you cannot just increase the array size), so use caution
  when picking the block size. */
template <class T>
class Queue
{
	public:
		//@ManDoc: Constructor with user define block size.
		Queue(int blocksize=0);

		//@ManDoc: Destructor.
		virtual ~Queue();

		//@ManDoc: Returns whether the Queue is empty or not.
		inline bool isEmpty();

		//@ManDoc: Returns the element at the head, or NULL if the Queue is empty.
		inline T* first();

		//@ManDoc: Places an item at the end of the queue.
		int enqueue(T& e);

		//@ManDoc: Remove \& return the item at the head of the array.
		T* dequeue();

		T* getItem(int i)
		{
			return(&q[i]);
		}

		int getHead(void)
		{
			return(head);
		}
		int getLength(void)
		{
			return(length);
		}
		int getRoom(void)
		{
			return(room);
		}

	protected:
		int head;			// Beginning of the circular array
		int tail;			// End of the circular array
		int length;			// Length of the array
		int block;			// Number of Ts to allocate memory for
		int room;			// Total number of Ts memory is allocated for
		T* q;				// Actual array that contains memory
};

template <class T>
bool Queue< T >::isEmpty()
{
	return (bool)(length == 0);
}

template <class T>
T* Queue< T >::first()
{
	if(length == 0)
	{
		return 0;
	}
	return &q[head];
}

template <class T>
Queue< T >::Queue(int blocksize)
{
	head   = 0;
	tail   = 0;
	length = 0;
	if(blocksize == 0)
	{
		block = (size_t(1) > size_t(4096/sizeof(T)) ?
				 size_t(1) : size_t(4096/sizeof(T)));
	}
	else
	{
		block = blocksize;
	}
	room   = block;
	q = (T*)malloc(room*sizeof(T));
}

template <class T>
Queue< T >::~Queue()
{
	free(q);
}

template <class T>
int Queue< T >::enqueue(T& e)
{
	int i, j;
	int loc;
	T* t;
	// enough room ?
	if(length == room)
	{
		t = (T*)malloc((room+block)*sizeof(T));
		for(j = head, i = 0; i < length; i++)
		{
			t[i] = q[j];
			j++;
			if(j == room)
			{
				j = 0;
			}
		}
		room += block;
		free(q);
		q = t;
		head = 0;
		tail = length;
	}
	// insert new element
	q[loc=tail] = e;
	length++;
	tail++;
	if(tail == room)
	{
		tail = 0;
	}
	return(loc);
}

template <class T>
T* Queue< T >::dequeue()
{
	T* t;
	if(length == 0)
	{
		return (T*)0;
	}
	else
	{
		t = &q[head];
		length--;
		head++;
		if(head == room)
		{
			head = 0;
		}
		return t;
	}
}

#endif

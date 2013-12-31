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
#ifndef __SHELF_H
#define __SHELF_H

#include <contour/basic.h>
#include <iostream>
#include <stddef.h>

#ifdef sparc
#define inline
#else
#define INLINE
#endif

template<class T> class Shelf;

template<class T> class ShelfItem
{
		friend class Shelf<T>;

		//friend std::ostream& operator<<(std::ostream& os, const Shelf<T>& s);
	private:
		T data;
		int prev;
		int next;
};

/** A templated double-linked list with direct access.
 Shelf is a simple templated container class. Objects of the
  templated class T can be added or removed from the Shelf. When an
  object is created, an integer index is returned, and this index is
  guaranteed not to change until deletion.  Memory for objects in the
  shelf is allocated upon need. Memory previously used by removed
  objects is reused. Notice however that a shelf never returns
  allocated memory to the free store. Objects are kept in a
  double-linked list, so sequencing through all objects in the shelf
  is efficiently supported.  Two Shelf objects to which exactly the
  same sequence of add/remove operations is applied, are guaranteed to
  return the same index i at the next add operation. */
template<class T>
class Shelf
{
	public:
		// Constructor.
		Shelf(int blocksize=0);

		// Copy constructor.
		Shelf(const Shelf<T> &shelf);

		// Destructor.
		~Shelf();

		// Assignment (returns a copy of the shelf and its contents).
		Shelf<T>& operator=(const Shelf<T> &s);

		// Add an item to the shelf (return index).
		int put(const T& data);

		// Add an uninitialized item to the shelf (return index).
		int put();

		// Remove item #index#.
		void remove(int index);

		// Reinitialize data structure.
		void cleanUp();

		// Move content of the shelf to a new shelf.
		void moveTo(Shelf<T>& newshelf);

		// Return the number of items currently in the shelf.
		int numItems() const;

		// Return item #index#.
		T& operator[](int index);

		// Return item #index#.
		const T& operator[](int index) const;

		// Iterator. Init #index# to first item.
		void first(int& index) const;

		// Iterator. Init #index# to last item.
		void last(int& index) const;

		// Iterator. Test last-item condition.
		int end(int index) const;

		// Iterator. Update #index# to next item.
		void next(int& index) const;

		// Iterator. Update #index# to previous item.
		void prev(int& index) const;

		// for debugging
		//friend std::ostream& operator<<(std::ostream& os, const Shelf<T>& s);

	private:
		void grow();				// allocate more memory
		void init(int);
		void copy(const Shelf<T> &s);
		void destroy();
		ShelfItem<T>& item(int index);
		const ShelfItem<T>& item(int index) const;

		ShelfItem<T>** _shelfbuffer;
		int _blockSize;
		int _first;
		int _last;
		int _freeList;
		int _nItems;
		int _nBuffers;
		int _lastBuffer;
};

template<class T>
Shelf<T>::Shelf(const Shelf<T> &shelf)
{
	copy(shelf);
}

template<class T>
Shelf<T>& Shelf<T>::operator=(const Shelf<T> &s)
{
	if(&s == this)
	{
		return *this;
	}
	destroy();
	copy(s);
	return *this;
}

template<class T>
ShelfItem<T>& Shelf<T>::item(int index)
{
	return _shelfbuffer[index / _blockSize][index % _blockSize];
}

template<class T>
const ShelfItem<T>& Shelf<T>::item(int index) const
{
	return _shelfbuffer[index / _blockSize][index % _blockSize];
}

template<class T>
T& Shelf<T>::operator[](int index)
{
	return _shelfbuffer[index / _blockSize][index % _blockSize].data;
}

template<class T>
const T& Shelf<T>::operator[](int index) const
{
	return _shelfbuffer[index / _blockSize][index % _blockSize].data;
}

template<class T>
int Shelf<T>::numItems() const
{
	return _nItems;
}

template<class T>
void Shelf<T>::first(int& index) const
{
	index = _first;
}

template<class T>
void Shelf<T>::last(int& index) const
{
	index = _last;
}

template<class T>
int Shelf<T>::end(int index) const
{
	return (index == INDEX_NULL);
}

template<class T>
void Shelf<T>::next(int& index) const
{
	index = item(index).next;
}

template<class T>
void Shelf<T>::prev(int& index) const
{
	index = item(index).prev;
}

template<class T>
Shelf<T>::Shelf(int bs)
{
	init(bs);
}

template<class T>
Shelf<T>::~Shelf()
{
	destroy();
}

template<class T>
void Shelf<T>::cleanUp()
{
	int bs = _blockSize;
	destroy();
	init(bs);
}

template<class T>
void Shelf<T>::moveTo(Shelf<T>& newshelf)
{
	memcpy(&newshelf, this, sizeof(Shelf<T>));     // do not use operator= !
	init(_blockSize);
}

template<class T>
void Shelf<T>::init(int bs)
{
	if(bs == 0)
	{
		_blockSize = (size_t(1) > size_t(4096/sizeof(T)) ?
					  size_t(1) : size_t(4096/sizeof(T)));
	}
	else
	{
		_blockSize = bs;
	}
	_nItems = 0;
	_shelfbuffer = 0;
	_first = INDEX_NULL;
	_last = INDEX_NULL;
	_freeList = INDEX_NULL;
	_nBuffers = 0;
	_lastBuffer = -1;
}

template<class T>
void Shelf<T>::copy(const Shelf<T> &shelf)
{
	int i;
	init(shelf._blockSize);
	for(shelf.first(i); !shelf.end(i); shelf.next(i))
	{
		put(shelf[i]);
	}
}

template<class T>
void Shelf<T>::destroy()
{
	int i;
	for(first(i); !end(i); next(i))
	{
		item(i).data.~T();
	}
	for(i = 0; i <= _lastBuffer; i++)
	{
		::operator delete(_shelfbuffer[i]);
	}
	delete[] _shelfbuffer;
}

template<class T>
int Shelf<T>::put(const T& data)
{
	// get one chunk of memory
	if(_freeList == INDEX_NULL)
	{
		grow();
	}
	int index = _freeList;
	_freeList = item(_freeList).next;
	// copy data
	item(index).data = data;
	// update pointers
	if(_nItems == 0)
	{
		_first = index;
	}
	else
	{
		item(_last).next = index;
	}
	item(index).prev = _last;
	item(index).next = INDEX_NULL;
	_last = index;
	++_nItems;
	return _last;
}

template<class T>
int Shelf<T>::put()
// leave data uninitialized
{
	// get one chunk of memory
	if(_freeList == INDEX_NULL)
	{
		grow();
	}
	int index = _freeList;
	_freeList = item(_freeList).next;
	// update pointers
	if(_nItems == 0)
	{
		_first = index;
	}
	else
	{
		item(_last).next = index;
	}
	item(index).prev = _last;
	item(index).next = INDEX_NULL;
	_last = index;
	++_nItems;
	return _last;
}

template<class T>
void Shelf<T>::remove(int index)
{
	if(item(index).prev == INDEX_NULL)
	{
		_first = item(index).next;
	}
	else
	{
		item(item(index).prev).next = item(index).next;
	}
	if(item(index).next == INDEX_NULL)
	{
		_last = item(index).prev;
	}
	else
	{
		item(item(index).next).prev = item(index).prev;
	}
	item(index).data.~T();
	item(index).next = _freeList;
	_freeList = index;
	--_nItems;
}

template<class T>
void Shelf<T>::grow()
{
	// do we have room for a new buffer pointer?
	if(++_lastBuffer == _nBuffers)
	{
		if(_shelfbuffer != 0)
		{
			ShelfItem<T>** tmp = _shelfbuffer;
			_shelfbuffer = new ShelfItem<T>*[_nBuffers+10];
			for(int i = 0; i < _nBuffers; i++)
			{
				_shelfbuffer[i] = tmp[i];
			}
			delete[] tmp;
			_nBuffers+=10;
		}
		else
		{
			_shelfbuffer = new ShelfItem<T>*[_nBuffers+=10];
		}
	}
	// allocate buffer and initialize free list
	_shelfbuffer[_lastBuffer] = (ShelfItem<T>*)
								::operator new(size_t(_blockSize*sizeof(ShelfItem<T>)));
	for(int i = 0; i < _blockSize-1; i++)
	{
		_shelfbuffer[_lastBuffer][i].next = _lastBuffer*_blockSize+i+1;
	}
	_shelfbuffer[_lastBuffer][_blockSize-1].next = _freeList;
	_freeList = _lastBuffer*_blockSize;
}

//Craig: Commented out this poorly written debugging function because it was causing template errors
/*
template<class T>
std::ostream& operator<<(std::ostream& os, const Shelf<T>& s)
{
	int i;
	os << "Data" << std::endl;
	for(s.first(i); !s.end(i); s.next(i))
	{
		os << s[i] << ' ';
	}
	os << std::endl << "Free" << std::endl;
	for(i = s._freeList; i != INDEX_NULL; i = s.item(i).next)
	{
		os << i << ' ';
	}
	return os << std::endl;
}
*/

#endif

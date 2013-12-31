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
// A Bin is basically a dynamic array.
// Elements can be added to the Bin using the add method.
// Memory is allocated automatically upon need.

#ifndef __BIN_H
#define __BIN_H

//#include <Utility/utility.h>

#include <cstdlib>
#include <assert.h>

template<class T> class Bin
{
	public:
		// Contructor. Creates a Bin with block size equal to blocksize.
		Bin(int blocksize=0);

		// Copy constructor. Copies the bin and its contents.
		Bin(const Bin<T> &bin)
		{
			*this = bin;
		}

		// Assignment. Returns a copy of the bin and its contents.
		Bin<T> &operator=(const Bin<T> &bin);

		// Destructor.
		~Bin();

		// Change the blocksize to be used from now on.
		void setBlocksize(int blocksize=0);

		// Set the offset.
		void setOffset(int off);

		// Add T to the end of the Bin.
		int add(const T& e);

		// Add an element to the Bin. The default constructor is used to create the element.
		int add()
		{
			return add(T());
		}

		// Remove last #num# items.
		int removeLast(int num=1);

		int removePos(int pos=0);

		// Returns the number of items in the Bin.
		int numItems() const;

		// Returns the T stored at position k in the bin.
		T& operator[](unsigned int k) const;

		// Resizes memory and returns a pointer to the array. The Bin is reinitialized
		T* done();

		// Returns a pointer without resizing memory.
		T* array() const;

		// Frees memory and reinitializes the Bin.
		void cleanUp();

		// Transfer the contents of this to newbin and clean this.
		void moveTo(Bin<T> &newbin);


	private:
		int offset;			// Position of zero-th element
		int nitems;			// Number of items in the Bin
		int block;			// Size of chunks for each realloc
		int room;			// Size of alloc'd memory
		T* a;				// Array of Ts

		void init(int blocksize);	// init the bin
		void destroy();			// destroy the bin and its contents
		void grow();		        // create additional space
};

template<class T>
void Bin<T>::init(int blocksize)
{
	if(blocksize == 0)
	{
		block = (size_t(1) > size_t(4096/sizeof(T)) ?
				 size_t(1) : size_t(4096/sizeof(T)));
	}
	else
	{
		block = blocksize;
	}
	offset = 0;
	nitems = 0;
	room = 0;
	a = 0;
}

template<class T>
void Bin<T>::destroy()
{
	for(int i = 0; i < nitems; i++)
	{
		a[i].~T();
	}
	free(a);
}

// Bin constructor
// blocksize: the number of Ts to allocate memory for
template<class T>
Bin<T>::Bin(int blocksize)
{
	init(blocksize);
}

template<class T>
Bin<T> &Bin<T>::operator=(const Bin<T> &bin)
{
	block = bin.block;
	offset = bin.offset;
	nitems = bin.nitems;
	room = bin.room;
	if(nitems > 0)
	{
		a = (T*)std::malloc(room*sizeof(T));
		for(int i = 0; i < nitems; i++)
		{
			a[i] = bin.a[i];
		}
	}
	else
	{
		a = NULL;
	}
	return *this;
}

template<class T>
Bin<T>::~Bin()
{
	destroy();
}

template<class T>
void Bin<T>::cleanUp()
{
	int bs = block;
	destroy();
	init(bs);
}

// returns the number of items in the bin
template<class T>
int Bin<T>::numItems() const
{
	return nitems;
}

// returns a data object from the bin
// k: the int in the bin of the T to return
template <class T>
T& Bin<T>::operator[](unsigned int k) const
{
	return a[k+offset];
}

// adds a new data object to the bin, returns the index
// where it was added to the bin
// &e reference to the T that is to be added to the bin
template<class T>
int Bin<T>::add(const T& e)
{
	if(nitems == room)
	{
		grow();
	}
	a[nitems] = e;
	return (nitems++)-offset;
}

// remove last num items from the Bin
// returns the index of the last element in the bin
//	num : number of items to remove
template<class T>
int Bin<T>::removeLast(int num)
{
	assert(num <= nitems);
	for(int i = nitems-num; i < nitems; i++)
	{
		a[i].~T();
	}
	nitems -= num;
	return nitems-offset;
}

template<class T>
int Bin<T>::removePos(int pos)
{
#ifdef SP2
	a[pos].~T();
#else
	delete a[pos];
#endif
	for(int i=pos+1; i < nitems; i++)
	{
		a[i-1] = a[i];
	}
	nitems--;
	return nitems-offset;
}

// returns a pointer to the data objects in the bin
template<class T>
T* Bin<T>::array() const
{
	return (T*)(&a[offset]);
}

template<class T>
void Bin<T>::setBlocksize(int blocksize)
{
	if(blocksize == 0)
	{
		block = (size_t(1) > size_t(4096/sizeof(T)) ?
				 size_t(1) : size_t(4096/sizeof(T)));
	}
	else
	{
		block = blocksize;
	}
}

template<class T>
void Bin<T>::setOffset(int off)
{
	offset = off;
}

// reallocates memory to exactly fit the data objects in the bin
// and returns a pointer to the data objects. The Bin is reinitialized
template<class T>
T* Bin<T>::done()
{
	T* tmp;
	if(nitems == 0)
	{
		free(a);
		tmp = 0;
	}
	else
	{
		tmp = (T*)std::realloc(a, nitems*sizeof(T));
		assert(tmp);
	}
	init(block);
	return tmp;
}


template<class T>
void Bin<T>::moveTo(Bin<T> &newbin)
// transfer the contents of this to newbin and cleanUp this
{
	memcpy(&newbin, this, sizeof(Bin<T>)); // do not use operator= 
	init(block);
}


template<class T>
void Bin<T>::grow()
{
	room += block;
	if(a == 0)
	{
		a = (T*)std::malloc(room*sizeof(T));
	}
	else
	{
		a = (T*)std::realloc(a, room*sizeof(T));
	}
	assert(a);
}

#endif

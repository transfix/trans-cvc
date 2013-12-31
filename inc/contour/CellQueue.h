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
// cellQueue.h - queue of cell identifiers.  The circular queue dyanmically
//               resizes itself when full.  Elements in the queue are of
//               type and size specified by the user.  Not even template'd
//               yet, just using a void *.

#ifndef CELL_QUEUE_H
#define CELL_QUEUE_H

//#include <Utility/utility.h>

class CellQueue
{
	public:
		// constructor/destructor
		CellQueue(int size=100);
		~CellQueue();
		// add item to the queue
		void Add(unsigned int cell);
		// remove and return the first item in queue
		int  Get(int& cell);
		// return the first item in queue
		int  Peek(int& cell);
		// remove the first item in queue
		void Pop();
		// reset to empty
		void Reset(void);
		// check if queue is empty
		int  Empty(void);
	private:
		int nel;
		int cellsize;  // # of elements in cell array
		int start;
		unsigned int* cells;
};

#endif

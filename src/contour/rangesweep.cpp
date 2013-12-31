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
#include <contour/rangesweep.h>
#include <contour/iqueue.h>

#define MAXQUEUESIZE
#ifdef MAXQUEUESIZE
static int maxqsize = 0;
#endif

extern int verbose;

class QueueRec
{
	public:
		QueueRec(int c=0)
		{
			init(c);
		}

		void init(int c)
		{
			cellid=c;
			range[0].MakeEmpty();
			range[1].MakeEmpty();
			range[2].MakeEmpty();
			range[3].MakeEmpty();
			range[4].MakeEmpty();
			range[5].MakeEmpty();
		}

		int operator ==(QueueRec& qr)
		{
			return(cellid == qr.cellid);
		}

		int cellid;
		Range fullrange;
		Range range[6];
};


void
rangeSweep::PropagateRegion(int cellid, float min, float max)
{
	static IndexedQueue<QueueRec, int> q;
	QueueRec qr, current, *old;
	RangeSweepRec rsr, *rsritem;
	u_int c;
	int adjc, cindex;
	qr.init(cellid);
	qr.fullrange.Set(min, max);
	q.enqueue(qr, cellid);
	while(!q.isEmpty())
	{
#ifdef MAXQUEUESIZE
		if(q.getLength() > maxqsize)
		{
			maxqsize = q.getLength();
			if(verbose)
				if(maxqsize%10000 == 0)
				{
					printf("qsize: %d\n", maxqsize);
				}
		}
#endif
		q.dequeue(current);
		// if cell is done, remove from heap and continue;
		if(plot.CellTouched(current.cellid))
		{
			rsr.cellid = current.cellid;
			if(queue.find(rsr.cellid) != NULL)
			{
				queue.remove(rsr.cellid);
			}
			continue;
		}
		// if cell not in heap, add it
		rsr.cellid = current.cellid;
		if((rsritem=queue.find(current.cellid)) == NULL)
		{
			data.getCellRange(rsr.cellid, min, max);
			rsr.range.Set(min, max);
			queue.insert(rsr, max-min, rsr.cellid);
			rsritem=queue.find(rsr.cellid);
		}
		if(rsritem->range.Disjoint(current.fullrange))
		{
			continue;
		}
		// subtract off the propagated range
		rsritem->range -= current.fullrange;
		// update priority, possibly removing item from queue
		if(rsritem->range.Empty())
		{
			queue.remove(current.cellid);
			plot.TouchCell(current.cellid);
		}
		else
		{
			queue.updatePriority(current.cellid, rsritem->range.MaxAll()-
								 rsritem->range.MinAll());
		}
		// don't use rsritem after this point.. may be deleted
		rsritem = NULL;
		// propagate ranges to shared faces
		for(c=0; c<data.getNCellFaces(); c++)
		{
			adjc = data.getCellAdj(current.cellid, c);
			if(adjc != -1 && !plot.CellTouched(adjc))
			{
				// get the range of the shared face
				data.getFaceRange(current.cellid, c, min, max);
				// find the index of this cell
				cindex = data.getAdjIndex(adjc, current.cellid);
				// propagate to this cell the intersection
				qr.init(adjc);
				qr.fullrange = current.fullrange;
				// need to take complement, but be careful not to remove
				// a constant cell, lest we not complete the propagation
				if(qr.fullrange.MinAll() != qr.fullrange.MaxAll())
				{
					qr.fullrange -= Range(-10000000, min);
					qr.fullrange -= Range(max, 10000000);
				}
				// don't propagate anything which came from this face
				qr.fullrange -= current.range[c];
				if((old = q.find(adjc)) != NULL)
				{
					// item already in queue
					old->fullrange += qr.fullrange;
					old->range[cindex] += qr.fullrange;
				}
				else if(!qr.fullrange.Empty())
				{
					qr.range[cindex] = qr.fullrange;
					q.enqueue(qr, adjc);
				}
			}
		}
	}
}

void rangeSweep::compSeeds(void)
{
	RangeSweepRec rsr, item;
	Range fullrange, added, outgoing;
	float min, max;
	// clear the array of mark bits
	plot.ClearTouched();
	seeds.Clear();
	// insert cell 0 into queue to begin
	rsr.cellid = 0;
	data.getCellRange(0, min, max);
	rsr.range.Set(min, max);
	queue.insert(rsr, max-min, rsr.cellid);
	// process queue of cells
	while(! queue.isEmpty())
	{
		// get the item
		queue.ipqmax(item);
		// cell is a seed cell
		seeds.AddSeed(item.cellid, item.range.MinAll(), item.range.MaxAll());
		// mark this cell as processed
		PropagateRegion(item.cellid, item.range.MinAll(), item.range.MaxAll());
	}
}

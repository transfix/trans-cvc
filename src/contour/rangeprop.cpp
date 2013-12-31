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
#include <string.h>
#include <contour/rangeprop.h>

extern int verbose;

void rangeProp::compSeeds(void)
{
	RangePropRec rpr, *item, *qitem;
	Range fullrange, added, outgoing;
	int index;
	float min, max;
	int adjc;
	int done;
	u_int c;
	if(verbose)
	{
		printf("------- computing seeds\n");
	}
	// clear the array of mark bits
	plot.ClearTouched();
	seeds.Clear();
	// insert cell 0 into queue to begin
	rpr.cellid = 0;
	data.getCellRange(0, min, max);
	rpr.resp.Set(min, max);
	rpr.comp.MakeEmpty();
	queue.enqueue(rpr);
	done = 0;
	// process queue of cells
	while(! queue.isEmpty())
	{
		// get the item
		item = queue.dequeue();
		done++;
		if(verbose)
			if(done % 1000 == 0)
			{
				printf("%g%% done\n", 100.0*done/(float)data.getNCells());
			}
		// mark this cell as processed
		plot.TouchCell(item->cellid);
		// compute the outgoing range (range which can be further propagated)
		outgoing.MakeEmpty();
		for(c=0; c<data.getNCellFaces(); c++)
		{
			adjc = data.getCellAdj(item->cellid, c);
			if(adjc != -1 && !plot.CellTouched(adjc))
			{
				// the range of the shared face may be propagated
				data.getFaceRange(item->cellid, c, min, max);
				outgoing += Range(min,max);
			}
		}
		// this is the full range of responsibility
		fullrange = (item->resp + outgoing) - item->comp;
		if(fullrange.Empty() || (!outgoing.Empty() &&
								 (outgoing.MinAll() <= fullrange.MinAll() &&
								  outgoing.MaxAll() >= fullrange.MaxAll())))
		{
			// propagate entire range
			for(c=0; c<data.getNCellFaces(); c++)
			{
				adjc = data.getCellAdj(item->cellid, c);
				if(adjc != -1 && !plot.CellTouched(adjc))
				{
					data.getFaceRange(item->cellid, c, min, max);
					// compute the range which should be propagated to this cell
					added = Range(min,max) - item->comp;
					item->comp += Range(min,max);
					rpr.cellid = adjc;
					if((index=queue.find(rpr)) != -1)
					{
						qitem=queue.getItem(index);
						qitem->resp += added;
						qitem->comp = Range(min,max) - added;
					}
					else
					{
						rpr.resp = added;
						rpr.comp = Range(min,max) - added;
						queue.enqueue(rpr);
					}
				}
			}
		}
		else
		{
			// cell is a seed cell
			seeds.AddSeed(item->cellid, fullrange.MinAll(), fullrange.MaxAll());
			// set range of faces to the _complement_ of adjacent cells
			for(c=0; c<data.getNCellFaces(); c++)
			{
				adjc = data.getCellAdj(item->cellid, c);
				if(adjc != -1 && !plot.CellTouched(adjc))
				{
					data.getFaceRange(item->cellid, c, min, max);
					rpr.cellid = adjc;
					if((index=queue.find(rpr)) != -1)
					{
						qitem=queue.getItem(index);
						qitem->comp.Set(min,max);
						// note : can we ADD to the comp here?
					}
					else
					{
						rpr.resp.MakeEmpty();
						rpr.comp.Set(min,max);
						queue.enqueue(rpr);
					}
				}
			}
		}
	}
	if(verbose)
	{
		printf("computed %d seeds\n", seeds.getNCells());
	}
}

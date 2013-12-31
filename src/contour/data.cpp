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
// data.C - class for scalar data
// Copyright (c) 1997 Dan Schikore

#include <stdio.h>
#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif

#include <string.h>
#include <contour/data.h>

int Data::funtopol1;
int Data::funtopol2;

// Data() - alternative constructor for the libcontour library
Data::Data(Data::DataType t, int _ndata, u_char* data)
{
	type     = t;
	ndata    = _ndata;
	filename = NULL;
	min = NULL;
	max = NULL;
	if(ndata > 1)
	{
		funcolor = 1;
		funcontour = 0;
		funtopol1 = 0;
		funtopol2 = 1;
	}
	else
	{
		funcontour = 0;
		funcolor = 0;
	}
}

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
// utilities.h -- some utility functions & macros
#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <contour/basic.h>
#include <fstream>
#include <limits.h>
#include <stdio.h>
#include <time.h>

// Issue a warning message.
void warning(const char* msg);

// Print an error message and exit.
void panic(const char* msg);

// Initialize random numebr generator (based on current clock).
void initrand(void);

// Print date, time, host.
void run_stamp(char* msg, int sz);

// read gzip'ed files
class gzifstream: public std::ifstream
{
	public:
		// Default c'tor
		gzifstream();
		/* Open file #name# for reading. If the file name ends with a ".gz"
		  extension, it is assumed to be a gzip'ed file and will be piped through zcat */
		gzifstream(const char* name);
		// D'tor
		~gzifstream();
		// Open file #name# for reading.
		void open(const char* name);

	private:
		int gzopen(const char* name);
		int gzclose();

		bool pipe;
		FILE* fp;
};

// timing
// Time elapsed between two instants of time.
double diffclock(clock_t t2, clock_t t1);

// brakes !!   (used to stop execution)
// When a <ctrl>-c is hit, the variable brakes is set to true.
extern bool brakes;

// Initialize #brake#.
void brake_init(void);

// Check Bool condition

// Check Boolean condition.
void check(int cond, const char* msg);

// Reset check counter.
void check_reset(void);

// Return number of failed checks.
int check_count(void);

#endif

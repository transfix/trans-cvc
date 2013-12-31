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
// contour.cpp - isocontouring library implementation
// Copyright (c) 1998 Emilio Camahort

#ifndef WIN32
#include <errno.h>
#endif
#include <time.h>
#include <contour/data.h>
#include <contour/Dataset.h>
#include <contour/Conplot.h>
#include <contour/contour.h>
#include <contour/conplot2d.h>
#include <contour/conplot3d.h>
#include <contour/conplotreg2.h>
#include <contour/conplotreg3.h>
#include <contour/datasetslc.h>
#include <contour/datasetvol.h>
#include <contour/datasetreg2.h>
#include <contour/datasetreg3.h>
//#include <Utility/utility.h>

// defaultHandler() - contour library default error handler
void defaultHandler(char* str, int fatal)
{
	if(fatal)
	{
		fprintf(stderr, "libcontour: fatal error: %s\n", str);
#ifndef WIN32
		if(errno)
#endif
			perror("libcontour");
		exit(0);
	}
	else
	{
		fprintf(stderr, "libcontour: error: %s\n", str);
#ifndef WIN32
		if(errno)
#endif
			perror("libcontour");
	}
}

// verbose, errorHandler - global variables
int verbose = 1;					// verbose level
void (*errorHandler)(char*, int) = defaultHandler;	// error handler

// setPreprocessing() - display the progress of preprocessing
void setPreprocessing(int percent, void* data)
{
	fprintf(stderr, "Preprocessing: %d done\n", percent);
}

// setVerboseLevel() - set the library's verbose level, default is 1
void setVerboseLevel(int level)
{
	verbose = level;
}

// setErrorHandler() - set the library's error handler
void setErrorHandler(void (*handler)(char*, int))
{
	errorHandler = handler;
}

// newDatasetUnstr() - create a new dataset for unstructured data
ConDataset* newDatasetUnstr(int datatype, int meshtype, int nvars,
			int ntime, int nverts, int ncells,
			double* verts, u_int* cells,
			int* celladj, u_char* data)
{
	int		t;		// a timestep index variable
	int		var;		// a variable index variable
	Data::DataType	type;		// data set parmeters
	Datasetslc*	slc;		// different dataset pointers
	Datasetvol*	vol;
	ConDataset*	dataset;
	type    = Data::DataType(datatype);	// miscellaneous initializations
	dataset = new ConDataset;
	dataset->vnames = NULL;
	switch(meshtype)		// create big objects: data and `plot'
	{
		case CONTOUR_2D :
			dataset->data = slc = new Datasetslc(type, nvars, ntime, nverts,
												 ncells, verts, cells, celladj, data);
			// create the plot
			dataset->plot = new Conplot2d(slc);
			break;
		case CONTOUR_3D :
			dataset->data = vol = new Datasetvol(type, nvars, ntime, nverts,
												 ncells, verts, cells, celladj, data);
			// create the plot
			dataset->plot = new Conplot3d(vol);
			break;
		default :
			errorHandler((char*)"newDatasetUnstr: incorrect mesh type", false);
			return NULL;
	};
	// allocate and init signature data
	dataset->sfun = new Signature** [dataset->data->nData()];
	for(var = 0; var < dataset->data->nData(); var++)
	{
		dataset->sfun[var] = new Signature*[dataset->data->nTime()];
		for(t = 0; t < dataset->data->nTime(); t++)
		{
			dataset->sfun[var][t] = NULL;
		}
	}
	// check for errors, return
	if(!dataset->data)
	{
		errorHandler((char*)"newDatasetUnstr: couldn't create dataset", false);
		return NULL;
	}
	if(!dataset->plot)
	{
		errorHandler((char*)"newDatasetUnstr: couldn't create plot", false);
		return NULL;
	}
	if(verbose)
	{
		printf("libcontour:newDatasetUnstr: data set created\n");
	}
	return dataset;
}

// newDatasetReg() - create a new dataset structure for a regular grid
ConDataset* newDatasetReg(int datatype, int meshtype, int nvars,
			int ntime, int* dim, u_char* data)
{
	int		t;			// a timestep index variable
	int		var;			// a variable index variable
	Data::DataType	type;		// data set parmeters
	Datasetreg2*	reg2;			// different dataset pointers
	Datasetreg3*	reg3;
	ConDataset*	dataset;
	type    = Data::DataType(datatype);	// miscellaneous initializations
	dataset = new ConDataset;
	dataset->vnames = NULL;
	switch(meshtype)			// create big objects: data and `plot'
	{
		case CONTOUR_REG_2D :
			dataset->data = reg2 = new Datasetreg2(type, nvars, ntime, dim, data);
			// create the plot
			dataset->plot = new Conplotreg2(reg2);
			break;
		case CONTOUR_REG_3D :
			dataset->data = reg3 = new Datasetreg3(type, nvars, ntime, dim, data);
			// create the plot
			dataset->plot = new Conplotreg3(reg3);
			break;
		default :
			errorHandler((char*)"newDatasetReg: incorrect mesh type", false);
			return NULL;
	};
	// allocate and init signature data
	dataset->sfun = new Signature** [dataset->data->nData()];
	for(var = 0; var < dataset->data->nData(); var++)
	{
		dataset->sfun[var] = new Signature*[dataset->data->nTime()];
		for(t = 0; t < dataset->data->nTime(); t++)
		{
			dataset->sfun[var][t] = NULL;
		}
	}
	// check for errors, return
	if(!dataset->data)
	{
		errorHandler((char*)"newDatasetReg: couldn't create dataset", false);
		return NULL;
	}
	if(!dataset->plot)
	{
		errorHandler((char*)"newDatasetReg: couldn't create plot", false);
		return NULL;
	}
	if(verbose)
	{
		printf("libcontour:newDatasetReg: data set created\n");
	}
	return dataset;
}

// loadDataset() - load data set from disk
ConDataset* loadDataset(int datatype, int meshtype, int nvars, int ntime, char** files)
{
	int		i;			// a variable name index variable
	int		t;			// a timestep index variable
	int		var;			// a variable index variable
	FILE*	vnfd;			// variable name file descriptor
	char	filename[256];		// variable name file name
	Data::DataType	type;		// data set parmeters
	Datasetslc*	slc;			// different dataset pointers
	Datasetvol*	vol;
	Datasetreg2*	reg2;
	Datasetreg3*	reg3;
	ConDataset*	dataset;
	type    = Data::DataType(datatype);	// miscellaneous initializations
	dataset = new ConDataset;
	switch(meshtype)			// create big objects: data and `plot'
	{
		case CONTOUR_2D :
			dataset->data = slc = new Datasetslc(type, nvars, ntime, files);
			// create the plot
			dataset->plot = new Conplot2d(slc);
			break;
		case CONTOUR_3D :
			dataset->data = vol = new Datasetvol(type, nvars, ntime, files);
			// create the plot
			dataset->plot = new Conplot3d(vol);
			break;
		case CONTOUR_REG_2D :
			dataset->data = reg2 = new Datasetreg2(type, nvars, ntime, files);
			// create the plot
			dataset->plot = new Conplotreg2(reg2);
			break;
		case CONTOUR_REG_3D :
			dataset->data = reg3 = new Datasetreg3(type, nvars, ntime, files);
			// create the plot
			dataset->plot = new Conplotreg3(reg3);
			break;
		default :
			errorHandler((char*)"loadDataset: incorrect mesh type", false);
			return NULL;
	};
	// read variable name file
	strcpy(filename, files[0]);
	strcpy(strrchr(filename, int('.')), ".var");
	vnfd = fopen(filename, "r");                // open file
	if(vnfd)
	{
		dataset->vnames = new char*[nvars]; // allocate memory to names
		for(i = 0; i < nvars; i++)
		{
			// read variable names
			dataset->vnames[i] = new char[80];
			fgets(dataset->vnames[i], 80, vnfd);
			dataset->vnames[i][strlen(dataset->vnames[i])-1] = 0;
		}
		fclose(vnfd);
	}
	else                                        // no variable names available
	{
		dataset->vnames = NULL;
	}
	// allocate and init signature data
	dataset->sfun = new Signature** [dataset->data->nData()];
	for(var = 0; var < dataset->data->nData(); var++)
	{
		dataset->sfun[var] = new Signature*[dataset->data->nTime()];
		for(t = 0; t < dataset->data->nTime(); t++)
		{
			dataset->sfun[var][t] = NULL;
		}
	}
	// check for errors, return
	if(!dataset->data)
	{
		errorHandler((char*)"loadDataset: couldn't create dataset", false);
		return NULL;
	}
	if(!dataset->plot)
	{
		errorHandler((char*)"loadDataset: couldn't create plot", false);
		return NULL;
	}
	if(verbose)
	{
		printf("libcontour:loadDataset: Data set loaded\n");
	}
	return dataset;
}

// getVariableNames() - get variable names for multi-variate data
char**	getVariableNames(ConDataset* dataset)
{
	if(!dataset)
	{
		errorHandler((char*)"getVariableNames: invalid dataset", false);
		return NULL;
	}
	return dataset->vnames;
}

// getDatasetInfo() - get the dataset's basic information
DatasetInfo*	getDatasetInfo(ConDataset* dataset)
{
	int		var;			// a variable index
	DatasetInfo*	datainfo;		// pointer to dataset's information
	if(!dataset || !dataset->data || !dataset->plot)
	{
		errorHandler((char*)"getDatasetInfo: invalid dataset", false);
		return NULL;
	}
	datainfo = new DatasetInfo;
	datainfo->datatype = int(dataset->data->dataType());
	datainfo->meshtype = dataset->data->meshType();
	datainfo->nvars = dataset->data->nData();
	datainfo->ntime = dataset->data->nTime();
	memset(datainfo->dim,  0, 3 * sizeof(u_int));
	memset(datainfo->orig, 0, 3 * sizeof(float));
	memset(datainfo->span, 0, 3 * sizeof(float));
	switch(dataset->data->meshType())		// not very orthodox, but ...
	{
		case CONTOUR_REG_2D :
			((Datareg2*)dataset->data->getData(0))->getDim(datainfo->dim);
			((Datareg2*)dataset->data->getData(0))->getOrig(datainfo->orig);
			((Datareg2*)dataset->data->getData(0))->getSpan(datainfo->span);
			break;
		case CONTOUR_REG_3D :
			((Datareg3*)dataset->data->getData(0))->getDim(datainfo->dim);
			((Datareg3*)dataset->data->getData(0))->getOrig(datainfo->orig);
			((Datareg3*)dataset->data->getData(0))->getSpan(datainfo->span);
			break;
		default :
			break;
	}
	dataset->data->getData(0)->getExtent(datainfo->minext, datainfo->maxext);
	datainfo->minvar = new float[dataset->data->nData()];
	datainfo->maxvar = new float[dataset->data->nData()];
	for(var = 0; var < dataset->data->nData(); var++)
	{
		datainfo->minvar[var] = dataset->data->getMinFun(var);
		datainfo->maxvar[var] = dataset->data->getMaxFun(var);
	}
	return datainfo;
}

// getSeedCells() - get seed cell data
SeedData* getSeedCells(ConDataset* dataset, int variable, int timestep)
{
	SeedData* seeddata; // pointer to seed data structure
	// sanity checks
	if(!dataset || !dataset->data || !dataset->plot)
	{
		errorHandler((char*)"getSeedCells: Couldn't find dataset", false);
		return NULL;
	}
	if(variable < 0 || variable >= dataset->data->nData())
	{
		errorHandler((char*)"getSeedCells: variable out of range", false);
		return NULL;
	}
	if(timestep < 0 || timestep >= dataset->data->nTime())
	{
		errorHandler((char*)"getSeedCells: timestep out of range", false);
		return NULL;
	}
	// extract seeds
	seeddata = new SeedData;
	dataset->data->getData(timestep)->setContourFun(variable);
	dataset->plot->setTime(timestep);
	// determine if seeds computed
	if(dataset->plot->getSeeds()->getNCells() == 0)
	{
		dataset->plot->Preprocess(timestep, setPreprocessing, NULL);
	}
	seeddata->nseeds = dataset->plot->getSeeds()->getNCells();
	seeddata->seeds  = (Seed*)dataset->plot->getSeeds()->getCellPointer();
	if(verbose > 1)
		for(int i = 0; i < seeddata->nseeds; i++)
			printf("seed cell %d --> min = %f max = %f  id = %d\n", i,
				   seeddata->seeds[i].min, seeddata->seeds[i].max,
				   seeddata->seeds[i].cell_id);
	if(verbose)
	{
		printf("libcontour:getSeedCells: seed data extracted\n");
	}
	return seeddata;
}

// getNumberOfSignatures() - get number of signature functions
int getNumberOfSignatures(ConDataset* dataset)
{
	if(!dataset)
	{
		errorHandler((char*)"getNumberOfSignatures: invalid dataset", false);
		return -1;
	}
	return dataset->data->getData(0)->getNFunctions();
}

// getSignatureFunctions() - get signature functions
Signature* getSignatureFunctions(ConDataset* dataset, int variable, int timestep)
{
	int t;		// a timestep index variable
	int fun;	// signature function index
	// sanity checks
	if(!dataset || !dataset->data || !dataset->plot)
	{
		errorHandler((char*)"getSignatureFunctions: Couldn't find dataset", false);
		return NULL;
	}
	if(variable < 0 || variable >= dataset->data->nData())
	{
		errorHandler((char*)"getSignatureFunctions: variable out of range", false);
		return NULL;
	}
	if(timestep < 0 || timestep >= dataset->data->nTime())
	{
		errorHandler((char*)"getSignatureFunctions: timestep out of range", false);
		return NULL;
	}
	// obtain signature functions
	dataset->data->getData(timestep)->setContourFun(variable);
	dataset->plot->setTime(timestep);
	// compute signature functions
	if(verbose)
	{
		printf("libcontour: computing signature functions ...\n");
	}
	dataset->nsfun = dataset->data->getData(0)->getNFunctions();
#ifdef COMPUTE_FUNCTIONS_FOR_ALL_TIMES
	for(t = 0; t < dataset->data->nTime(); t++)	// per time step t
	{
		dataset->sfun[variable][t] = new Signature[dataset->nsfun];
		for(fun = 0; fun < dataset->nsfun; fun++)
		{
			dataset->sfun[variable][t][fun].name =
				strdup(dataset->data->getData(0)->fName(fun));
			dataset->sfun[variable][t][fun].fy   =
				dataset->data->getData(t)->compFunction(fun,
						dataset->sfun[variable][t][fun].nval,
						&dataset->sfun[variable][t][fun].fx);
		}
	}
#endif /* of COMPUTE_FUNCTIONS_FOR_ALL_TIMES */
	if(!dataset->sfun[variable][timestep])	// have signatures already?
	{
		t = timestep;
		dataset->sfun[variable][t] = new Signature[dataset->nsfun];
		for(fun = 0; fun < dataset->nsfun; fun++)
		{
			dataset->sfun[variable][t][fun].name =
				strdup(dataset->data->getData(0)->fName(fun));
			dataset->sfun[variable][t][fun].fy   =
				dataset->data->getData(t)->compFunction(fun,
						dataset->sfun[variable][t][fun].nval,
						&dataset->sfun[variable][t][fun].fx);
		}
	}
	if(verbose)
	{
		printf("libcontour::getSignatureData: signature data computed \n");
	}
	return dataset->sfun[variable][timestep];
}

// getSignatureValues() - get signature values for isovalue
float*	getSignatureValues(ConDataset* dataset, int variable, int timestep, float isovalue)
{
	int t;	 // a timestep index variable
	int fun; // signature function index
	float*	svalues; // signature values
	// sanity checks
	if(!dataset || !dataset->data || !dataset->plot)
	{
		errorHandler((char*)"getSignatureValues: Couldn't find dataset", false);
		return NULL;
	}
	if(variable < 0 || variable >= dataset->data->nData())
	{
		errorHandler((char*)"getSignatureValues: variable out of range", false);
		return NULL;
	}
	if(timestep < 0 || timestep >= dataset->data->nTime())
	{
		errorHandler((char*)"getSignatureValues: timestep out of range", false);
		return NULL;
	}
	// obtain signature values
	dataset->data->getData(timestep)->setContourFun(variable);
	dataset->plot->setTime(timestep);
	dataset->nsfun = dataset->data->getData(0)->getNFunctions();
	// do we have signatures for timestep?
	if(!dataset->sfun[variable][timestep])
	{
		t = timestep;
		dataset->sfun[variable][t] = new Signature[dataset->nsfun];
		for(fun = 0; fun < dataset->nsfun; fun++)
		{
			dataset->sfun[variable][t][fun].name =
				strdup(dataset->data->getData(0)->fName(fun));
			dataset->sfun[variable][t][fun].fy   =
				dataset->data->getData(t)->compFunction(fun,
						dataset->sfun[variable][t][fun].nval,
						&dataset->sfun[variable][t][fun].fx);
		}
	}
	svalues = new float[dataset->nsfun];
	for(fun = 0; fun < dataset->nsfun; fun++)
	{
		int l, r, m; // binary search from SIoXtSpectrum.h
		m = 0;
		l = 0;
		r = dataset->sfun[variable][timestep][fun].nval;
		while(l < r)
		{
			m = (l + r) >> 1;
			if(isovalue < dataset->sfun[variable][timestep][fun].fx[m])
			{
				r = m - 1;
			}
			else
			{
				l = m + 1;
			}
		}
		svalues[fun] = dataset->sfun[variable][timestep][fun].fy[m];
		if(verbose > 1)
			printf("function %d %s\t --> %d values: (55, %f)\n", fun,
				   dataset->sfun[variable][timestep][fun].name,
				   dataset->sfun[variable][timestep][fun].nval,
				   dataset->sfun[variable][timestep][fun].fy[55]);
	}
	if(verbose)
	{
		printf("libcontour:getSignatureValues: signature values computed\n");
	}
	return svalues;
}

// getSlice() - extract a 2d slice from a 3d regular data grid
SliceData* getSlice(ConDataset* dataset, int variable, int timestep, char axis, u_int index)
{
	u_int	dim[3];		// dataset dimensions
	SliceData* slice;	// slice data
	Data::datatypes	buffer;	// buffer to hold actual slice
	// sanity checks
	if(!dataset || !dataset->data || !dataset->plot)
	{
		errorHandler((char*)"getSlice: Couldn't find dataset", false);
		return NULL;
	}
	if(dataset->data->meshType() != CONTOUR_REG_3D)
	{
		errorHandler((char*)"getSlice: invalid mesh type: must be 3D regular", false);
		return NULL;
	}
	if(variable < 0 || variable >= dataset->data->nData())
	{
		errorHandler((char*)"getSlice: variable out of range", false);
		return NULL;
	}
	if(timestep < 0 || timestep >= dataset->data->nTime())
	{
		errorHandler((char*)"getSlice: timestep out of range", false);
		return NULL;
	}
	if(axis != 'x' && axis != 'y' && axis != 'z')
	{
		errorHandler((char*)"getSlice: invalid slice axis", false);
		return NULL;
	}
	slice = new SliceData;
	((Datareg3*)dataset->data->getData(0))->getDim(dim);
	// check index range and
	switch(axis)				// determine width and height
	{
		case 'x' :
			if(index >= dim[0])
			{
				errorHandler((char*)"getSlice: x-index out of range", false);
				return NULL;
			}
			else
			{
				slice->width  = dim[1];
				slice->height = dim[2];
			}
			break;
		case 'y' :
			if(index >= dim[1])
			{
				errorHandler((char*)"getSlice: y-index out of range", false);
				return NULL;
			}
			else
			{
				slice->width  = dim[2];
				slice->height = dim[0];
			}
			break;
		case 'z' :
			if(index >= dim[2])
			{
				errorHandler((char*)"getSlice: z-index out of range", false);
				return NULL;
			}
			else
			{
				slice->width  = dim[0];
				slice->height = dim[1];
			}
			break;
	}
	dataset->data->getData(timestep)->setContourFun(variable);
	dataset->plot->setTime(timestep);
	slice->datatype = int(dataset->data->dataType());
	// allocate memory for slice
	switch(slice->datatype)
	{
		case CONTOUR_UCHAR :
			buffer.ucdata = new u_char[slice->width * slice->height];
			break;
		case CONTOUR_USHORT :
			buffer.usdata = new u_short[slice->width * slice->height];
			break;
		case CONTOUR_FLOAT :
			buffer.fdata = new float[slice->width * slice->height];
			break;
	}
	// extract slice from dataset
	if(((Datareg3*)dataset->data->getData(timestep))->getSlice(
				variable, axis, index, &buffer))
	{
		errorHandler((char*)"Datareg3::getSlice(): Couldn't extract slice", false);
		return NULL;
	}
	if(verbose)
	{
		printf("libcontour::extractSlice: slice %d along axis %c \n", index, axis);
	}
	switch(slice->datatype)		     // assign buffer to slice data
	{
		case CONTOUR_UCHAR :
			slice->ucdata = buffer.ucdata;
			break;
		case CONTOUR_USHORT :
			slice->usdata = buffer.usdata;
			break;
		case CONTOUR_FLOAT :
			slice->fdata = buffer.fdata;
			break;
	}
	if(verbose)
	{
		printf("libcontour::extractSlice: slice extracted\n");
	}
	return slice;
}

// getContour2d() - extract a 2d isocontour from a 2d data set
Contour2dData*	getContour2d(ConDataset* dataset, int variable, int timestep, float isovalue)
{
	Contour2d*	isocontour;	// new isocontour
	Contour2dData*	contour2d;	// 2d isocontour data structure
	// sanity checks
	if(!dataset || !dataset->data || !dataset->plot)
	{
		errorHandler((char*)"getContour2d: Couldn't find dataset", false);
		return NULL;
	}
	if(dataset->data->meshType() != CONTOUR_2D &&
			dataset->data->meshType() != CONTOUR_REG_2D)
	{
		errorHandler((char*)"getContour2d: invalid mesh type: must be 2D", false);
		return NULL;
	}
	if(variable < 0 || variable >= dataset->data->nData())
	{
		errorHandler((char*)"getContour2d: variable out of range", false);
		return NULL;
	}
	if(timestep < 0 || timestep >= dataset->data->nTime())
	{
		errorHandler((char*)"getContour2d: timestep out of range", false);
		return NULL;
	}
	dataset->data->getData(timestep)->setContourFun(variable);
	dataset->plot->setTime(timestep);
	contour2d = new Contour2dData;
	if(verbose)
	{
		printf("libcontour:getContour2d: isovalue = %f\n", isovalue);
	}
	// determine if seeds computed
	if(dataset->plot->getSeeds()->getNCells() == 0)
	{
		dataset->plot->Preprocess(timestep, setPreprocessing, NULL);
	}
	// extract isocontour
	dataset->plot->ResetAll();
	dataset->plot->Extract(isovalue);
	isocontour = dataset->plot->getContour2d();
	contour2d->nvert = isocontour->getNVert();
	contour2d->nedge = isocontour->getNEdge();
	contour2d->vert  = isocontour->vert;
	contour2d->edge  = isocontour->edge;
	if(verbose)
	{
		printf("libcontour:getContour2d: nr of vertices: %d\n", contour2d->nvert);
		printf("libcontour:getContour2d: nr of edges: %d\n", contour2d->nedge);
	}
	return contour2d;
}

// getContour3d() - extract a 3d isocontour from a 3d data set
Contour3dData*	getContour3d(ConDataset* dataset, int variable, int timestep, float isovalue, int colorvar)
{
	Contour3d*		isocontour;	// new isocontour
	Contour3dData*	contour3d;	// 3d isocontour data structure
	// sanity checks
	if(!dataset || !dataset->data || !dataset->plot)
	{
		errorHandler((char*)"getContour3d: Couldn't find dataset", false);
		return NULL;
	}
	if(dataset->data->meshType() != CONTOUR_3D &&
			dataset->data->meshType() != CONTOUR_REG_3D)
	{
		errorHandler((char*)"getContour3d: invalid mesh type: must be 3D", false);
		return NULL;
	}
	if(variable < 0 || variable >= dataset->data->nData())
	{
		errorHandler((char*)"getContour3d: variable out of range", false);
		return NULL;
	}
	if(colorvar != NO_COLOR_VARIABLE)
		if(colorvar < 0 || colorvar >= dataset->data->nData())
		{
			errorHandler((char*)"getContour3d: invalid color variable", false);
			return NULL;
		}
	if(timestep < 0 || timestep >= dataset->data->nTime())
	{
		errorHandler((char*)"getContour3d: timestep out of range", false);
		return NULL;
	}
	dataset->data->getData(timestep)->setContourFun(variable);
	dataset->data->getData(timestep)->setColorFun(colorvar);
	dataset->plot->setTime(timestep);
	contour3d = new Contour3dData;
	if(verbose)
	{
		printf("libcontour::getContour3d: isovalue = %f\n", isovalue);
	}
	// determine if seeds computed
	if(dataset->plot->getSeeds()->getNCells() == 0)
	{
		dataset->plot->Preprocess(timestep, setPreprocessing, NULL);
	}
	// extract isocontour
	dataset->plot->ResetAll();
	dataset->plot->Extract(isovalue);
	isocontour = dataset->plot->getContour3d();
	contour3d->nvert    = isocontour->getNVert();
	contour3d->ntri     = isocontour->getNTri();
	contour3d->vert     = isocontour->vert;
	contour3d->vnorm    = isocontour->vnorm;
	contour3d->vfun     = isocontour->vfun;
	contour3d->tri      = isocontour->tri;
	contour3d->colorvar = colorvar;
	contour3d->fmin     = isocontour->fmin;
	contour3d->fmax     = isocontour->fmax;
	return contour3d;
}

// saveContour2d() - extract a 2d isocontour and save it to a file
void saveContour2d(ConDataset* dataset, int variable, int timestep, float isovalue, char* filename)
{
	// sanity checks
	if(!dataset || !dataset->data || !dataset->plot)
	{
		errorHandler((char*)"saveContour2d: Couldn't find dataset", false);
		return;
	}
	if(dataset->data->meshType() != CONTOUR_2D &&
			dataset->data->meshType() != CONTOUR_REG_2D)
	{
		errorHandler((char*)"saveContour2d: invalid mesh type: must be 2D", false);
		return;
	}
	if(variable < 0 || variable >= dataset->data->nData())
	{
		errorHandler((char*)"saveContour2d: variable out of range", false);
		return;
	}
	if(timestep < 0 || timestep >= dataset->data->nTime())
	{
		errorHandler((char*)"saveContour2d: timestep out of range", false);
		return;
	}
	dataset->data->getData(timestep)->setContourFun(variable);
	dataset->plot->setTime(timestep);
	if(verbose)
	{
		printf("libcontour:saveContour2d: isovalue = %f\n", isovalue);
	}
	// determine if seeds computed
	if(dataset->plot->getSeeds()->getNCells() == 0)
	{
		dataset->plot->Preprocess(timestep, setPreprocessing, NULL);
	}
	// extract isocontour
	dataset->plot->ResetAll();
	dataset->plot->Extract(isovalue);
	// save contour to iPoly file
	if(dataset->plot->getContour2d()->write(filename))
	{
		char	str[256];
		sprintf(str, "saveContour2d: couldn't save to file: %s\n", filename);
		errorHandler(str, false);
	}
	else if(verbose)
	{
		fprintf(stderr, "libcontour:saveContour2d: saved to: %s\n", filename);
	}
}

// saveContour3d() - extract a 3d isocontour and save it to a file
void saveContour3d(ConDataset* dataset, int variable, int timestep, float isovalue, int colorvar, char* filename)
{
	// sanity checks
	if(!dataset || !dataset->data || !dataset->plot)
	{
		errorHandler((char*)"saveContour3d: Couldn't find dataset", false);
		return;
	}
	if(dataset->data->meshType() != CONTOUR_3D &&
			dataset->data->meshType() != CONTOUR_REG_3D)
	{
		errorHandler((char*)"saveContour3d: invalid mesh type: must be 3D", false);
		return;
	}
	if(variable < 0 || variable >= dataset->data->nData())
	{
		errorHandler((char*)"saveContour3d: variable out of range", false);
		return;
	}
	if(colorvar != NO_COLOR_VARIABLE)
		if(colorvar < 0 || colorvar >= dataset->data->nData())
		{
			errorHandler((char*)"saveContour3d: invalid color variable", false);
			return;
		}
	if(timestep < 0 || timestep >= dataset->data->nTime())
	{
		errorHandler((char*)"saveContour3d: timestep out of range", false);
		return;
	}
	dataset->data->getData(timestep)->setContourFun(variable);
	dataset->data->getData(timestep)->setColorFun(colorvar);
	dataset->plot->setTime(timestep);
	if(verbose)
	{
		printf("libcontour::saveContour3d: isovalue = %f\n", isovalue);
	}
	// determine if seeds computed
	if(dataset->plot->getSeeds()->getNCells() == 0)
	{
		dataset->plot->Preprocess(timestep, setPreprocessing, NULL);
	}
	// extract isocontour
	dataset->plot->ResetAll();
	dataset->plot->Extract(isovalue);
	// save contour to iPoly file
	if(dataset->plot->getContour3d()->write(filename))
	{
		char	str[256];
		sprintf(str, "saveContour3d: couldn't save to file: %s\n", filename);
		errorHandler(str, false);
	}
	else if(verbose)
	{
		fprintf(stderr, "libcontour:saveContour3d: saved to: %s\n", filename);
	}
}

// writeIsoComponents() - extract and write isocontour components to disk
void writeIsoComponents(ConDataset* dataset, int variable, int timestep, float isovalue, int colorvar, char* fprefix)
{
	// sanity checks
	if(!dataset || !dataset->data || !dataset->plot)
	{
		errorHandler((char*)"writeIsoComponents: Couldn't find dataset", false);
		return;
	}
	if(variable < 0 || variable >= dataset->data->nData())
	{
		errorHandler((char*)"writeIsoComponents: variable out of range", false);
		return;
	}
	if(colorvar != NO_COLOR_VARIABLE)
		if(colorvar < 0 || colorvar >= dataset->data->nData())
		{
			errorHandler((char*)"writeIsoComponents: invalid color variable", false);
			return;
		}
	if(timestep < 0 || timestep >= dataset->data->nTime())
	{
		errorHandler((char*)"writeIsoComponents: timestep out of range", false);
		return;
	}
	dataset->data->getData(timestep)->setContourFun(variable);
	dataset->data->getData(timestep)->setColorFun(colorvar);
	dataset->plot->setTime(timestep);
	if(verbose)
	{
		printf("libcontour::writeIsoComponents: isovalue = %f\n", isovalue);
	}
	// determine if seeds computed
	if(dataset->plot->getSeeds()->getNCells() == 0)
	{
		dataset->plot->Preprocess(timestep, setPreprocessing, NULL);
	}
	// extract isocontour
	dataset->plot->ResetAll();
	dataset->plot->BeginWrite(fprefix);
	dataset->plot->Extract(isovalue);
	dataset->plot->EndWrite();
	if(verbose)
	{
		printf("libcontour:writeIsoComponents: components saved\n");
	}
}

// clearDataset() - clear (remove) dataset from memory
void clearDataset(ConDataset* dataset)
{
	int t;	// timestep index variable
	int v;	// variable index variable
	if(dataset == NULL)
	{
		return;
	}
	if(dataset->data && dataset->plot)		// sanity check
	{
		for(v = 0; v < dataset->data->nData(); v++)	// delete signatures
		{
			for(t = 0; t < dataset->data->nTime(); t++)
				if(dataset->sfun[v][t])
				{
					delete [] dataset->sfun[v][t];
				}
			delete [] dataset->sfun[v];
		}
		delete [] dataset->sfun;
		delete dataset->data;			// delete data, set to NULL
		delete dataset->plot;
	}
}

/// \file util.c
/// \brief Utility functions.
//
//	Copyright (c) 2014, Christian B. Mendl
//	All rights reserved.
//	http://christian.mendl.net
//
//	This program is free software; you can redistribute it and/or
//	modify it under the terms of the Simplified BSD License
//	http://www.opensource.org/licenses/bsd-license.php
//
//	Reference:
//	  Christian B. Mendl, Francesc Malet, Paola Gori-Giorgi
//	  Wigner localization in quantum dots from Kohn-Sham density functional theory without symmetry breaking
//	  Physical Review B 89, 125106 (2014)
//	  (preprint http://arxiv.org/abs/1311.6011)
//_______________________________________________________________________________________________________________________
//

#include "util.h"
#include <math.h>
#include <stdio.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Euclidean norm of a vector
///
double Norm(const double *x, int n)
{
	int i;

	double maxElem = fabs(x[0]);
	for (i = 1; i < n; i++)
	{
		if (fabs(x[i]) > maxElem) maxElem = fabs(x[i]);
	}

	if (maxElem == 0) {
		return 0;
	}

	double r = 0;
	for (i = 0; i < n; i++)
	{
		r += square(x[i]/maxElem);
	}
	r = maxElem*sqrt(r);

	return r;
}



//_______________________________________________________________________________________________________________________
///
/// \brief Read 'n' doubles from file 'fileName', expecting the file size to be exactly n*sizeof(double)
///
int ReadDoubleData(const char *filename, double *data, int n)
{
	FILE *fd = fopen(filename, "rb");
	if (fd == NULL)
	{
		perror("'fopen' failed.");
		return -1;
	}

	// obtain the file size
	fseek(fd, 0 , SEEK_END);
	long filesize = ftell(fd);
	rewind(fd);
	// printf("file size: %d\n", filesize);
	if (filesize != (long)(n*sizeof(double)))
	{
		perror("File size does not match.");
		return -2;
	}

	// copy the file into the data array
	if (fread(data, sizeof(double), n, fd) != (size_t)n)
	{
		perror("'fread' failed.");
		return -3;
	}

	fclose(fd);

	return 0;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Write 'n' doubles to file 'fileName'
///
int WriteDoubleData(const char *filename, double *data, int n, bool append)
{
	const char *mode = append ? "ab" : "wb";

	FILE *fd = fopen(filename, mode);
	if (fd == NULL)
	{
		perror("'fopen' failed.");
		return -1;
	}

	// write data array to file
	if (fwrite(data, sizeof(double), n, fd) != (size_t)n)
	{
		perror("'fwrite' failed.");
		return -3;
	}

	fclose(fd);

	return 0;
}

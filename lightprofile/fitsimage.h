#include "fitsio.h"

#ifndef FITSIMAGE_H
#define FITSIMAGE_H

class Fitsimage
{
public:
	Fitsimage( char *filename, int iomode );
	~Fitsimage();
	void getDim( long naxis[] );
	void readSubset( long cpixel[], long subsize, double data[] );
	void writeSubset( long cpixel[], long subsize, double data[] );
	void printError();
private:
	fitsfile *fptr;
	int status;
	long naxes[2];
};

#endif
#include <iostream>
#include "fitsimage.h"
#include "fitsio.h"
using namespace std;


Fitsimage::Fitsimage( char *filename, int iomode )
{
	int nfound;

	status = 0;

	fits_open_file( &fptr, filename, iomode, &status );
	fits_read_keys_lng( fptr, "NAXIS", 1, 2, naxes, &nfound, &status );
	printError();
}

Fitsimage::~Fitsimage()
{
	fits_close_file( fptr, &status );
	printError();
}

void Fitsimage::getDim( long naxis[] )
{
	naxis[0] = naxes[0];
	naxis[1] = naxes[1];
}

void Fitsimage::readSubset( long cpixel[], long subsize, double data[] )
{
	int anynul = 0;
	double nulval = 0;
	bool overflow = 0;
	long fpixel[2], lpixel[2], inc[2] = {1, 1};

	fpixel[0] = cpixel[0] - subsize / 2; 
	fpixel[1] = cpixel[1] - subsize / 2;
	lpixel[0] = fpixel[0] + subsize - 1;
	lpixel[1] = fpixel[1] + subsize - 1;

	if ( fpixel[0] < 0 ) { fpixel[0] = 0; overflow = 1; }
	if ( fpixel[1] < 0 ) { fpixel[1] = 0; overflow = 1; }
	if ( lpixel[0] > naxes[0] ) { lpixel[0] = naxes[0]; overflow = 1; }
	if ( lpixel[1] > naxes[1] ) { lpixel[1] = naxes[1]; overflow = 1; }

	if ( overflow == 1 )
	{
		cerr << "Reading subset failed! An attempt putting a subset outside the image detected!" << endl;
		exit(1);
	}

	fits_read_subset( fptr, TDOUBLE, fpixel, lpixel, inc, &nulval, data, &anynul, &status );
}

void Fitsimage::writeSubset( long cpixel[], long subsize, double data[] )
{
	bool overflow = 0;
	long fpixel[2], lpixel[2];

	fpixel[0] = cpixel[0] - subsize / 2; 
	fpixel[1] = cpixel[1] - subsize / 2;
	lpixel[0] = fpixel[0] + subsize - 1;
	lpixel[1] = fpixel[1] + subsize - 1;

	if ( fpixel[0] < 0 ) { fpixel[0] = 0; overflow = 1; }
	if ( fpixel[1] < 0 ) { fpixel[1] = 0; overflow = 1; }
	if ( lpixel[0] > naxes[0] ) { lpixel[0] = naxes[0]; overflow = 1; }
	if ( lpixel[1] > naxes[1] ) { lpixel[1] = naxes[1]; overflow = 1; }

	if ( overflow )
	{
		cerr << "Writing subset failed! An attempt putting a subset outside the image detected!" << endl;
		exit(1);
	}

	fits_write_subset( fptr, TDOUBLE, fpixel, lpixel, data, &status );
}

void Fitsimage::printError()
{
	if (status)
	{
		fits_report_error(stderr, status);
		exit(status);
	}
}
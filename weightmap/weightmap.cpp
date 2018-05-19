#include <iostream>
#include <cmath>
#include "fitsio.h"
using namespace std;

int main(int argc, char* argv[])
{	
	fitsfile *fptr_read, *fptr_write;
	int status, nfound, anynull;
	long naxes[2], fpixel, npixels;
	float nullval;
	float *buffer;

	status = 0;
	fpixel = 1;	
	nullval = 0;

	fits_open_file( &fptr_read, argv[1], READONLY, &status);
	fits_read_keys_lng( fptr_read, "NAXIS", 1, 2, naxes, &nfound, &status);

	npixels = naxes[0] * naxes[1];
	buffer = new float [npixels];

	fits_read_img( fptr_read, TFLOAT, fpixel, npixels, &nullval, buffer, &anynull, &status);
	
	remove( argv[2]);
	fits_create_file( &fptr_write, argv[2], &status);
	fits_create_img( fptr_write, FLOAT_IMG, 2, naxes, &status);

	for ( int ii = 0 ; ii < npixels; ii ++ )
		buffer[ii] = 1 / sqrt( buffer[ii] );

	fits_write_img( fptr_write, TFLOAT, fpixel, npixels, buffer, &status);
	
	delete [] buffer;

	if ( fits_close_file( fptr_read, &status) )
		fits_report_error(stderr, status);
	if ( fits_close_file( fptr_write, &status) )
		fits_report_error(stderr, status);
	return (status);
}


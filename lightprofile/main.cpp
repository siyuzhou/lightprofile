#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include "fitsio.h"
#include "catalog.h"
#include "fitsimage.h"
#include "lightprofile.h"
#define XOBJ "X_IMAGE"
#define YOBJ "Y_IMAGE"
#define ROBJ "KRON_RADIUS"
#define AOBJ "A_IMAGE"
#define MAGOBJ "MAG_ISO"
#define pi 3.1415926
using namespace std;


void mag_minmax( int nvar, long nobjects, double data[], double *mag_min, double *mag_max )
{
	for ( long nobj = 0; nobj < nobjects; nobj ++ )
	{
		if ( data[nobj * nvar + nvar - 1] < *mag_min )
			*mag_min = data[nobj * nvar + nvar - 1];
		if ( data[nobj * nvar + nvar - 1] > *mag_max )
			*mag_max = data[nobj * nvar + nvar - 1];
	}
}

int overlay( int nvar, long nobjects, double data[], long subsize, long cpixel[] )
{
	long overlaysize;
	overlaysize = (subsize > 50 )? 50: subsize;
	for (long nobj = 0; nobj < nobjects; nobj ++ )
		if ( (cpixel[0] - data[nvar*nobj]) * (cpixel[0] - data[nvar*nobj]) + (cpixel[1] - data[nvar*nobj+1]) * (cpixel[1] - data[nvar*nobj+1])
			 < (data[nvar*nobj+2] * data[nvar*nobj+3] + overlaysize / 2.) * (data[nvar*nobj+2] * data[nobj*nvar+3] + overlaysize / 2.) )
			return 1;
	return 0;
}

int indeadzones( int ndeadzones, double deadzones[], long subsize, long cpixel[] )
{
	for (int n = 0; n < ndeadzones; n ++ )
		if ( overlay( 4, ndeadzones, deadzones, subsize, cpixel ) )
			return 1;
	return 0;
}

void readdeadzones( char *deadname, int *pndeadzones, double deadzones[] )
{
	string header;

	ifstream rdead( deadname );

	if ( rdead.fail() )
	{
		cerr << "Reading \"dead zones\" file failed." << endl;
		exit(1);
	}

	*pndeadzones = 0;
	
	getline( rdead, header );
	while ( rdead >> deadzones[*pndeadzones*4] >> deadzones[*pndeadzones*4+1] >> deadzones[*pndeadzones*4+2] >> deadzones[*pndeadzones*4+3] && *pndeadzones < 20 )
		++ (*pndeadzones);

	if ( !(*pndeadzones < 20) )
	{
		cerr << "Max number of dead zones is 20" << endl;
		exit(1);
	}

	rdead.close();
}

void addsub( long subsize, double sub1[], double sub2[] )
{
	for (long i = 0; i < subsize * subsize; i ++ )
		sub1[i] += sub2[i];
}

int main( int argc, char *argv[] )
{
	double *vardata;
	long nobjects;
	double rand_mag, rand_re;
	double *subimage1, *subimage2;
	long subsize;
	long cpixel[2], dim[2];
	double deadzones[80];
	int ndeadzones;
	double *coordarts;

	const int nvar = 5;
	string varname[nvar] = { XOBJ, YOBJ, ROBJ, AOBJ, MAGOBJ };
	string comment;
	double SERSIC_N, DA, RE_MIN, RE_MAX, PIXSCALE, MAG_MIN, MAG_MAX, MAG_ZERO, SUBRATIO;
	int NARTOBJECTS;

	ifstream readpara(argv[1]);
		getline( readpara, comment );
		readpara >> comment >> SERSIC_N;
		readpara >> comment >> DA;
		readpara >> comment >> RE_MIN;
		readpara >> comment >> RE_MAX;
		readpara >> comment >> PIXSCALE;
		readpara >> comment >> MAG_ZERO;
		readpara >> comment >> MAG_MIN;
		readpara >> comment >> MAG_MAX;
		readpara >> comment >> SUBRATIO;
		readpara >> comment >> NARTOBJECTS;
	readpara.close();

	
	coordarts = new double [4*NARTOBJECTS];	

	srand(time(NULL));

	Catalog drz( argv[2] );										/* Open catalog file */	
	
	nobjects = drz.getNumberObjs();								//Get number of objects in the catalog

	vardata = new double [ nobjects * nvar ];
	drz.readData( nvar, varname, vardata );						//Read data from catalog

	Fitsimage drz_fits( argv[3], READWRITE );					//Open fits file for reading and writing

	drz_fits.getDim( dim );
																

	Galaxy art_gal( MAG_ZERO, PIXSCALE );						//Create galaxy light profiles and write it on to the image
	
	ofstream art_record( argv[4] );								/* Open a file to record the positions 
																and magnitudes of artificial objects */
	art_record << '#' << setw(8) << "X_IMAGE" << setw(8) << "Y_IMAGE" << setw(8) << "RE_RAND" << setw(10) << "MAG_RAND" << setw(10) << "MAG_TOTAL" << setw(10) << "MAG_ERR" << endl;
	
	ndeadzones = 0;
	if ( argc == 5 )
		readdeadzones( argv[5], &ndeadzones, deadzones );

	for ( int nn = 0; nn < NARTOBJECTS; nn ++ )
	{
		rand_mag = rand() % 1000 / 1000. * ( MAG_MAX - MAG_MIN ) + MAG_MIN;	//Random magnitude
		rand_re = rand() % 100 / 100. * (RE_MAX - RE_MIN ) + RE_MIN;

		art_gal.makeProfile( SERSIC_N, rand_re, DA, rand_mag );	//Make a light profile according to re and mag

		subsize = 2 * (long) ( rand_re * SUBRATIO ) + 1;		//Assign subimage size
		subimage1 = new double [subsize * subsize];
		subimage2 = new double [subsize * subsize];
			
		art_gal.createSubimage( subsize, subimage1 );			//Create a subimage according to the light profile

		do														//Random position for the subimage
		{
			cpixel[0] = rand() % ( dim[0] - 2 * subsize ) + subsize;
			cpixel[1] = rand() % ( dim[1] - 2 * subsize ) + subsize;

			coordarts[4*nn] = cpixel[0];						//Record the center of the artificial object
			coordarts[4*nn+1] = cpixel[1];					
			coordarts[4*nn+2] = subsize / 2;					//Size of the artificial object is the subimage size
			coordarts[4*nn+3] = 1;								//Size modifier
		}
		while( overlay( nvar, nobjects, vardata, subsize, cpixel ) || overlay( 4, nn, coordarts, subsize, cpixel ) || indeadzones( ndeadzones, deadzones, subsize, cpixel ) ); 

		drz_fits.readSubset( cpixel, subsize, subimage2 );		//Read background where the subimage are to be placed
			
		addsub( subsize, subimage1, subimage2 );				//Add background to the subimage and write it to the sci image
		drz_fits.writeSubset( cpixel, subsize, subimage1 );

		delete [] subimage1;
		delete [] subimage2;
		
		art_record << setw(8) << cpixel[0] << setw(8) << cpixel[1] << setw(8) << rand_re << setw(10) << rand_mag << setw(10) << art_gal.getMag() << setw(10) << art_gal.getMag() - rand_mag << endl;
																//Record the position and magnitude of the artificial object
		
	}

	delete [] vardata;
	art_record.close();

	return 0;
}


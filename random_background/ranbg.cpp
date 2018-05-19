#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <string>
#include <iomanip>
#include "fitsio.h"
#define BACKMAX 7.0				//maximum reach of background in number of times of ROBJ
#define BACKMIN 1.2				//minimum reach of background in number of times of ROBJ
#define BACKSIZE 1				//background size in unit of object radius
#define NBACK 100				//number of apertures for background calculation
#define MAXNNEIGHBOUR 1000
#define XOBJ "X_IMAGE"
#define YOBJ "Y_IMAGE"
#define ROBJ "KRON_RADIUS"
#define AOBJ "A_IMAGE"
using namespace std;

int get_neighbours( char *catalogname, long objectid, float neighbour_axes[] );
int get_background( char *filename, int nneighbours, float neighbour_axes[], int backpixwidth, long backarray[NBACK][2][2], float backvalue[] );
float cal_background( float *avgback, float backvalue[] );
int check_image( char *fileread, char *filewrite, int nneighbours, float neighbour_axes[], int backpixwidth, long backarray[NBACK][2][2] );
int check_background( char *fileread_back, char *fileread_rms, long objx, long objy, float *check_back, float *check_rms );

int main( int argc, char *argv[] )
{
	long objectid, backarray[NBACK][2][2];
	int nneighbours, backpixwidth;
	float neighbour_axes[4*MAXNNEIGHBOUR + 4], backvalue[NBACK]={0};
	float back_avg, back_rms, check_back, check_rms;

	if ( !( argc == 5 || argc == 7 ) )
	{
		cerr << "4(or 6) arguments expected. [Image name] [Catalog name] [Object ID list] [Output file]"
			 << "([Background image name] [Background rms image name])."<< endl;
		exit(1);
	}

	if ( atol(argv[3]) == 0 )
	{
		ifstream objectlist(argv[3]);
		
		remove(argv[4]);
		ofstream outputfile(argv[4]);

		outputfile << '#' << setiosflags(ios::left) << setw(4) << 1 << setw(15) << "OBJECTID" << '\n'
				   << '#' << setw(4) << 2 << setw(15) << "NNEIGHBOURS" << "\tNumber of neighbours\n"
				   << '#' << setw(4) << 3 << setw(15) << "BACK_AVG" << "\tAverage background\n"
				   << '#' << setw(4) << 4 << setw(15) << "BACK_RMS" << "\tBackground RMS\n"
				   << '#' << setw(4) << 5 << setw(15) << "SE_BACK" << "\tBackground given by SE\n"
				   << '#' << setw(4) << 6 << setw(15) << "SE_RMS" << "\tBackground RMS given by SE" << endl;
		while ( objectlist >> objectid )
		{
			nneighbours = get_neighbours( argv[2], objectid, neighbour_axes );
	
			backpixwidth = ((int)neighbour_axes[3] * BACKSIZE ) / 2 * 2;     //aperture of background scale with size of object

			get_background( argv[1], nneighbours, neighbour_axes, backpixwidth, backarray, backvalue );
			back_rms = cal_background( &back_avg, backvalue );
			outputfile << setiosflags(ios::right) << setw(5) << objectid << setw(8) << nneighbours << setw(15) << back_avg << setw(15) << back_rms;	

//			check_image( argv[1], /*name*/, nneighbours, neighbour_axes, backpixwidth, backarray );

			if ( argc == 7 )
			{
				check_background( argv[5], argv[6], neighbour_axes[1], neighbour_axes[2], &check_back, &check_rms );
				outputfile << setw(15) << check_back << setw(15) << check_rms;
			}	
			outputfile << endl;
			cout << objectid << " finished" << endl;
		}
		
		objectlist.close();
		outputfile.close();
	}
	return 0;
}

int get_neighbours( char *catalogname, long objectid, float neighbour_axes[] )
{
	char hh;
	long nx, ny, nr, na, nn, objectn;
	float *databuffer, rdis, rx, ry, rback;
	int nneighbours;
	string label,comment;

	ifstream catalog_read( catalogname );
	if ( catalog_read.fail() )
	{	
		cerr << "Opening catalog failed." << endl;
		exit(1);
	}
	
	nx = 0;
	ny = 0;
	nr = 0;
	na = 0;
	nn = 0;

	catalog_read >> hh;    //read header and find columes for x, y and r of all objects
	while ( hh == '#' )
	{
		catalog_read >> nn >> label;
		getline( catalog_read, comment );

		if      ( label == XOBJ ) nx = nn;
		else if ( label == YOBJ ) ny = nn;
		else if ( label == ROBJ ) nr = nn;
		else if ( label == AOBJ ) na = nn;
			
		catalog_read >> hh;
	}
	catalog_read.seekg( -1, ios::cur ); //set file pointer ready to read data
	
	databuffer = new float [nn];

	/*looking for the target object*/
	catalog_read >> objectn;
	while ( (objectn <= objectid) && !(catalog_read.eof()) )
	{
		databuffer[0] = objectn;		
		for (int ii = 1; ii < nn; ii++)
			catalog_read >> databuffer[ii];
		catalog_read >> objectn;	
	}
	if ( objectn < objectid )
	{
		cerr << "Object ID not found!" << endl;
		exit(1);
	}

	neighbour_axes[0] = objectid;
	neighbour_axes[1] = databuffer[nx-1];
	neighbour_axes[2] = databuffer[ny-1];
	neighbour_axes[3] = databuffer[nr-1] * databuffer[na-1];
	
	rback = neighbour_axes[3] * BACKMAX;
	
	/*looking for neighbours*/
	catalog_read.seekg( 0, ios::beg );
	for ( int ii = 0; ii < nn; ii++ )		//skip header
		getline( catalog_read, comment );

	nneighbours = 0;

	catalog_read >> objectn;
	while ( !(catalog_read.eof()) )
	{
		databuffer[0] = objectn;

		for (int ii = 1; ii < nn; ii++)
			catalog_read >> databuffer[ii];
		
		rx = databuffer[nx-1] - neighbour_axes[1];
		ry = databuffer[ny-1] - neighbour_axes[2];
		rdis = sqrt( rx * rx + ry * ry ) - databuffer[nr-1] * databuffer[na-1]; //databuffer[nr-1] is the kron radius of the object
		
		if ( rdis < rback )				//found a neighbour
		{
			++nneighbours;
			if ( nneighbours > MAXNNEIGHBOUR )
			{
				cerr << "Too many neighbours" << endl;
				exit(1);
			}
			neighbour_axes[4*nneighbours] = databuffer[0];
			neighbour_axes[4*nneighbours + 1] = databuffer[nx-1];
			neighbour_axes[4*nneighbours + 2] =	databuffer[ny-1];
			neighbour_axes[4*nneighbours + 3] = databuffer[nr-1] * databuffer[na-1];
		}
		catalog_read >> objectn;
	}
	delete [] databuffer;			
	catalog_read.close();
	
	return nneighbours;
}

int get_background( char *filename, int nneighbours, float neighbour_axes[], int backpixwidth, long backarray[NBACK][2][2], float backvalue[] )
{
	fitsfile *fptr_read;
	int status, anynul, naper;	
	long xx, yy, fpixel[2], lpixel[2], inc[2]={1,1};
	float rr, angle;
	float dis, safedis, nulval, *backaper, apersum;
	bool overlap;
	status = 0;
	nulval = 0;

	backaper = new float [backpixwidth * backpixwidth];

	if ( fits_open_file( &fptr_read, filename, READONLY, &status ) )
	{
		fits_report_error( stderr, status );
		cerr << "File Open failed" << endl;
		exit(1);
	}
	naper = 0;

	srand( time(NULL) );
	while ( naper < NBACK )
	{
		rr = rand() / (float)(RAND_MAX) * (BACKMAX - BACKMIN) + BACKMIN;
		angle = rand() % 360;
	
		xx = (long) (cos( angle ) * (rr * neighbour_axes[3] - 0.75 * backpixwidth) + neighbour_axes[1]);
		yy = (long) (sin( angle ) * (rr * neighbour_axes[3] - 0.75 * backpixwidth) + neighbour_axes[2]);
		
		overlap = 0;
		for ( int ii = 1; ii < nneighbours+1; ii++ )		//check overlap of background aperture with neighbours
		{
			dis = sqrt( (xx - neighbour_axes[4*ii+1]) * (xx - neighbour_axes[4*ii+1])
					  + (yy - neighbour_axes[4*ii+2]) * (yy - neighbour_axes[4*ii+2]) );
			safedis = neighbour_axes[4*ii+3] * BACKMIN + 0.75 * backpixwidth;     //safe distance between background apertures and neighbours
			if ( dis < safedis )
			{
				overlap = 1;
				break;
			}
		} 
		if (!(overlap))
		{
			fpixel[0] = xx - backpixwidth / 2 + 1;				//whether background aperture is on an eage of the image is not checked
			fpixel[1] = yy - backpixwidth / 2 + 1;				//"backpixwidth" must be an even number
			lpixel[0] = xx + backpixwidth / 2;
			lpixel[1] = yy + backpixwidth / 2;

			backarray[naper][0][0] = fpixel[0];				//store background apertures' cooridinates
			backarray[naper][0][1] = fpixel[1];
			backarray[naper][1][0] = lpixel[0];
			backarray[naper][1][1] = lpixel[1];

			apersum = 0;	
			fits_read_subset( fptr_read, TFLOAT, fpixel, lpixel, inc, &nulval, backaper, &anynul, &status );
			for ( int jj = 0; jj < backpixwidth * backpixwidth; jj++ )
				apersum += backaper[jj];
			backvalue[naper] = apersum / (backpixwidth * backpixwidth);						//record background value, can be used along with backarray
			
			++naper;
		}
	}

	fits_close_file( fptr_read, &status );
	delete [] backaper;

	return status;
}

int check_image( char *fileread, char *filewrite, int nneighbours, float neighbour_axes[], int backpixwidth, long backarray[NBACK][2][2] )
{
	fitsfile *fptr_read, *fptr_write;
	int status, anynul;
	float *backmask;
	float *object_section, nulval;
	long imagesize, naxis[2], fpixel[2], lpixel[2], inc[2]={1,1};
	
	status = 0;
	fits_open_file( &fptr_read, fileread, READONLY, &status );    //copy object section
	if ( status )	 
	{
		fits_report_error( stderr, status );
		exit(1);
	}	
	naxis[0] = 2 * ((long)neighbour_axes[3] * BACKMAX);
	naxis[1] = 2 * ((long)neighbour_axes[3] * BACKMAX);
	imagesize = naxis[0] * naxis[1];							  //image size

	object_section = new float [ imagesize ];
	fpixel[0] = neighbour_axes[1] - naxis[0] / 2 + 1;			  //find the section with object in the center
	fpixel[1] = neighbour_axes[2] - naxis[1] / 2 + 1;
	lpixel[0] = neighbour_axes[1] + naxis[0] / 2;
	lpixel[1] = neighbour_axes[2] + naxis[1] / 2;

	fits_read_subset( fptr_read, TFLOAT, fpixel, lpixel, inc, &nulval, object_section, &anynul, &status );
	fits_close_file( fptr_read, &status );

	remove(filewrite);
	fits_create_file( &fptr_write, filewrite, &status );		  //copy object section	
	fits_create_img( fptr_write, FLOAT_IMG, 2, naxis, &status );

	fits_write_img( fptr_write, TFLOAT, 1, imagesize, object_section, &status );
	
	backmask = new float [backpixwidth * backpixwidth];
	for ( int ii = 0; ii < backpixwidth*backpixwidth; ii++ )              //create uniformed background masks
		backmask[ii] = 0.3;
	for ( int ii = 0; ii < NBACK; ii++ )						  //write background masks into check image
	{
		backarray[ii][0][0] = backarray[ii][0][0] - fpixel[0] + 1;
		backarray[ii][0][1] = backarray[ii][0][1] - fpixel[1] + 1;
		backarray[ii][1][0] = backarray[ii][1][0] - fpixel[0] + 1;
		backarray[ii][1][1] = backarray[ii][1][1] - fpixel[1] + 1;
		fits_write_subset( fptr_write, TFLOAT, backarray[ii][0], backarray[ii][1], backmask, &status );
	}

	fits_close_file( fptr_write, &status );
	delete [] backmask;
	delete [] object_section;
	return status;
}

float cal_background( float *avgback, float backvalue[] )
{
	float backsum, rms;

	backsum = 0.;	
	for ( int ii = 0; ii < NBACK; ii++ )
		backsum += backvalue[ii];
	*avgback = backsum / NBACK;

	backsum = 0.;	
	for ( int jj = 0; jj < NBACK; jj++ )
		backsum += ( backvalue[jj] - *avgback ) * ( backvalue[jj] - *avgback );
	rms = sqrt( 1./(NBACK-1) * backsum );

	return rms;
}

int check_background( char *fileread_back, char *fileread_rms, long objx, long objy, float *check_back, float *check_rms )
{
	fitsfile *fptr_read_back, *fptr_read_rms;
	int status, anynul;
	long fpixel[2];
	float nulval;

	fpixel[0] = objx;
	fpixel[1] = objy;
	status = 0;
	nulval = 0.;

	fits_open_file( &fptr_read_back, fileread_back, READONLY, &status );
	fits_open_file( &fptr_read_rms, fileread_rms, READONLY, &status );
	
	fits_read_pix( fptr_read_back, TFLOAT, fpixel, 1, &nulval, check_back, &anynul, &status );
	fits_read_pix( fptr_read_rms, TFLOAT, fpixel, 1, &nulval, check_rms, &anynul, &status );

	fits_close_file( fptr_read_rms, &status );
	fits_close_file( fptr_read_back, &status );
	
	return status;
}

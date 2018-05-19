#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include "catalog.h"
#define MAG "MAG_ISO"
using namespace std;

int sameobject( double x1, double y1, double x2, double y2, double r )
{
	return ( (x1 < x2+r/2) && (x1 > x2-r/2) && (y1 < y2+r/2) && (y1 > y2-r/2) );
}

int main( int argc, char *argv[] )
{
	double *artdata, temp;
	double *objdata;
	int nobj, narts;
	int artcount;
	bool found;

	string varname[5] = {"X_IMAGE", "Y_IMAGE", "KRON_RADIUS", "A_IMAGE", MAG};
	string comment;
	
	Catalog sci_art( argv[1] );

	nobj = sci_art.getNumberObjs();
	objdata = new double [5*nobj];
	sci_art.readData( 5, varname, objdata );
	
	ifstream readart( argv[2] );
	getline( readart, comment );

	narts = atol( argv[3] );
	artdata = new double [5*narts];
	
	for ( int i = 0; i < narts; i ++ )
	{
		for ( int j = 0; j < 5; j ++ )
			readart >> artdata[5*i+j];
		getline( readart, comment );
	}
	readart.close();

	ofstream writeart( argv[4] );
	writeart << '#' << setw(8) << "X_RAND" << setw(8) << "Y_RAND" << setw(10) << "RE_RAND" << setw(10) << "MAG_RAND" << setw(10) << "MAG_TOTAL" 
			 << setw(10) << "X_SE" << setw(10) << "Y_SE" << setw(10) << "KRON_SE" << setw(10) << "MAG_SE" << setw(15) << "MAG_SE_ERR" 
			 << setw(8) << "FOUND" << endl;

	artcount = 0;
	for ( int art = 0; art < narts; art ++ )
	{
		writeart << setw(8) << artdata[5*art] << setw(8) << artdata[5*art+1] << setw(10) << artdata[5*art+2] << setw(10) << artdata[5*art+3] << setw(10) << artdata[5*art+4];
		found = 0;
		for ( int obj = 0; obj < nobj; obj ++ )
			if ( sameobject( objdata[5*obj], objdata[5*obj+1], artdata[5*art], artdata[5*art+1], artdata[5*art+2] ) )
			{
				found = 1;
				writeart << setw(10) << objdata[5*obj] << setw(10) << objdata[5*obj+1] << setw(10) << objdata[5*obj+2]*objdata[5*obj+3] << setw(10) << objdata[5*obj+4] << setw(15) << objdata[5*obj+4] - artdata[5*art+4]
						 << setw(8) << found << endl;
				artcount ++;
				break;
			}
		
		if ( !found )
			writeart << setw(10) << -1 << setw(10) << -1 << setw(10) << -1 << setw(10) << -1 << setw(15) << -1
					 << setw(8) << found << endl;
	}

	cout << "Artificial objects extracted " << artcount << endl;
	writeart.close();

	delete [] artdata;

	return 0;
}

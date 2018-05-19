#include "catalog.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

Catalog::Catalog( char* cat_file )
{
	char hh;
	int nn;
	string label[30], comment;
	
	catRead.open( cat_file );
	if ( catRead.fail() )
	{
		cerr << "Opening catalog failed!" << endl;
		exit(1);
	}
			
	catRead >> hh;
	while ( hh == '#' )
	{
		catRead >> nn;
		catRead >> label[nn-1];
		getline( catRead, comment );
		
		catRead >> hh;
	}
	
	numberCols = nn;							//Initialise : numberCols
	nameCols = new string [numberCols];
	
	for ( nn = 0; nn < numberCols; nn ++ )
		nameCols[nn] = label[nn];				//Initialise : nameCols
	
	catRead.seekg( -1, ios::cur );
	
	numberObjs = 0;
	while ( !catRead.eof() )
	{
		getline( catRead, comment );
		numberObjs ++;
	}											//Find out number of Objects : numberObjects
	numberObjs --;
	catRead.clear();							//Clear eof() status
}

void Catalog::readData( const int nvar, string varname[], double vardata[] )
{
	int *varCol;
	double *databuffer;
	string comment;

	varCol = new int [nvar];
	databuffer = new double [numberCols];

	for ( int var = 0; var < nvar; var ++ )
	{
		for ( int col = 0; col < numberCols; col ++ )
			if ( varname[var] == nameCols[col] )
			{
				varCol[var] = col;		//Locate colume number of desired variable
				break;
			}
		if ( varCol[var] == 0 )
		{
			cerr << "Column \"" << varname[var] << "\" not found!" << endl;
			exit(1);
		}
	}
	
	catRead.seekg( 0, ios::beg );
	for ( int hh = 0; hh < numberCols; hh ++ )
		getline( catRead, comment );
	
	for ( int obj = 0; obj < numberObjs; obj ++ )
	{
		for ( int col = 0; col < numberCols; col ++ )
			catRead >> databuffer[col];
		for ( int var = 0; var < nvar; var ++ )
			vardata[obj*nvar + var] = databuffer[varCol[var]];
	}
	
	delete [] varCol;
	delete [] databuffer;

	catRead.clear();
}

long Catalog::getNumberObjs()
{
	return numberObjs;	
}

int Catalog::getNumberCols()
{
	return numberCols;
}

Catalog::~Catalog()
{
	catRead.close();
	delete [] nameCols;
}
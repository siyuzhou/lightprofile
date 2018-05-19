#include <string>
#include <fstream>
using namespace std;

#ifndef CATALOG_H
#define CATALOG_H

class Catalog
{
public:
	Catalog( char* cat_file );
	void readData( const int nvar, string varname[], double vardata[]);
	long getNumberObjs();
	int getNumberCols();
	~Catalog();
private:
	ifstream catRead;
	int numberCols;
	string *nameCols;
	long numberObjs;
};

#endif
#include <cmath>

#ifndef LIGHTPROFILE_H
#define LIGHTPROFILE_H

class Galaxy
{
public:
	Galaxy( double zeropoint, double pixscale );
	
	void makeProfile( double n, double re, double da, double mag );
	void createSubimage( long subsize, double subdata[] );
	double getMag();
	double getAlpha();
private:
	double alpha;
	double nSersic;
	double I0;
	double totalMag;
	const double zeroMag, pixScale;
	double profile( double rr );
};

#endif
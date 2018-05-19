#include "lightprofile.h"
//#include "Gamma.h"
#define pi 3.1415926

Galaxy::Galaxy( double zeropoint, double pixscale ): zeroMag( zeropoint ), pixScale( pixscale )
{
	nSersic = 0;				//Initialize a galaxy
	I0 = 0;
	totalMag = 99;
}

void Galaxy::makeProfile( double n, double re, double da, double mag )
{
	double F;

	nSersic = n;
	alpha = re / da * 180 / pi * 3600 / pixScale / 3459;
	F = pow( 10, (zeroMag - mag) * 0.4 );
	I0 = F / (2 * pi * (alpha * alpha) * nSersic * tgamma(2 * nSersic));
}

void Galaxy::createSubimage( long subsize, double subdata[] )
{
	double rr, xx, yy;
	double totalflux, pixflux;
	int nbins;
	
	totalflux = 0;

	nbins = 8;
	for ( long ii = 0; ii < subsize; ii ++ )
		for ( long jj = 0; jj < subsize; jj ++)
		{
			subdata[ii*subsize+jj] = 0;
			for ( int m = 0; m < nbins; m ++ )
				for ( int n = 0; n < nbins; n ++ )
				{
					xx = ii + (m + 0.5)/nbins;
					yy = jj + (n + 0.5)/nbins;
					rr = sqrt( (xx-subsize/2)*(xx-subsize/2) + (yy-subsize/2)*(yy-subsize/2) );
					subdata[ii*subsize+jj] += profile(rr) / (nbins * nbins);
				}

			totalflux += subdata[ii*subsize+jj];
		}

	totalMag = -2.5 * log10( totalflux ) + zeroMag;
}

double Galaxy::getMag()
{
	return totalMag;
}

double Galaxy::getAlpha()
{
	return alpha;
}

double Galaxy::profile( double rr )
{
	double Ir;
	Ir = I0 * exp(- pow((rr / alpha), 1 / nSersic));
	return Ir;
}

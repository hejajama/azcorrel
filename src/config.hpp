/*
 * Azimutal angle correlations
 *
 * Program-wide configurations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#ifndef AZCORREL_CONFIG
#define AZCORREL_CONFIG

//#define USE_MPI
//#define USE_FFTW
const double LAMBDAQCD = 0.241;
const int Nc=3;
const double Cf=(Nc*Nc-1.0)/(2.0*Nc);
const int Nf=3;

#endif


/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _CONFIG_HPP
#define _CONFIG_HPP 


#include <iostream>
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
#include <string>
using std::string;
#include <cmath>

typedef double REAL;
typedef unsigned int uint;

// Physical constants
//const double LAMBDAQCD2 = 0.21416*0.21416;   // GeV^2
const double LAMBDAQCD2 = 0.241*0.241;    // 0902.1112
const double LAMBDAQCD = 0.241;
//const double LAMBDAQCD2 = 0.2*0.2;    // 0704.0612, AN06 model
const int Nf=3;
const int Nc=3;
const double Cf = (Nc*Nc-1.0)/(2.0*Nc);
const double ALPHA_e = 1.0/137.035999679; 
const double e = sqrt(4.0*M_PI*ALPHA_e);

// Reqularization of the running coupling
// Berger&Stasto 1010.0671: 0.3
// Fit to HERA: 0.7 0902.1112
//const double MAXALPHA = 0.7;  
const double MAXALPHA = 0.5;  // 0704.0612
const double ALPHABAR_s = 0.05;       // \alphabar_s if RC=constant
const double ALPHAS = 0.2;    // compare with ALPHABAR_s.... ok, doesn't make sense

// Other constants

const double eps=0.000001;

// Inline functions

#define LINEINFO __FILE__ << ":" << __LINE__

inline const double SQR(const double x) { return x*x; }

enum Parton
{
    UVAL,DVAL,USEA,DSEA,U,D,S,C,B,G    // Valence quarks, sea quarks,all quarks, gluons
};

#endif

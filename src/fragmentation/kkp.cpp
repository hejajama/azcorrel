/*
 * KKP Fragmentation function
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "kkp.hpp"
#include <tools/config.hpp>
#include "fragmentation.hpp"

extern "C"
{
    // parton labels in output array:
    //  0    1    2    3    4    5    6    7    8     9    10
    //  g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
    
    void kkp_(int& h, int& set, REAL& x, REAL &qs, REAL output[]);
}

// D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
REAL KKP::Evaluate(Parton p, Hadron h, REAL x, REAL qs)
{
    int hadron=0; 
    REAL partons[11];
    partons[0]=0; partons[1]=1; partons[2]=2; partons[3]=3; partons[10]=10;
    if (h==PI) hadron=1;
    else if (h==K) hadron=2;
    else if (h==K0) hadron=3;
    else if (h==P) hadron=4;
    else if (h==PI0) hadron=5;
    else if (h==N) hadron=6;
    else if (h==H) hadron=7;

    int set=1;  // NLO
    kkp_(hadron, set, x, qs, partons);
    REAL result=0;
    if (p==G) result = partons[0];
    else if (p==U) result = partons[1];
    else if (p==D) result = partons[3];
    else if (p==S) result = partons[5];
    else if (p==C) result = partons[7];
    else if (p==B) result = partons[9];
    // Check: partons[2] is ubar etc...?
    else
        cerr << "Parton " << p << " is not supported! " << LINEINFO << endl;
    
    return result;

}



std::string KKP::GetString()
{
    return "KKP";
}

KKP::KKP()
{

}
 

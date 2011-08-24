#ifndef _FRAGMENTATION_H
#define _FRAGMENTATION_H

/*
 * Virtual class to hide different fragmentation functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <string>
#include "../config.hpp"

enum Hadron
{
    PI,   // pi+, pi-
    K,      // k+, k-
    K0,      // K^0, \bar K^0
    P,      // P, \bar P
    PI0,    // \pi^0
    N,      // N, \bar N
    H       // sum of pions, kaons and protons
};

class FragmentationFunction
{
    public:
        // D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
        virtual double Evaluate(Parton p, Hadron h, double x, double qs)=0;
        FragmentationFunction();
        virtual std::string GetString();
    private:


};


#endif

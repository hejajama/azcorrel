#ifndef _KKP_H
#define _KKP_H

/*
 * KKP fragmentation function
 * Uses fragmentation_kkp.f downloaded from
 * http://www.desy.de/~poetter/kkp.html
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <string>
#include "../config.hpp"
#include "fragmentation.hpp"

class KKP : public FragmentationFunction
{
    public:
        REAL Evaluate(Parton p, Hadron h, REAL x, REAL qs);
        KKP();
        std::string GetString();
    private:


};


#endif

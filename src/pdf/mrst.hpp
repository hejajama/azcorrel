#ifndef _MRST_H
#define _MRST_H

/*
 * Wrapper class to use mrst99 PDF
 * Real code is in files mrst99.{cpp,h} downloaded from
 * http://durpdg.dur.ac.uk/HEPDATA/PDF
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include "../pdf.hpp"
#include <tools/config.hpp>
#include "mrst99.h"
#include <string>

class MRST : public PDF
{
    public:
        ~MRST();
        REAL xq(REAL x, REAL q, Parton p);    // return x*q(x,q)
        void Initialize(int param=-1);
        std::string GetString();
    private:
        c_mrst* mrst;
        int set;
};

#endif

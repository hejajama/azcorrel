#ifndef _CTEQ_H
#define _CTEQ_H

/*
 * Wrapper calss for CTEQ-TEA parton distribution functions CT10
 * Real code is in file CT10Pdf.f downloaded from
 * http://hep.pa.msu.edu/cteq/public/index.html
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <tools/config.hpp>
#include "../pdf.hpp"
#include <string>

class CTEQ : public PDF
{
        REAL xq(REAL x, REAL q, Parton p);    // return x*q(x,q)
        void Initialize(int param=-1);
        std::string GetString();

};

// The following functions are implemented in CT10Pdf.f file
extern "C"
{
    void setct10_(int& iset_);
    REAL ct10pdf_(int& iparton, REAL& x, REAL& q);
}

#endif

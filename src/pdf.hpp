#ifndef _PDF_H
#define _PDF_H

/*
 * Virtual class to hide different parton distribution functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include "config.hpp"
#include <string>

class PDF
{
    public:
        ~PDF();
        virtual REAL xq(REAL x, REAL q, Parton p)=0;    // return x*q(x,q)
        virtual void Initialize(int param=-1);
        virtual std::string GetString();

        void PlotPdf(REAL Q);
};

#endif

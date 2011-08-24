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
        virtual double xq(double x, double q, Parton p)=0;    // return x*q(x,q)
        virtual void Initialize(int param=-1);
        virtual std::string GetString();

        void PlotPdf(double Q);
};

#endif

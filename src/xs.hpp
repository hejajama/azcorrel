#ifndef _XS_H
#define _XS_H

/*
 * Class to calculate cross sections related to two-body correlations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "config.hpp"
#include "pdf.hpp"
#include "amplitudelib/amplitudelib.hpp"

class CrossSection2
{
    public:
        CrossSection2(AmplitudeLib* N_, PDF* pdf_);
        REAL dSigma(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL theta, REAL sqrts);
        REAL Sigma(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL sqrts);

        REAL z(REAL pt1, REAL pt2, REAL y1, REAL y2);
        REAL xa(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL sqrts);
        REAL xh(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL sqrts);
        REAL Delta(REAL pt1, REAL pt2, REAL theta);

    private:
        AmplitudeLib* N;
        PDF* pdf;


        
};


#endif

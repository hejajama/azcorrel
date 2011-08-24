#ifndef _XS_H
#define _XS_H

/*
 * Class to calculate cross sections related to two-body correlations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <tools/config.hpp>
#include "pdf.hpp"
#include <amplitudelib/amplitudelib.hpp>
#include "fragmentation/fragmentation.hpp"

class CrossSection2
{
    public:
        CrossSection2(AmplitudeLib* N_, PDF* pdf_, FragmentationFunction* frag);
        double dSigma_lo(double pt1, double pt2, double y1, double y2, double theta, double sqrts
            , bool multiply_pdf=true);
        double dSigma(double pt1, double pt2, double y1, double y2, double theta, double sqrts,
            bool pdf=true);
        double CorrectionTerm(double pt1, double pt2, double ya, double phi);
        double Sigma(double pt1, double pt2, double y1, double y2, double sqrts);

        double NPair(double theta, double sqrts);

        double z(double pt1, double pt2, double y1, double y2);
        double xa(double pt1, double pt2, double y1, double y2, double sqrts);
        double xh(double pt1, double pt2, double y1, double y2, double sqrts);
        double Delta(double pt1, double pt2, double theta);

        // \int dr S(r) J_1(k*r)
        double G(double kt, double x);

        PDF* Pdf();
        FragmentationFunction* FragFun(); 

    private:
        AmplitudeLib* N;
        PDF* pdf;
        FragmentationFunction* fragfun;
        
        double gcacheval, gcachek;
        double fcacheval, fcachek;
        
};


#endif

#ifndef _XS_H
#define _XS_H

/*
 * Class to calculate cross sections related to two-body correlations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <tools/config.hpp>
#include <complex>
#include "pdf.hpp"
#include <amplitudelib/amplitudelib.hpp>
#include "fragmentation/fragmentation.hpp"

class CrossSection2
{
    public:
        CrossSection2(AmplitudeLib* N_, PDF* pdf_, FragmentationFunction* frag);
        ~CrossSection2();
        double dSigma_lo(double pt1, double pt2, double y1, double y2, double theta, double sqrts
            , bool multiply_pdf=true);
        double dSigma(double pt1, double pt2, double y1, double y2, double theta, double sqrts,
            bool pdf=true);
        double CorrectionTerm(double pt1, double pt2, double ya, double phi, double z);
        double CorrectionTerm_nomc(double pt1, double pt2, double ya, double phi);
        double CorrectionTerm_fft(double pt1, double pt2, double ya, double phi, double z);
        void CalculateCorrection_fft(double ya, double z);
        double Sigma(double pt1, double pt2, double y1, double y2, double sqrts);

        double NPair(double theta, double sqrts);

        double z(double pt1, double pt2, double y1, double y2);
        double xa(double pt1, double pt2, double y1, double y2, double sqrts);
        double xh(double pt1, double pt2, double y1, double y2, double sqrts);
        double Delta(double pt1, double pt2, double theta);

        // \int dr r*S(r)*K_1(mzr)*J1(kr)
        double G(double kt, double x, double z=0);
        // \int dr r*S(r)*K_0(mzr)*J0(kr)
        double H(double kt, double x, double z);

        PDF* Pdf();
        FragmentationFunction* FragFun(); 

    private:
        AmplitudeLib* N;
        PDF* pdf;
        FragmentationFunction* fragfun;
        
        double gcacheval, gcachek;
        double fcacheval, fcachek;

        //FFT
        std::complex<double> *data,*transform;
        int Nd;
        double delta;

        double V2int(double rx, double ry, double v1x, double v1y, double y, double z);
        
};
const double M_Q =1; //0.14;  // GeV

#endif

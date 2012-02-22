#ifndef _XS_H
#define _XS_H

/*
 * Class to calculate cross sections related to two-body correlations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <tools/config.hpp>
#include <complex>
#include <pdf/pdf.hpp>
#include <amplitudelib/amplitudelib.hpp>
#include <fragmentation/fragmentation.hpp>
#include <tools/interpolation2d.hpp>

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
        double CorrectionTerm_fft(double pt1, double pt2, double ya, double phi);
        void CalculateCorrection_fft(double ya, double z);
        double Sigma(double pt1, double pt2, double y1, double y2, double sqrts);


        double z(double pt1, double pt2, double y1, double y2);
        double xa(double pt1, double pt2, double y1, double y2, double sqrts);
        double xh(double pt1, double pt2, double y1, double y2, double sqrts);
        double Delta(double pt1, double pt2, double theta);

        
        double dSigma_full(double pt1, double pt2, double y1, double y2, double phi,
            double sqrts, bool deuteron=false);
        double dSigma_integrated(double minpt1, double minpt2, double miny, double maxy,
            double phi, double sqrts, bool deuteron=false);

        // \int dr r*S(r)*K_1(mzr)*J1(kr)
        double G(double kt, double x, double z=0);
        // \int dr r*S(r)*K_0(mzr)*J0(kr)
        double H(double kt, double x, double z);

        PDF* Pdf();
        FragmentationFunction* FragFun();

        AmplitudeLib* GetN();

        Interpolator2D *Ptinterpolator2d();
        Interpolator2D *Ptinterpolator2d_rev();

        void SetMCIntPoints(unsigned long long points);

        int LoadPtData(double y1, double y2);
        void Prepare2DInterpolators(double phi);    // initialize ptinerpolator2d{,_rev}

        void SetM_Q(double mq);
        double M_Q();

        double MaxPt(); // ptvals[ptvals.size()-1]

    private:
        double m_q; // Quark mass
        AmplitudeLib* N;
        PDF* pdf;
        FragmentationFunction* fragfun;
        
        double gcacheval, gcachek, gcachex, gcachez;
        double hcacheval, hcachek, hcachex, hcachez;

        //FFT
        std::complex<double> *data,*transform;
        int Nd;
        double delta;

        double V2int(double rx, double ry, double v1x, double v1y, double y, double z);
        unsigned long long mcintpoints;

        // [pt1ind][pt2ind]
        std::vector<std::vector< Interpolator*> > ptinterpolators;
        std::vector<std::vector< Interpolator*> > ptinterpolators_rev;
        Interpolator2D *ptinterpolator2d;
        Interpolator2D *ptinterpolator2d_rev;
        std::vector<std::vector< Interpolator*> > ptinterpolators_correction;
        std::vector<std::vector< Interpolator*> > ptinterpolators_rev_correction;
        Interpolator2D *ptinterpolator2d_correction;
        Interpolator2D *ptinterpolator2d_rev_correction;

        std::vector<double> ptvals;
        std::vector<std::string> ptstrings;

        bool apply_corrections; 
        
};

#endif

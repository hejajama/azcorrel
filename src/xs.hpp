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

enum Channel
{
	Q_QG,		// q->qg
	G_QQ,		// g->qq
	G_GG		// g->gg
};

class CrossSection2
{
    public:
        CrossSection2(AmplitudeLib* N_, PDF* pdf_, FragmentationFunction* frag);
        ~CrossSection2();
        double dSigma_lo(double pt1, double pt2, double y1, double y2, double theta, double sqrts
            , bool multiply_pdf=false);
        double dSigma(double pt1, double pt2, double y1, double y2, double theta, double sqrts);
        double CorrectionTerm(double pt1, double pt2, double ya, double phi, double z);
       
        double CorrectionTerm_nomc(double pt1, double pt2, double ya, double phi);
        double CorrectionTerm_fft(double pt1, double pt2, double ya, double phi);
        void CalculateCorrection_fft(double ya, double z);
        double Sigma(double pt1, double pt2, double y1, double y2, double sqrts);
        
        double GluonGluonGluon(double pt1, double pt2, double y1, double y2, double theta, double sqrts);
        double GluonCorrectionTerm(double pt1, double pt2, double ya, double phi, double z);

		double GluonQuarkQuark(double pt1, double pt2, double y1, double y2, double theta, double sqrts);


        double Z(double pt1, double pt2, double y1, double y2);
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
        
        double WavefSqr_k(double k1, double k2, double k1_dot_k2, double zfrac);	// \sum \psi(k1)\psi^*(k2)

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
        
        bool FiniteNc();
        void SetFiniteNc(bool fnc);
        
        void SetChannel(Channel c_);
        
        Interpolator2D *ptinterpolator2d;
        Interpolator2D *ptinterpolator2d_rev;
        
        void SetGluon(bool gluon_);

		double Alpha_s(double Qsqr, double scaling=1.0);
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
        unsigned long long mcintpoints, mcintpoints_orig;

        // [pt1ind][pt2ind]
        std::vector<std::vector< Interpolator*> > ptinterpolators;
        std::vector<std::vector< Interpolator*> > ptinterpolators_rev;
        
        
        std::vector<std::vector< Interpolator*> > ptinterpolators_correction;
        std::vector<std::vector< Interpolator*> > ptinterpolators_rev_correction;
        Interpolator2D *ptinterpolator2d_correction;
        Interpolator2D *ptinterpolator2d_rev_correction;
        
        // Where to read ptdata
		string fileprefix; 
		string fileprefix_cor; 
		string postfix;

        std::vector<double> ptvals;
        std::vector<std::string> ptstrings;
        
        bool finite_nc;

        bool apply_corrections; 
        
        
        Channel channel;
        
};


struct Inthelper_correction
{
    double pt1,pt2,ya;
    double phi;
    AmplitudeLib* N;
    size_t calln; int monte;
    CrossSection2* xs;
    bool finite_nc;
    bool gluon;

    double u1,u2,r,theta1,theta2,thetar,z;
};

double step(double x);

#endif

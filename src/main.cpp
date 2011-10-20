#include <pdf/pdf.hpp>
#include <pdf/mrst.hpp>
#include <pdf/cteq.hpp>
#include <tools/config.hpp>
#include <tools/tools.hpp>
#include <amplitudelib/amplitudelib.hpp>
#include "xs.hpp"
#include <fragmentation/kkp.hpp>
#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <sstream>


using std::cout; using std::endl;


void PlotPdf(double Q);

int main(int argc, char* argv[])
{
    if (string(argv[1])=="-help")
    {
        cout << "-y1 y: set rapidity of particle 1" << endl;
        cout << "-y2 2: set rapidity of particle 2" << endl;
        cout << "-pt1 pt: rapidity of particle 1" << endl;
        cout << "-pt2 pt: rapidit yof particle 2" << endl;
        cout << "-sqrts sqrts: \\sqrts(s)" << endl;
        cout << "-pdfmode [MRST, CTEQ]: choose pdf" << endl;
        cout << "-pdf Q: plot PDF as a function of x" << endl;
        cout << "-amplitude filename: where to load amplitude data" << endl;
        cout << "-mcintpoints points" << endl;
        return 0;
    }

    std::stringstream s;
    for (int i=0; i<argc; i++)
        s << argv[i] << " ";
    cout << "# " << s.str() << endl;

    gsl_set_error_handler(&ErrHandler);
    std::string filename="amplitude.dat";
    

    double y1=3.5; double y2=2; double pt1=5; double pt2=3; double sqrts=200;
    double Q=-1;
    PDF *pdf=0;
    size_t mcintpoints=1e7;
    
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-y1")
            y1 = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-y2")
            y2= StrToReal(argv[i+1]);
        else if (string(argv[i])=="-pt1")
            pt1 = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-pt2")
            pt2 = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-sqrts")
            sqrts = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-amplitude")
            filename = argv[i+1];
        else if (string(argv[i])=="-pdf")
            Q=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-pdfmode")
        {
            if (string(argv[i+1])=="MRST")
            {
                pdf = new MRST();
                pdf->Initialize();
            } else if (string(argv[i+1])=="CTEQ")
            {
                pdf = new CTEQ();
                pdf->Initialize();
            }
            else
            {
                cerr << "Unrecognized PDF " << argv[i+1] << endl;
                return -1;
            }
        }
        else if (string(argv[i])=="-mcintpoints")
        {
            std::stringstream tmpstr;
            tmpstr << argv[i+1];
            tmpstr >> mcintpoints;
            
        }
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }

    if (!pdf)   // Not initializd, use default
    {
        pdf = new CTEQ();
        pdf->Initialize();
    }
   
    /*if (pt1*std::exp(y1) < pt2*std::exp(y2))
    {
        cerr <<"pt1*exp(y1) < pt2*exp(y2) " << endl;
        return 0;
    }*/
    
    cout <<"# Parton distribution function used: " << pdf->GetString() << endl;

    cout <<"# Quark mass: " << M_Q << " GeV" << endl;

    
    if (Q>=0)   // Plot PDF and exit
    {
        pdf->PlotPdf(Q);
        delete pdf;
        return 0;
    }
    AmplitudeLib amplitude(filename);
    KKP fragmentation;
/*
    for (double x=1e-4; x<1; x*=1.1)
    {
        cout << x << " " << fragmentation.Evaluate(G, H, x, 80.2) << endl;
    }
    delete pdf; return 0;
  */  
    CrossSection2 cross_section(&amplitude, pdf, &fragmentation);
    cross_section.SetMCIntPoints(mcintpoints);
    double ya=std::log(0.01 / cross_section.xa(pt1,pt2,y1,y2,sqrts));
    cout << "# pt1=" << pt1 <<", pt2=" << pt2 <<", y1=" << y1 <<", y2=" << y2 <<
    " y_A=" << ya << endl;
    cout << "# x1=" << pt1*std::exp(y1)/sqrts << ", x2=" << pt2*std::exp(y2)/sqrts
        << " sqrts=" << sqrts << endl;
    cout << "# z=" << cross_section.z(pt1,pt2,y1,y2) <<", 1-z=" << cross_section.z(pt2,pt1,y2,y1)
    << " xa=" << cross_section.xa(pt1,pt2,y1,y2,sqrts)
    << " xh=" << cross_section.xh(pt1,pt2,y1,y2,sqrts) << endl;
    cout << "# Q_s = " << 1.0/amplitude.SaturationScale(
        std::log(0.01/cross_section.xa(pt1,pt2,y1,y2,sqrts)), 0.5 ) << " GeV " << endl;
    cout << "# MC Integration points " << mcintpoints << endl;


    amplitude.InitializeInterpolation(
        std::log(0.01 / cross_section.xa(pt1,pt2,y1,y2,sqrts)) );
    double normalization = 1;//cross_section.Sigma(pt1, pt2, y1, y2, sqrts);
    //cout << "# Normalization totxs " << normalization << endl;
    //cout << "# Theta=2.5 " << cross_section.dSigma(pt1,pt2,y1,y2,2.5,sqrts) << endl;
    int points=50;
    bool fftw=false;

    // FFTW
    if (fftw)
        cross_section.CalculateCorrection_fft(ya, cross_section.z(pt1,pt2,y1,y2));
        
    //#pragma omp parallel for
    for (int i=0; i<points; i++)
    {
        //double theta = 0.2*2.0*M_PI*(i+1.0);
        double theta = M_PI/10.0 + (2.0*M_PI-M_PI/10.0)*i/((double)points-1.0);
        double result = cross_section.dSigma(pt1,pt2,y1,y2,theta,sqrts,!fftw);
        if (result<-0.5)
        {
            #pragma omp critical
            cout <<"# " << theta <<" MC integral failed " << endl;
            continue;
        }
        
        #pragma omp critical
        {
            cout << theta << " " << result/normalization << " "
                //<< cross_section.CorrectionTerm_fft(pt1,pt2, ya, theta)
                << endl;
        }
        
    }


    delete pdf;
    



    return 0;
}





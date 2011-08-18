#include "pdf.hpp"
#include "pdf/mrst.hpp"
#include "pdf/cteq.hpp"
#include "config.hpp"
#include "tools.hpp"
#include "amplitudelib/amplitudelib.hpp"
#include "xs.hpp"
#include "fragmentation/kkp.hpp"
#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>


using std::cout; using std::endl;


void PlotPdf(REAL Q);

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
        return 0;
    }

    gsl_set_error_handler(&ErrHandler);
    std::string filename="amplitude.dat";
    

    REAL y1=3.5; REAL y2=2; REAL pt1=5; REAL pt2=3; REAL sqrts=200;
    REAL Q=-1;
    PDF *pdf=0;
    
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
   
    if (pt1*std::exp(y1) < pt2*std::exp(y2))
    {
        cerr <<"pt1*exp(y1) < pt2*exp(y2) " << endl;
        return 0;
    }
    
    cout <<"# Parton distribution function used: " << pdf->GetString() << endl;

    
    if (Q>=0)   // Plot PDF and exit
    {
        pdf->PlotPdf(Q);
        delete pdf;
        return 0;
    }
    AmplitudeLib amplitude(filename);
    KKP fragmentation;
/*
    for (REAL x=1e-4; x<1; x*=1.1)
    {
        cout << x << " " << fragmentation.Evaluate(G, H, x, 80.2) << endl;
    }
    delete pdf; return 0;
  */  
    CrossSection2 cross_section(&amplitude, pdf, &fragmentation);

    cout << "# pt1=" << pt1 <<", pt2=" << pt2 <<", y1=" << y1 <<", y2=" << y2 <<
    " y_A=" << std::log(0.01 / cross_section.xa(pt1,pt2,y1,y2,sqrts)) << endl;
    cout << "# z=" << cross_section.z(pt1,pt2,y1,y2) <<", 1-z=" << cross_section.z(pt2,pt1,y2,y1)
    << " xa=" << cross_section.xa(pt1,pt2,y1,y2,sqrts)
    << " xh=" << cross_section.xh(pt1,pt2,y1,y2,sqrts) << endl;
    cout << "# Q_s = " << 1.0/amplitude.SaturationScale(
        std::log(0.01/cross_section.xa(pt1,pt2,y1,y2,sqrts)), 0.5 ) << " GeV " << endl;

    

    ///DEBUG    
    /*amplitude.InitializeInterpolation(y1);
    
    
    for (int kind=0; kind<80; kind++)
    {
        REAL k = 0.3*std::pow(1.1,kind);
        REAL result = 1.0-k*cross_section.G(k, 0.01*std::exp(-y1));
        //REAL result = SQR(k)*amplitude.S_k(k, y1);
        //REAL result  =amplitude.N_k(k,y);
        cout << SQR(k) << " " << result << endl;

    }
    delete pdf; return 0;
    */

    amplitude.InitializeInterpolation(
        std::log(0.01 / cross_section.xa(pt1,pt2,y1,y2,sqrts)) );
    REAL normalization =1.0;//= cross_section.Sigma(pt1, pt2, y1, y2, sqrts);
    cout << "# Normalization totxs " << normalization << endl;

    for (REAL theta=0.15; theta<2.0*M_PI-0.15; theta+=0.15)
    {
        //REAL result = cross_section.dSigma(pt1,pt2,y1,y2,theta,sqrts)
        //    + cross_section.dSigma(pt2,pt1,y2,y1,theta,sqrts);
        REAL result = cross_section.NPair(theta, sqrts);
        cout << theta << " " << result/normalization << endl;
    }


    delete pdf;
    



    return 0;
}





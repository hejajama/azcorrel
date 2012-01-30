#include <pdf/pdf.hpp>
#include <pdf/mrst.hpp>
#include <pdf/cteq.hpp>
#include <tools/config.hpp>
#include <tools/tools.hpp>
#include <amplitudelib/amplitudelib.hpp>
#include <fragmentation/kkp.hpp>
#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <sstream>
#include <fstream>
#include <limits>

#include "config.hpp"
#ifdef USE_MPI
    #include <mpi.h>
#endif

#include "xs.hpp"



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
        cout << "-phi phi: calculate only one angle" << endl;
        cout << "-pdfmode [MRST, CTEQ]: choose pdf" << endl;
        cout << "-pdf Q: plot PDF as a function of x" << endl;
        cout << "-nopdf: don't multiply by pdf, prints ~d\\sigma/d^3kd^3q" << endl;
        cout << "-amplitude filename: where to load amplitude data" << endl;
        cout << "-mcintpoints points" << endl;
        cout << "-mq quark_mass in GeV" << endl;
        cout << "-deuteron: use deuteron as a probe" << endl;
        cout << "-output_file filename" << endl;
        return 0;
    }

    #ifdef USE_MPI
    int rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        cerr << "MPI Initialization failed" << endl;
        return -1;
    }
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    #endif

    std::stringstream s;
    for (int i=0; i<argc; i++)
        s << argv[i] << " ";
    #ifdef USE_MPI
    if (id==0)
    #endif
    cout << "# " << s.str() << endl;


    gsl_set_error_handler(&ErrHandler);
    std::string amplitude_filename="amplitude.dat";

    std::string output_file="";
    

    double y1=3.5; double y2=2; double pt1=5; double pt2=3; double sqrts=200;
    double Q=-1;
    double phi=-1;
    double mq=0.14;
    bool multiply_pdf=true;
    PDF *pdf=0;
    unsigned long long mcintpoints=1e7;
    
    bool deuteron=false;   

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
        else if (string(argv[i])=="-phi")
            phi = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-amplitude")
            amplitude_filename = argv[i+1];
        else if (string(argv[i])=="-pdf")
            Q=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-nopdf")
            multiply_pdf=false;
        else if (string(argv[i])=="-mq")
            mq = StrToReal(argv[i+1]);
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
        else if (string(argv[i])=="-deuteron")
            deuteron=true;
        else if (string(argv[i])=="-output_file")
            output_file = argv[i+1];
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
    
   
    if (Q>=0)   // Plot PDF and exit
    {
        pdf->PlotPdf(Q);
        delete pdf;
        return 0;
    }
    AmplitudeLib amplitude(amplitude_filename);
    KKP fragmentation;

    CrossSection2 cross_section(&amplitude, pdf, &fragmentation);
    cross_section.SetMCIntPoints(mcintpoints);
    cross_section.SetM_Q(mq);
    double ya=std::log(amplitude.X0() / cross_section.xa(pt1,pt2,y1,y2,sqrts));

    std::ofstream output;
    #ifdef USE_MPI
    if (id==0)
    {
    #endif
    if (output_file!="")
        output.open(output_file.c_str());

    std::stringstream infostr;

    if (multiply_pdf)
        infostr << "# Parton distribution function used: " << pdf->GetString() << endl;
    else
        infostr <<"# Not multiplying by PDF" << endl;
    if (deuteron) infostr << "# Probe is deuteron" << endl;
    else infostr <<"# Probe is proton" << endl;
    infostr <<"# Quark mass: " << cross_section.M_Q() << " GeV" << endl;

    infostr << "# pt1=" << pt1 <<", pt2=" << pt2 <<", y1=" << y1 <<", y2=" << y2 <<
    " y_A=" << ya << endl;
    infostr << "# x1=" << pt1*std::exp(y1)/sqrts << ", x2=" << pt2*std::exp(y2)/sqrts
        << " sqrts=" << sqrts << endl;
    infostr << "# z=" << cross_section.z(pt1,pt2,y1,y2) <<", 1-z=" << cross_section.z(pt2,pt1,y2,y1)
    << " xa=" << cross_section.xa(pt1,pt2,y1,y2,sqrts)
    << " xh=" << cross_section.xh(pt1,pt2,y1,y2,sqrts) << endl;
    infostr << "# Q_s = " << 1.0/amplitude.SaturationScale(ya, 0.22) << " GeV " << endl;
    infostr << "# MC Integration points " << mcintpoints << " / supported: "
        << std::numeric_limits<unsigned long long>::max() << endl;
    infostr << "# Max pt when loading data: " << cross_section.MaxPt() << endl;

    cout << infostr.str();
    if (output_file!="")
        output << infostr.str();

    #ifdef USE_MPI
    }
    #endif

    amplitude.InitializeInterpolation(ya);
    double normalization = 1;//cross_section.Sigma(pt1, pt2, y1, y2, sqrts);
    //cout << "# Normalization totxs " << normalization << endl;
    //cout << "# Theta=2.5 " << cross_section.dSigma(pt1,pt2,y1,y2,2.5,sqrts) << endl;
    //int points=10;
    int points=40;
    if (phi>-0.5) points=1;    // calculate only given angle
    //double maxphi=2.0*M_PI-0.5;
    double maxphi=M_PI;
    double minphi = 0.5;
    bool fftw=false;

    // FFTW
    if (fftw)
        cross_section.CalculateCorrection_fft(ya, cross_section.z(pt1,pt2,y1,y2));

    if (fftw and multiply_pdf)
        cerr <<"Can't calculate FFT and multiply by PDF!" << endl;

    /*if (!fftw and false)
        cross_section.LoadPtData(y1,y2);
    */
    int ready=0;
    #pragma omp parallel for
    for (int i=0; i<points; i++)
    {
        //double theta = 0.2*2.0*M_PI*(i+1.0);
        double theta = minphi + (maxphi-minphi)*i/((double)points-1.0);
        if (phi>-0.5) theta=phi;    // calculate given value
        double result=0;
        
        #pragma omp critical
        {
            ready++;
            #ifdef USE_MPI
            if (id==0)
            #endif
            cout << "# Starting index " << ready << "/" << points << " angle " << theta << endl;
            
        }
        
        if (!fftw)
            result = cross_section.dSigma_integrated(2, 1, 2.4, 4, theta, sqrts, deuteron);
                /*+ cross_section.dSigma_integrated(2,1,3.2,2.4, theta, sqrts, deuteron); 
            */
            //result = cross_section.dSigma_full(pt1,pt2,y1,y2,theta,sqrts, deuteron);
            //result = cross_section.dSigma(pt1,pt2,y1,y2,theta,sqrts,multiply_pdf);
            //    + cross_section.dSigma(pt2,pt1,y2,y1,theta,sqrts,multiply_pdf);
        else
            result = cross_section.CorrectionTerm_fft(pt1, pt2, ya, theta);
        /*if (result<-0.5)
        {
            //cout <<"# " << theta <<" MC integral failed " << endl;
            continue;
        }*/
        
        #pragma omp critical
        {
            #ifdef USE_MPI
            if (id==0)
            {
            #endif
            std::stringstream tmp_out;
            tmp_out << theta << " " << result/normalization << " "
                //<< cross_section.dSigma_lo(pt1, pt2, y1, y2, theta, sqrts, multiply_pdf) << " "
                //<< cross_section.CorrectionTerm_fft(pt1,pt2, ya, theta)
                << endl;
            cout << tmp_out.str();
            if (output_file!="")
                output << tmp_out.str();
            #ifdef USE_MPI
            }
            #endif
        }
        
    }


    delete pdf;

    #ifdef USE_MPI
    //if (id==0)
        MPI_Finalize();
    if (id==0)
    {
    #endif
        if (output_file!="")
            output.close();
    #ifdef USE_MPI
    }
    #endif

    return 0;
}





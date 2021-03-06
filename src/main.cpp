#include <pdf/pdf.hpp>
#include <pdf/mrst.hpp>
#include <pdf/cteq.hpp>
#include <tools/config.hpp>
#include <tools/tools.hpp>
#include <amplitudelib/amplitudelib.hpp>
#include <fragmentation/kkp.hpp>
#include <fragmentation/dss.hpp>
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

enum MODE
{
    MODE_CP,
    MODE_FIXED_PT_Y,
    MODE_DSIGMA
};

void PlotPdf(double Q);

int main(int argc, char* argv[])
{
    if (string(argv[1])=="-help")
    {
        cout << "-mode [CP,FIXED,DSIGMA]" << endl;
        cout << "-y1 y: set rapidity of particle 1" << endl;
        cout << "-y2 2: set rapidity of particle 2" << endl;
        cout << "-pt1 pt: rapidity of particle 1" << endl;
        cout << "-pt2 pt: rapidit yof particle 2" << endl;
        cout << "-sqrts sqrts: \\sqrts(s)" << endl;
        cout << "-phi phi: calculate only one angle" << endl;
        cout << "-minphi phi, -maxphi phi" << endl;
        cout << "-pdfmode [MRST, CTEQ]: choose pdf" << endl;
        cout << "-nopdf: don't multiply by pdf, prints ~d\\sigma/d^3kd^3q" << endl;
        cout << "-amplitude filename: where to load amplitude data" << endl;
        cout << "-mcintpoints points" << endl;
        cout << "-mq quark_mass in GeV" << endl;
        cout << "-deuteron: use deuteron as a probe" << endl;
        cout << "-output_file filename" << endl;
        cout << "-points points" << endl;
        cout << "-finite_nc: use finite Nc when calculating quadrupole" << endl;
        cout << "-channel Q_QG [default], G_QQ, G_GG" << endl;
        return 0;
    }

    MODE mode = MODE_DSIGMA;
    int points=10;
    double minphi=0.5, maxphi=M_PI;
    bool fftw=false;
    Channel ch = Q_QG;
    
    bool finite_nc=false;

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
    unsigned long long mcintpoints=1e8;
    
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
        else if (string(argv[i])=="-minphi")
            minphi = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxphi")
            maxphi = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-amplitude" or string(argv[i])=="-data")
            amplitude_filename = argv[i+1];
        else if (string(argv[i])=="-nopdf")
            multiply_pdf=false;
        else if (string(argv[i])=="-mq")
            mq = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-fftw")
			fftw=true;
        else if (string(argv[i])=="-mode")
        {
            if (string(argv[i+1])=="CP")
                mode=MODE_CP;
            else if (string(argv[i+1])=="DSIGMA")
                mode = MODE_DSIGMA;
            else if (string(argv[i+1])=="FIXED")
                mode = MODE_FIXED_PT_Y;
            else
            {
                cerr << "Unknown mode " << argv[i+1] << endl;
                return -1;
            }
            
        }
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
        else if (string(argv[i])=="-points")
            points = StrToInt(argv[i+1]);
        else if (string(argv[i])=="-deuteron")
            deuteron=true;
        else if (string(argv[i])=="-output_file")
            output_file = argv[i+1];
        else if (string(argv[i])=="-finite_nc")
			finite_nc=true;
		else if (string(argv[i])=="-channel")
		{
			if (string(argv[i+1])=="Q_QG")
				ch=Q_QG;
			else if (string(argv[i+1])=="G_QQ")
				ch=G_QQ;
			else if (string(argv[i+1])=="G_GG")
				ch=G_GG;
			else
			{
				cerr << "Unknown channel " << argv[i+1] << endl;
				return -1;
			}
		}
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }

    if (!pdf and mode != MODE_DSIGMA)   // Not initializd, use default
    {
        pdf = new CTEQ();
		pdf->SetOrder(Amplitude::LO);
        pdf->Initialize();
    }
    

    AmplitudeLib amplitude(amplitude_filename);
    DSS *fragmentation = NULL;
    
    if (mode != MODE_DSIGMA)
	{
		fragmentation = new DSS();
    	fragmentation -> SetOrder(Amplitude::LO);
	}
    CrossSection2 cross_section(&amplitude, pdf, fragmentation);
    cross_section.SetMCIntPoints(mcintpoints);
    cross_section.SetM_Q(mq);
    cross_section.SetChannel(ch);
    cross_section.SetFiniteNc(finite_nc);
    double ya=std::log(amplitude.X0() / cross_section.xa(pt1,pt2,y1,y2,sqrts));
    
    /*
    amplitude.InitializeInterpolation(y1);
    for (double pt=1e-4; pt<400; pt*=1.1)
    {
        cout << pt << " " << amplitude.S_k(pt, y1) << " " << pt*pt*amplitude.S_k(pt, y1) << endl;
        //cout << pt << " " << 1.0 - pt*cross_section.G(pt, 0.02*std::exp(-y1), 0.5) << endl;
    }
    return 0;
    */
    
    std::ofstream output;
    #ifdef USE_MPI
    if (id==0)
    {
    #endif
    if (output_file!="")
        output.open(output_file.c_str());

    std::stringstream infostr;
    infostr <<"# Using AmplitudeLib v. " << amplitude.Version() << endl;

    if (multiply_pdf)
    {
        infostr << "# Parton distribution function used: " << pdf->GetString() << endl;
        infostr << "# Fragfun: " << fragmentation->GetString() << endl;
    }
    else
        infostr <<"# Not multiplying by PDF" << endl;
    if (deuteron) infostr << "# Probe is deuteron" << endl;
    else infostr <<"# Probe is proton" << endl;
    infostr <<"# Quark mass: " << cross_section.M_Q() << " GeV" << endl;
    infostr << "# pt1=" << pt1 <<", pt2=" << pt2 <<", y1=" << y1 <<", y2=" << y2 <<
    " y_A=" << ya << endl;
    infostr << "# x1=" << pt1*std::exp(y1)/sqrts << ", x2=" << pt2*std::exp(y2)/sqrts
        << " sqrts=" << sqrts << endl;
    infostr << "# z=" << cross_section.Z(pt1,pt2,y1,y2) <<", 1-z=" << 1.0-cross_section.Z(pt1,pt2,y1,y2)
        << " xa=" << cross_section.xa(pt1,pt2,y1,y2,sqrts)
        << " xh=" << cross_section.xh(pt1,pt2,y1,y2,sqrts)  
        << " x_0=" << amplitude.X0() << endl;
    infostr << "# Q_s = " << 1.0/amplitude.SaturationScale(amplitude.X0()*std::exp(-ya), 0.22) << " GeV " << endl;
    infostr << "# MC Integration points " << mcintpoints << " / supported: "
        << std::numeric_limits<unsigned long long>::max() << endl;
    infostr << "# Max pt when loading data: " << cross_section.MaxPt() << endl;
    infostr << "# Channel: "; if (ch==Q_QG) infostr << " Q->QG"; else if (ch==G_QQ) infostr << "G->QQ"; else if (ch==G_GG) infostr << "G->GG"; infostr << endl;
    
    if (cross_section.xh(pt1,pt2,y1,y2,sqrts)>1  and mode==MODE_DSIGMA)
    {
		infostr << "# Kinematically forbidden parton level process!" << endl;
	}

    cout << infostr.str();
    if (output_file!="")
        output << infostr.str();

    #ifdef USE_MPI
    }
    #endif    
    
    if (ya<0)
    {
        cerr << "Negative rapidity " << ya << " at " << LINEINFO << endl;
       // return -1;
    }   
    
    amplitude.InitializeInterpolation(amplitude.X0()*std::exp(-ya));
    double normalization = 1;//cross_section.Sigma(pt1, pt2, y1, y2, sqrts);
    if (phi>-0.5) points=1;    // calculate only given angle

    // FFTW
    if (fftw)
        cross_section.CalculateCorrection_fft(ya, cross_section.Z(pt1,pt2,y1,y2));

    if (fftw and multiply_pdf)
        cerr <<"Can't calculate FFT and multiply by PDF!" << endl;

    if (mode == MODE_FIXED_PT_Y)
        cross_section.LoadPtData(y1,y2);
    int ready=0;
    
    //#pragma omp parallel for
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

        if (mode == MODE_FIXED_PT_Y)
        {   
            cross_section.Prepare2DInterpolators(theta);
            result = cross_section.dSigma_full(pt1,pt2,y1,y2,theta,sqrts, deuteron);
        }
        if (mode == MODE_CP)
            result = cross_section.dSigma_integrated(2, 1, 2.4, 4, theta, sqrts, deuteron);
        
        if (mode==MODE_DSIGMA)
        {
            result = cross_section.dSigma(pt1,pt2,y1,y2,theta,sqrts);
            if (result < -999999)	// convergence problem
				return -1;
        }

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

	if (pdf != NULL)
		delete pdf;
	if (fragmentation != NULL)
		delete fragmentation;

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





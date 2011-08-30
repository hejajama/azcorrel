/*
 * Calculate correction term from 4-point function using fft
 * Idea: we integrate over r, v and v' (v1, v2 in the code). Integrations over
 * r and v are 2D FFT's, so perform them by using FFTW. The last integral over
 * v2 depends thus on 4 d.o.f. Basically, we have one 4-dimensional FT, so
 * we first fill a 4-dimensional array and then FT it.
 */

#include <complex>
#include <fftw3.h>
#include <gsl/gsl_integration.h>
#include "xs.hpp"
#include <amplitudelib/amplitudelib.hpp>
#include <tools/tools.hpp>
#include <fstream>
#include <string>
#include <sstream>

const double R_RANGE = 50.0;

const bool READFILE = false; // Read integrand from a file

double CrossSection2::CorrectionTerm_fft(double pt1, double pt2,
    double ya, double phi)
{
    const int Nd = 4;   // Number of datapoints for each dimension
                            // must be 2^N for some N
    const double maxr = R_RANGE;   // Maximum length of vector component
    const double delta = maxr/Nd;

    std::complex<double> *data = new std::complex<double>[Nd*Nd*Nd*Nd];
    N->SetOutOfRangeErrors(false);
    for (int i=0; i<Nd*Nd*Nd*Nd; i++)
    {
        if (READFILE) data[i]=0;
        else data[i]=-9999;
    }

        // Fill array
        // Coordinates: 0=r_x, 1=r_y, 2=v1_x, 3=v1_y

    if (READFILE)
    {
        // Read from file
        std::ifstream file("integrand.txt");
        std::string line;
        while (getline(file, line))
        {
            string id,val;
            std::istringstream liness(line);
            getline( liness, id, ' ');
            getline( liness, val, ' ');
            int i=StrToInt(id);; double v=StrToReal(val);
            data[i] = v;
        }
        file.close();
    }
    else
    {
        // Calculate
        int done=0;
        #pragma omp parallel for
        for (int rxind=0; rxind<Nd; rxind++)
        {
            double rx = rxind*delta - maxr/2.0;
            int rxi = rxind-Nd/2;    // table index
            if (rx < 0) rxi = rxind + Nd/2;
            
            for (int ryind=0; ryind<Nd; ryind++)
            {
                double ry = ryind*delta - maxr/2.0;
                double ryi = ryind-Nd/2;
                if (ry<0) ryi = ryind + Nd/2;
                
                for (int vxind=0; vxind<Nd; vxind++)
                {
                    double vx = vxind*delta - maxr/2.0;
                    double vxi = vxind - Nd/2;
                    if (vx<0) vxi = vxind + Nd/2;
                    
                    for (int vyind=0; vyind<Nd; vyind++)
                    {
                        double vy = vyind*delta - maxr/2.0;
                        double vyi = vyind-Nd/2;
                        if (vy<0) vyi = vyind + Nd/2;

                        int index = vyi + Nd*(vxi + Nd*(ryi + Nd*rxi) );                        
                        
                        // Mirror: x components -> -x
                        
                        int rxi_xm = rxi;
                        int vxi_xm = vxi;
                        if (rxi != 0) rxi_xm = Nd-rxi;
                        if (vxi != 0) vxi_xm = Nd-vxi;
                        int index2 = vyi + Nd*(vxi_xm + Nd*(ryi + Nd*rxi_xm) );

                        // Mirror y
                        int ryi_ym = ryi;
                        int vyi_ym = vyi;
                        if (ryi != 0) ryi_ym = Nd-ryi;
                        if (vyi != 0) vyi_ym = Nd-vyi;

                        int index3 = vyi_ym + Nd*(vxi + Nd*(ryi_ym + Nd*rxi) );

                        // Mirror both
                        int index4 = vyi_ym + Nd*(vxi_xm + Nd*(ryi_ym + Nd*rxi_xm) );
                        
                        if (data[index].real()>-9990)
                        {
                            if (data[index2].real()<-9000 or data[index3].real()<-9000
                                or data[index4].real()<-9000)
                                cerr << "Something wrong with mirrors...\n";
                            continue;
                        }

                        double result = V2int(rx, ry, vx, vy, ya);

                        
                        
                        #pragma omp critical
                        {
                            data[index]=result;
                            data[index2]=result;
                            data[index3]=result;
                            data[index4]=result;
                            done++;
                            if (data[index].real() != 0 or data[index].imag()!=0)
                            {
                                cerr << " # done " << done << "/" << std::pow(Nd, 4)/4 << " prevres " << data[index]
                                << " rx " << rx << " ry " << ry << " vx " << vx << " vy " << vy << endl;
                                cout << index << " " << data[index].real() <<
                                " " << rxi << " " << ryi << " " << vxi << " " << vyi << endl;
                                cout << index2 << " " << data[index].real() << " " << rxi_xm << " "
                                 << ryi << " " << vxi_xm << " " << vyi << endl;
                                cout << index3 << " " << data[index].real() <<
                                " " << rxi << " " << ryi_ym << " " << vxi << " " << vyi_ym << endl;
                                cout << index4 << " " << data[index].real() <<
                                " " << rxi_xm << " " << ryi_ym << " " << vxi_xm << " " << vyi_ym << endl;
                            }
                        }
                    }
                }
            } // end for ryind
        } // end for rxind
    }   // end if(READFILE)
    







    delete[] data;

    return 0;
}

struct Inthelper_v2int
{
    double rx,ry,v1x,v1y;
    double r;
    double v2x;
    double y;
    AmplitudeLib* N;
};
const int VINTINTERVALS = 100;
const double VINTACCURACY = 0.1;
double Inthelperf_v2xint(double v2x, void* p);
double Inthelperf_v2yint(double v2y, void* p);

// Integrate over v2
double CrossSection2::V2int(double rx, double ry, double v1x, double v1y, double y)
{
    Inthelper_v2int helper;
    helper.rx=rx; helper.ry=ry; helper.v1x=v1x; helper.v1y=v1y; helper.y=y;
    helper.N=N; helper.r = std::sqrt(SQR(rx) + SQR(ry) );

    gsl_function f;
    f.function = Inthelperf_v2xint;
    f.params = &helper;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(VINTINTERVALS);
    double result,abserr;

    int status = gsl_integration_qag(&f, -R_RANGE/2.0, R_RANGE/2.0,
            0, VINTACCURACY, VINTINTERVALS, GSL_INTEG_GAUSS41, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "v2xint failed at " << LINEINFO <<", result " << result
        << " relerror " << std::abs(abserr/result) << " rx: " << rx << " ry "
        << ry << " v1x " << v1x << " v1y " << v1y << endl;

    return result;
}

double Inthelperf_v2xint(double v2x, void* p)
{
    Inthelper_v2int* par = (Inthelper_v2int*)p;
    par->v2x = v2x;

    gsl_function f;
    f.function = Inthelperf_v2yint;
    f.params = par;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(VINTINTERVALS);
    double result,abserr;

    int status = gsl_integration_qag(&f, -R_RANGE/2.0, R_RANGE/2.0,
            0, VINTACCURACY, VINTINTERVALS, GSL_INTEG_GAUSS41, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status and std::abs(result)>1e-10)
        cerr << "v2yint failed at " << LINEINFO <<", result " << result
        << " relerror " << std::abs(abserr/result) << " rx: " << par->rx << " ry "
        << par->ry << " v1x " << par->v1x << " v1y " << par->v1y << " "
        << " v2x " << v2x << endl;

    return result;
}

double Inthelperf_v2yint(double v2y, void* p)
{
    Inthelper_v2int* par = (Inthelper_v2int*)p;
    double rx = par->rx; double ry=par->ry;
    double v1x = par->v1x; double v1y = par->v1y;
    double v2x = par->v2x; double r = par->r;
    AmplitudeLib* N = par->N;
    double y = par->y;

    // vector lengths
    // r - v2 + 1/2v1
    double r_m_v2_p_05v1 = std::sqrt(
        SQR(rx - v2x + v1x/2.0) + SQR(ry - v2y + v1y/2.0) );
    if (r_m_v2_p_05v1 < N->MinR()) return 0;

    // r + 1/2v + v2
    double r_p_05v1_p_v2 = std::sqrt(
        SQR(rx + v1x/2.0 + v2x) + SQR(ry + v1y/2.0 + v2y) );
    if (r_p_05v1_p_v2 < N->MinR()) return 0;

    // r+v1
    double r_p_v1 = std::sqrt( SQR(rx+v1x) + SQR(ry+v1y) );
    if (r_p_v1 < N->MinR()) return 0;

    // 1/2v1 + v2
    double v2_p_05v1 = std::sqrt( SQR(v1x/2.0 + v2x) + SQR(v1y/2.0+v2y) );
    if (v2_p_05v1 < N->MinR()) return 0;

    // v2-1/2*v1
    double v2_m_05v1 = std::sqrt( SQR(v1x/2.0 - v2x) + SQR(v1y/2.0 - v2y) );
    if (v2_m_05v1 < N->MinR()) return 0;

    double minsprod = 1e-15;    /// TODO: result shouldn't depend on this
    
    REAL nom1 = N->S(r_m_v2_p_05v1, y)*N->S(r_p_05v1_p_v2, y);
    REAL denom1 = N->S(r,y)*N->S(r_p_v1 ,y);
    if (denom1 < minsprod) return 0;

    REAL nom2 = N->S(v2_p_05v1, y)*N->S(v2_m_05v1, y);
    if (nom1 < minsprod or nom2<minsprod or std::log(nom2/denom1)==0) return 0;
    REAL denom2 = denom1; // optimize, same as before

    double result;

    // F/F
    result = std::log(nom1/denom1) / std::log(nom2/denom2);

    result *= N->S(r_p_v1, y) *( N->S(v2_p_05v1, y) * N->S(v2_m_05v1, y)
        - N->S(r, y)*N->S(r_p_v1, y) );

    // wave function product
    result *= -16.0*SQR(M_PI)
        * ( SQR(v2x) + SQR(v2y) - 0.25*(SQR(v1x) + SQR(v1y)) )
        / ( SQR(v2_m_05v1) * SQR(v2_p_05v1) );

    if (isnan(result))
        cerr << "NaN at " << LINEINFO << endl;
    
    return result;

}

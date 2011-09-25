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
#include <gsl/gsl_sf_bessel.h>
#include "xs.hpp"
#include <amplitudelib/amplitudelib.hpp>
#include <tools/tools.hpp>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

const double R_RANGE = 8;

const bool READFILE = false; // Read integrand from a file

void CrossSection2::CalculateCorrection_fft(double ya, double z)
{
    Nd = 36;   // Number of datapoints for each dimension
                            // must satisfy Nd%2==0
    const double maxr = R_RANGE;   // Maximum length of vector component
    delta = maxr/Nd*2.0;

    data = new std::complex<double>[Nd*Nd*Nd*Nd];
    transform = new std::complex<double>[Nd*Nd*Nd*Nd];

    // Initialize FFTW
    const int dim[4] = {Nd, Nd, Nd, Nd};
    fftw_plan plan = fftw_plan_dft(4, dim, reinterpret_cast<fftw_complex*>(data),
        reinterpret_cast<fftw_complex*>(transform), -1, FFTW_ESTIMATE);

    
    N->SetOutOfRangeErrors(false);
    for (int i=0; i<Nd*Nd*Nd*Nd; i++)
    {
        if (READFILE) data[i]=0;
        else data[i]=-99999;
    }

        // Fill array
        // Coordinates: 0=r_x, 1=r_y, 2=v1_x, 3=v1_y
    
    if (READFILE)
    {
        cerr <<"# Loading file..." << endl;
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
        cerr << "# File loaded" << endl;


        ///DEBUG
        // set v1=(4*delta, -10*delta)
        //int v1i_x=4; int v1i_y=-10+Nd;

        ///DEBUG
        {
        int vyind=3; int vxind=48;
        
        double vx = vxind*delta;
        if (vxind > Nd/2)
            vx = -(Nd-vxind)*delta;
        double vy = vyind*delta;
        if (vyind > Nd/2)
            vy = -(Nd-vyind)*delta;

        cout << "# vx: "  << vx <<" vy: " << vy << endl;
        for 
        (int ryind=0; ryind<Nd; ryind++)
        {
        for (int rxind=0; rxind<Nd; rxind++)
         {   
            double rx = rxind*delta;
            if (rxind > Nd/2)
                rx = -(Nd-rxind)*delta;
            double ry = ryind*delta;
            if (ryind > Nd/2)
                ry = -(Nd-ryind)*delta;


            int index = vyind + Nd*(vxind + Nd*(ryind + Nd*rxind) );
            double result=data[index].real();
            //double result = V2int(rx, ry, vx, vy, ya, z);

            cout << rx << " " << ry << " " << result << " " << data[index].imag() << endl;    
        }
        } }
        exit(1);
        
    }   // end if datafile==true
    
    else
    {
        // Calculate
        int done=0;
        cout <<"#Nd=" << Nd << " delta=" << delta << " max_component = " << maxr << endl;
        cout <<"# index   result   rx  ry  vx  vy " << endl;
        
        #pragma omp parallel for schedule(dynamic,50)
        for (int index=0; index<Nd*Nd*Nd*Nd; index++)
        {
            int vyind = index%Nd;
            int vxind = ((index-vyind)/Nd)%Nd;
            int ryind = ((index-vyind-Nd*vxind)/(Nd*Nd))%Nd;
            int rxind = ((index-vyind-Nd*vxind-Nd*Nd*ryind)/(Nd*Nd*Nd))%Nd;

            if (index != vyind + Nd*vxind + Nd*Nd*ryind + Nd*Nd*Nd*rxind)
                cerr << "Indexes don't match! " << LINEINFO << endl;
        //for (int rxind=0; rxind<Nd; rxind++)
        //{
            double rx = rxind*delta;
            if (rxind > Nd/2)
                rx = -(Nd-rxind)*delta;
            
            //for (int ryind=0; ryind<Nd; ryind++)
            //{
                double ry = ryind*delta;
                if (ryind > Nd/2)
                    ry = -(Nd-ryind)*delta;
                
               // for (int vxind=0; vxind<Nd; vxind++)
                //{
                    double vx = vxind*delta;
                    if (vxind > Nd/2)
                        vx = -(Nd-vxind)*delta;
                    
                   // for (int vyind=0; vyind<Nd; vyind++)
                    //{
                        double vy = vyind*delta;
                        if (vyind > Nd/2)
                            vy = -(Nd-vyind)*delta;
                        
                        

                        //int index = vyind + Nd*(vxind + Nd*(ryind + Nd*rxind) );
                                            
                        
                        // Mirror: x components -> -x
                        
                        int rxi_xm = rxind;
                        int vxi_xm = vxind;
                        if (rxind != 0) rxi_xm = Nd-rxind;
                        if (vxind != 0) vxi_xm = Nd-vxind;
                        int index2 = vyind + Nd*(vxi_xm + Nd*(ryind + Nd*rxi_xm) );

                        // Mirror y
                        int ryi_ym = ryind;
                        int vyi_ym = vyind;
                        if (ryind != 0) ryi_ym = Nd-ryind;
                        if (vyind != 0) vyi_ym = Nd-vyind;

                        int index3 = vyi_ym + Nd*(vxind + Nd*(ryi_ym + Nd*rxind) );

                        // Mirror both
                        int index4 = vyi_ym + Nd*(vxi_xm + Nd*(ryi_ym + Nd*rxi_xm) );
                        
                        if (data[index].real()>-99900 and data[index2].real()>-99900
                            and data[index3].real()>-99900 and data[index4].real()>-99900)
                        {
                            continue;
                        }

                        double result = V2int(rx, ry, vx, vy, ya, z);

                        
                        
                        #pragma omp critical
                        {
                            data[index]=result;
                            data[index2]=result;
                            data[index3]=result;
                            data[index4]=result;
                            done++;
                            if (data[index].real() != 0 or data[index].imag()!=0)
                            {
                                if (done % 1000 == 0)
                                cerr << " # done " << done << "/" << std::pow(Nd, 4)/4
                                    << " (approximation)" << endl;
                    
                                //<< " prevres " << data[index]
                                //<< " rx " << rx << " ry " << ry << " vx " << vx << " vy " << vy << endl;
                               /*
                                cout << index << std::setprecision(6) << " " << data[index].real() <<
                                " " << rx << " " << ry << " " << vx << " " << vy << endl;
                                cout << index2 << std::setprecision(6)<< " " << data[index].real() << " " << -rx << " "
                                 << ry << " " << -vx << " " << vy << endl;
                                cout << index3 << std::setprecision(6)<< " " << data[index].real() <<
                                " " << rx << " " << -ry << " " << vx << " " << -vy << endl;
                                cout << index4 << std::setprecision(6)<< " " << data[index].real() <<
                                " " << -rx << " " << -ry << " " << -vx << " " << -vy << endl;
                                */
                            }
                        } // end omp critical
                   // }
                //}
            //} // end for ryind
        } 


        // Consistency check
        for (int i=0; i<Nd*Nd*Nd*Nd; i++)
        {
            if (data[i].real() < -90000)
                cerr << "data[" << i <<"] = " << data[i] << endl;
        }
    }   // end if(READFILE)


    cout << "# Transforming...." << endl;
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    delete[] data;

}
double CrossSection2::CorrectionTerm_fft(double pt1, double pt2,
    double ya, double phi, double z)
{
    // Coordinates: 0=r_x, 1=r_y, 2=v1_x, 3=v1_y
    // k (=gluon, pt1) couples with v and r couples with \Delta = k+q = pt1+pt2
    // Let's fix x-axis to be along the k vector, k=(pt1, 0)
    // q = pt2(cos \delphi, sin \delphi)
    // \Delta = (pt1 + pt2 cos \delphi, pt2 sin \delphi)

    if (pt1 > 2.0*M_PI/(2.0*delta) or pt2*std::sin(phi) > 2.0*M_PI/(2.0*delta)
        or pt1 + pt2*std::cos(phi) > 2.0*M_PI/(2.0*delta) )
        cerr << "Too large momenta!!" << endl;
        
    int kyind = 0; double ky=0;
    int kxind = round(Nd*delta*pt1/(2.0*M_PI));
    double kx = (double)kxind / (Nd*delta) * 2.0*M_PI;
    //cout <<"# kx: " << pt1 << " -> " << kx << " kxind " << kxind <<endl;

    int dxind = round(Nd*delta*(pt1 + pt2*std::cos(phi)) / (2.0*M_PI));
    if (dxind<0) dxind += Nd;
    double dx = dxind/(Nd*delta) * 2.0*M_PI;
    if (dxind-1 > Nd/2)
        dx = -(Nd-dxind)/(Nd*delta) * 2.0*M_PI;
    //cout << "#dx: " << pt1 + pt2*std::cos(dphi) << " -> " << dx << " dxind " << dxind <<endl;

    int dyind = round(Nd*delta*(pt2*std::sin(phi)) / (2.0*M_PI) );
    if (dyind<0) dyind +=Nd;
    double dy = dyind/(Nd*delta) * 2.0*M_PI;
    if (dyind-1 > Nd/2)
        dy = -(Nd-dyind)/(Nd*delta) * 2.0*M_PI;
    //cout << "#dy: " << pt2*std::sin(dphi) << " -> " << dy << " dyind " << dyind << endl;


    //int index = vyind + Nd*(vxind + Nd*(ryind + Nd*rxind) );
    int index = kyind + Nd*(kxind + Nd*(dyind + Nd*dxind));

    cerr<< "# phi " << phi << " result " << transform[index]*std::pow(delta,4)/std::pow(2.0*M_PI,6)
        << " |delta|: " << std::sqrt(SQR(pt1)+SQR(pt2) + 2.0*pt1*pt2*std::cos(phi))
        << " -> " << std::sqrt(SQR(dx) + SQR(dy)) << ", k: " << pt1 << " -> " << kx
        << " index " << index << " / " << Nd*Nd*Nd*Nd-1 << endl;
    //cout << phi << " " << transform[index].real()*std::pow(delta,4)/std::pow(2.0*M_PI,6) << endl;
            
    
    

    return transform[index].real()*std::pow(delta,4)/std::pow(2.0*M_PI,6) ;
}

struct Inthelper_v2int
{
    double rx,ry,v1x,v1y;
    double r,s_r;
    double s_r_p_v1;
    double r_p_v1;
    double v2x;
    double y;
    double v2_r;
    double theta;
    double z;
    AmplitudeLib* N;
};
const int VINTINTERVALS = 100;
const double VINTACCURACY = 0.01;
double Inthelperf_v2thetaint(double theta, void* p);
double Inthelperf_v2rint(double lnv2r, void* p);

// Integrate over v2
double CrossSection2::V2int(double rx, double ry, double v1x, double v1y, double y,
    double z)
{
    Inthelper_v2int helper;
    helper.rx=rx; helper.ry=ry; helper.v1x=v1x; helper.v1y=v1y; helper.y=y;
    helper.N=N; helper.r = std::sqrt(SQR(rx) + SQR(ry) );
    helper.z=z;
    helper.s_r = N->S(helper.r, y);
    helper.r_p_v1 = std::sqrt( SQR(rx+v1x) + SQR(ry+v1y) );
    helper.s_r_p_v1 = N->S(helper.r_p_v1, y);

    ///DEBUG
    /*
    for (double v2=0.001; v2<1e6; v2*=1.2)
    {
        for (double theta=0; theta < 2.0*M_PI; theta+=0.03)
        { 
        helper.rx=-3; helper.ry=4;
        helper.v1x=3.5; helper.v1y=-2.3; helper.r = std::sqrt(SQR(helper.rx) + SQR(helper.ry) );
        helper.v2_r = v2; helper.s_r = N->S(helper.r, y);
        cout <<  v2 << "  " << theta << " " << Inthelperf_v2thetaint(theta, &helper)*v2 << endl;
        //cout << v2 << " " << Inthelperf_v2rint(std::log(v2), &helper)/v2 << endl;
        }   
    }
    exit(1);*/
    
    
    
    gsl_function f;
    f.function = Inthelperf_v2rint;
    f.params = &helper;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(VINTINTERVALS);
    double result,abserr;

    double min = std::log(1e-5);
    double max = std::log(1e8);


    int status = gsl_integration_qag(&f, min, max,
            0, VINTACCURACY, VINTINTERVALS, GSL_INTEG_GAUSS41, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);
    result *= -8.0*SQR(M_PI)*SQR(M_Q)*SQR(z);

    if (status)
        cerr << "v2rint failed at " << LINEINFO <<", result " << result
        << " relerror " << std::abs(abserr/result) << " rx: " << rx << " ry "
        << ry << " v1x " << v1x << " v1y " << v1y << endl;

    return result*helper.s_r_p_v1;
}

double Inthelperf_v2rint(double lnv2, void* p)
{
    Inthelper_v2int* par = (Inthelper_v2int*)p;
    par->v2_r = std::exp(lnv2);

    gsl_function f;
    f.function = Inthelperf_v2thetaint;
    f.params = par;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(VINTINTERVALS);
    double result,abserr;

    double min = 0;
    double max = 2.0*M_PI;

    int status = gsl_integration_qag(&f, min, max,
            0, VINTACCURACY, VINTINTERVALS, GSL_INTEG_GAUSS51, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status and std::abs(result*SQR(par->v2_r))>1e-6)
        cerr << "v2thetaint failed at " << LINEINFO <<", result " << result
        << " relerror " << std::abs(abserr/result) << " rx: " << par->rx << " ry "
        << par->ry << " v1x " << par->v1x << " v1y " << par->v1y << " "
        << " v2 " << std::exp(lnv2) << endl;

    return result*SQR(par->v2_r);
}

double Inthelperf_v2thetaint(double theta, void* p)
{
    Inthelper_v2int* par = (Inthelper_v2int*)p;
    double rx = par->rx; double ry=par->ry;
    double v1x = par->v1x; double v1y = par->v1y;
    double v2x = par->v2_r * std::cos(theta);
    double v2y = par->v2_r * std::sin(theta);
    double rsqr = SQR(rx)+SQR(ry);
    
  
    AmplitudeLib* N = par->N;
    double y = par->y;

    // vector lengths
    // r - v2 + 1/2v1
    double r_m_v2_p_05v1 = std::sqrt(
        SQR(rx - v2x + v1x/2.0) + SQR(ry - v2y + v1y/2.0) );

    // r + 1/2v + v2
    double r_p_05v1_p_v2 = std::sqrt(
        SQR(rx + v1x/2.0 + v2x) + SQR(ry + v1y/2.0 + v2y) );

    // r+v1
    //double r_p_v1 = std::sqrt( SQR(rx+v1x) + SQR(ry+v1y) );

    // 1/2v1 + v2
    double v2_p_05v1 = std::sqrt( SQR(v1x/2.0 + v2x) + SQR(v1y/2.0+v2y) );

    // v2-1/2*v1
    double v2_m_05v1 = std::sqrt( SQR(v1x/2.0 - v2x) + SQR(v1y/2.0 - v2y) );

    // cache the values of S's that appear twise
    // S(r+v)
    //double s_r_p_v1 = N->S(r_p_v1,y);
    // S(v'+1/2v)
    double s_v2_p_05v1 = N->S(v2_p_05v1,y);
    // S(v'-1/2v)
    double s_v2_m_05v1 = N->S(v2_m_05v1,y);
    // S(r)
    double s_r = par->s_r; //N->S(r,y);
    // S(r+v)
    double s_r_p_v1 = par->s_r_p_v1;
    
    double nom1 = N->S(r_m_v2_p_05v1, y)*N->S(r_p_05v1_p_v2, y);
    double denom1 = s_r*s_r_p_v1; //N->S(r,y)*N->S(r_p_v1 ,y);

    REAL nom2 = s_v2_p_05v1 * s_v2_m_05v1; //N->S(v2_p_05v1, y)*N->S(v2_m_05v1, y);
    if (nom2<=0) return 0;
    REAL denom2 = denom1; // optimize, same as before

    double result;

    // F/F
    if (rsqr < SQR(par->N->MinR())) result=1;
    else
        result = std::log(nom1/denom1) / std::log(nom2/denom2);


    result *= ( s_v2_p_05v1*s_v2_m_05v1 - s_r*s_r_p_v1 );
        //N->S(r_p_v1, y) *( N->S(v2_p_05v1, y) * N->S(v2_m_05v1, y)
        // - N->S(r, y)*N->S(r_p_v1, y) );

    // wave function product

    // Factor -8pi^2 m^2 z^2 is outside the integral
    result *= (1.0+SQR(1.0-par->z))
             *(SQR(par->v2_r)-0.25*(SQR(v1x) + SQR(v1y)))
                / (v2_m_05v1 * v2_p_05v1)
             * gsl_sf_bessel_K1(par->z*M_Q*v2_m_05v1)
             * gsl_sf_bessel_K1(par->z*M_Q*v2_p_05v1)
            - SQR(par->z)*gsl_sf_bessel_K0(par->z*M_Q*v2_m_05v1)
             * gsl_sf_bessel_K0(par->z*M_Q*v2_p_05v1)  ;
    
/*    
    result *= -8.0*SQR(M_PI)*(2.0-par->z)*SQR(par->z)*SQR(M_Q)
                * gsl_sf_bessel_K1(par->z*M_Q*v2_m_05v1)
                * gsl_sf_bessel_K1(par->z*M_Q*v2_p_05v1)
                * (SQR(par->v2_r)-0.25*(SQR(v1x) + SQR(v1y)) )
                  / (v2_m_05v1 * v2_p_05v1);
  */      

    // m,z=0
    /*
    result *= -16.0*SQR(M_PI)
        * ( SQR(v2x) + SQR(v2y) - 0.25*(SQR(v1x) + SQR(v1y)) )
        / ( SQR(v2_m_05v1) * SQR(v2_p_05v1) );
    */

    /*
    if (isnan(result))
    {
        cerr << "Result is NaN at " << LINEINFO << endl;
        return 0;
    }
    if (isinf(result))
        cerr << result << " at " << LINEINFO << ", v2m05v1 " << v2_m_05v1 <<
         " " << " v2p05v1 " << v2_p_05v1 << " r " << std::sqrt(rsqr)
         << " v2r " << par->v2_r << " theta " << theta <<  endl;
    */
    
    return result;

}

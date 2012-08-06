/*
 * Calculate correction term from 4-point function using fft
 * Idea: we integrate over r, v1, v2. Integrations over
 * r and v2 are 2D FFT's, so perform them by using FFTW.
 * Basically, we have one 4-dimensional FT, so
 * we first fill a 4-dimensional array and then FT it.
 */

#include "config.hpp"
#ifndef USE_FFTW
	#include "xs.hpp"
	void CrossSection2::CalculateCorrection_fft(double ya, double z)
		{ return; }
	double CrossSection2::CorrectionTerm_fft(double pt1, double pt2,
    double ya, double phi)
		{ return 0; }
#else
	

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

const double R_RANGE = 23;

const bool READFILE = false; // Read integrand from a file

bool fftw_cyrille=true;
bool fftw_correction=false;

void CrossSection2::CalculateCorrection_fft(double ya, double z)
{
    Nd = 42;  // Number of datapoints for each dimension
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
            if (std::abs(data[i])>1e-10 and std::abs(data[i].real()/v-1.0) > 0.001)
                cerr <<"trying to set data[" << i <<"] to " << v << " but it is " << data[i] << endl;
            data[i] = v;
        }
        file.close();
        cerr << "# File loaded" << endl;


        ///DEBUG
        // set v1=(4*delta, -10*delta)
        //int v1i_x=4; int v1i_y=-10+Nd;

        ///DEBUG
        /*
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
        */
    }   // end if datafile==true
    
    else
    {
        // Calculate
        int done=0;
        cout <<"#Nd=" << Nd << " delta=" << delta << " max_component = " << maxr << endl;
        cout <<"# index   result   rx  ry  vx  vy " << endl;
        cout <<"# Marquet:" << fftw_cyrille << " correction: " << fftw_correction << endl;
        
        #pragma omp parallel for //schedule(dynamic,50)
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
                        if (rxind != 0 and vxind != 0 and rxind != Nd/2 and ryind != Nd/2)
                        {
                            rxi_xm = Nd-rxind;
                            vxi_xm = Nd-vxind;
                        }
                        int index2 = vyind + Nd*(vxi_xm + Nd*(ryind + Nd*rxi_xm) );

                        // Mirror y
                        int ryi_ym = ryind;
                        int vyi_ym = vyind;
                        if (ryind != 0 and vyind != 0 and ryind != Nd/2 and vyind != Nd/2)
                        {
                            ryi_ym = Nd-ryind;
                            vyi_ym = Nd-vyind;
                        }

                        int index3 = vyi_ym + Nd*(vxind + Nd*(ryi_ym + Nd*rxind) );

                        // Mirror both
                        int index4=index;
                        if (index3 != index and index2 != index)
                        {
                            index4 = vyi_ym + Nd*(vxi_xm + Nd*(ryi_ym + Nd*rxi_xm) );
                        }

                        
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
                            if (done % 100 == 0)
                                cerr << " # done " << done << "/" << std::pow((float)Nd, 4)/4
                                    << " (approximation)" << endl;
                            if (data[index].real() != 0 or data[index].imag()!=0)
                            {
                                
                    
                                //<< " prevres " << data[index]
                                //<< " rx " << rx << " ry " << ry << " vx " << vx << " vy " << vy << endl;
                               
                                cout << index << std::setprecision(6) << " " << data[index].real() <<
                                " " << rx << " " << ry << " " << vx << " " << vy << endl;
                                if (index2 != index)
                                    cout << index2 << std::setprecision(6)<< " " << data[index].real() << " " << -rx << " "
                                    <<  ry << " " << -vx << " " << vy << endl;
                                if (index3 != index)
                                cout << index3 << std::setprecision(6)<< " " << data[index].real() <<
                                " " << rx << " " << -ry << " " << vx << " " << -vy << endl;
                                if (index4 != index)
                                    cout << index4 << std::setprecision(6)<< " " << data[index].real() <<
                                    " " << -rx << " " << -ry << " " << -vx << " " << -vy << endl;
                                
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
    double ya, double phi)
{
    // Coordinates: 0=r_x, 1=r_y, 2=v1_x, 3=v1_y
    // k (=gluon, pt1) couples with r
    // q (=quark, pt2) couples with (v1-v2-r)
    // Let's fix x axis to be along the k_x => k_y=0
    // Thus q = (q cos \delta, q \sin \delta)

    if (pt1 > 2.0*M_PI/(2.0*delta) or pt2*std::sin(phi) > 2.0*M_PI/(2.0*delta)
        or pt2*std::cos(phi) > 2.0*M_PI/(2.0*delta) )
        cerr << "Too large momenta!!" << endl;
        
    int kyind = 0; 
    int kxind = round(Nd*delta*pt1/(2.0*M_PI));
    double kx = (double)kxind / (Nd*delta) * 2.0*M_PI;
    cout <<"# kx: " << pt1 << " -> " << kx << " kxind " << kxind <<endl;

    int qxind = round(Nd*delta*(pt2*std::cos(phi) + pt2*std::sin(phi)) / (2.0*M_PI));
    if (qxind<0) qxind += Nd;
    double qx = qxind/(Nd*delta) * 2.0*M_PI;
    if (qxind-1 > Nd/2)
        qx = -(Nd-qxind)/(Nd*delta) * 2.0*M_PI;
    cout << "#qx: " << pt2*std::cos(phi)  << " -> " << qx << " qxind " << qxind <<endl;

    int qyind = round(Nd*delta*(pt2*std::sin(phi)) / (2.0*M_PI) );
    if (qyind<0) qyind +=Nd;
    double qy = qyind/(Nd*delta) * 2.0*M_PI;
    if (qyind-1 > Nd/2)
        qy = -(Nd-qyind)/(Nd*delta) * 2.0*M_PI;
    cout << "#qy: " << pt2*std::sin(phi) << " -> " << qy << " qyind " << qyind << endl;


    //int index = vyind + Nd*(vxind + Nd*(ryind + Nd*rxind) );
    int index = kyind + Nd*(kxind + Nd*(qyind + Nd*qxind));

    cerr<< "# phi " << phi << " result " << transform[index]*std::pow(delta,4)/std::pow(2.0*M_PI,6)
        << " |delta|: " << std::sqrt(SQR(pt1)+SQR(pt2) + 2.0*pt1*pt2*std::cos(phi))
        << " -> " << std::sqrt(SQR(qx + kx) + SQR(qy)) << ", k: " << pt1 << " -> " << kx
        << " index " << index << " / " << Nd*Nd*Nd*Nd-1 << endl;
    //cout << phi << " " << transform[index].double()*std::pow(delta,4)/std::pow(2.0*M_PI,6) << endl;
            
    
    

    return transform[index].real()*std::pow(delta,4)/std::pow(2.0*M_PI,6) ;
}

struct Inthelper_v2int
{
    double rx,ry,v2x,v2y;
    double s_v2;
    double r,s_r;
    double y;
    double v1;
    double theta;
    double z;
    AmplitudeLib* N;
    CrossSection2 *xs;
};
const int VINTINTERVALS = 10;
const int VTHETAINERVALS = 13;
const double VINTACCURACY = 0.04;
double Inthelperf_v1thetaint(double theta, void* p);
double Inthelperf_v1rint(double lnv1r, void* p);

// Integrate over v1
double CrossSection2::V2int(double rx, double ry, double v2x, double v2y, double y,
    double z)
{
    Inthelper_v2int helper;
    helper.rx=rx; helper.ry=ry; helper.v2x=v2x; helper.v2y=v2y; helper.y=y;
    helper.s_v2 = N->S( sqrt(SQR(v2x) + SQR(v2y)), y);
    helper.N=N; helper.r = std::sqrt(SQR(rx) + SQR(ry) );
    helper.z=z;
    helper.s_r = N->S(helper.r, y);
    helper.xs=this;

    ///DEBUG
    
    /*for (double v2=0.001; v2<100*R_RANGE; v2*=1.2)
    {
        for (double theta=0; theta < 2.0*M_PI; theta+=0.03)
        { 
        helper.rx=-4; helper.ry=-1;
        helper.v1x=3.5; helper.v1y=-2; helper.r = std::sqrt(SQR(helper.rx) + SQR(helper.ry) );
        helper.v2_r = v2; helper.s_r = N->S(helper.r, y);
        cout <<  v2 << "  " << theta << " " << Inthelperf_v2thetaint(theta, &helper)*v2 << endl;
        //cout << v2 << " " << Inthelperf_v2rint(std::log(v2), &helper)/v2 << endl;
        }   
    }
    exit(1);
    */
    
    
    gsl_function f;
    f.function = Inthelperf_v1rint;
    f.params = &helper;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(VINTINTERVALS);
    double result,abserr;

    double min = std::log(1e-3);
    double max = std::log(10.0*R_RANGE); // 1gev: 8*


    int status = gsl_integration_qag(&f, min, max,
            0, VINTACCURACY, VINTINTERVALS, GSL_INTEG_GAUSS51, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);
    result *= 8.0*SQR(M_PI)*SQR(M_Q())*SQR(z);

    /*if (status)
        cerr << "v2rint failed at " << LINEINFO <<", result " << result
        << " relerror " << std::abs(abserr/result) << " rx: " << rx << " ry "
        << ry << " v1x " << v1x << " v1y " << v1y << endl;
*/
    return result;
}

double Inthelperf_v1rint(double lnv1, void* p)
{
    Inthelper_v2int* par = (Inthelper_v2int*)p;
    par->v1 = std::exp(lnv1);

    gsl_function f;
    f.function = Inthelperf_v1thetaint;
    f.params = par;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(VTHETAINERVALS);
    double result,abserr;

    double min = 0;
    double max = 2.0*M_PI;

    int status = gsl_integration_qag(&f, min, max,
            0, VINTACCURACY, VTHETAINERVALS, GSL_INTEG_GAUSS51, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);
/*
    if (status and std::abs(result*SQR(par->v2_r))>1e-6)
        cerr << "v2thetaint failed at " << LINEINFO <<", result " << result
        << " relerror " << std::abs(abserr/result) << " rx: " << par->rx << " ry "
        << par->ry << " v1x " << par->v1x << " v1y " << par->v1y << " "
        << " v2 " << std::exp(lnv2) << endl;
*/
    return result*SQR(par->v1);
}

// smooth step function
double step(double x);

double Inthelperf_v1thetaint(double theta, void* p)
{
    Inthelper_v2int* par = (Inthelper_v2int*)p;
    double M_Q = par->xs->M_Q();
    double rx = par->rx; double ry=par->ry;
    double v2x = par->v2x; double v2y = par->v2y;
    double v1x = par->v1 * std::cos(theta);
    double v1y = par->v1 * std::sin(theta);
    double v1 = std::sqrt(SQR(v1x)+SQR(v1y));
    
    //double rsqr = SQR(rx)+SQR(ry);
    
  
    AmplitudeLib* N = par->N;
    double y = par->y;


	double result=0;
	
	// u'+u-r
	double v1_p_v2_m_r = std::sqrt(SQR(v1x + v2x - rx) + SQR(v1y + v2y - ry));
	double s_v1_p_v2_m_r = par->N->S(v1_p_v2_m_r, y);
	
	// u1 + u2
	double v1_p_v2 = std::sqrt(SQR(v1x+v2x) + SQR(v1y+v2y));
	double s_v1_p_v2 = par->N->S(v1_p_v2, y);
	
	// r - u1
	double r_m_v1 = std::sqrt(SQR(rx - v1x) + SQR(rx - v1y));
	double s_r_m_v1 = par->N->S(r_m_v1, y);
	
	double s_r = par->s_r;
	double s_v1 = par->N->S(v1, y);
	double s_v2 = par->s_v2;
	
	
	double f1 = s_r_m_v1 * s_v1_p_v2 / (s_v2 * s_v1);
	double f2 = s_v1 * s_v1_p_v2_m_r / (s_v2 * s_v1);
	
	double loglog = f1/f2;
	if (isnan(loglog))
	{
		if (s_v1 * s_v1_p_v2_m_r > 0 and s_r*s_v2<1e-10)
			loglog=1;
		else
			loglog=0;
	}
	else if (isinf(loglog)) loglog=0;
	result = s_r * s_v1 * s_v1_p_v2_m_r - s_r*loglog*( s_v1 * s_v1_p_v2_m_r - s_r*s_v2);
	
	result -= s_v1 * s_v1_p_v2_m_r;
	
	result -= step(v1_p_v2_m_r - 1.0/LAMBDAQCD)*step(v1 - 1.0/LAMBDAQCD)
			* s_r * s_v2;
	


    // wave function product

    // Factor 8pi^2 m^2 z^2 is outside the integral
    // bessels for efficiency
    double v1_bessel[2];
    gsl_sf_bessel_Kn_array(0,1,par->z*M_Q*v1, v1_bessel);
    double v2_bessel[2];
    gsl_sf_bessel_Kn_array(0,1,par->z*M_Q*v1_p_v2_m_r, v2_bessel);
    
    result *= (1.0+SQR(1.0-par->z))
             *(v1x*(v1x + v2x - rx) + v1y*(v1y + v2y - ry)) / (v1 * v1_p_v2_m_r)
             * v1_bessel[1] // gsl_sf_bessel_K1(par->z*M_Q*v1)
             * v2_bessel[1] //gsl_sf_bessel_K1(par->z*M_Q*(v1+v2-r))
            + SQR(par->z)* v1_bessel[0] //gsl_sf_bessel_K0(par->z*M_Q*(v1))
             * v2_bessel[0]; //gsl_sf_bessel_K0(par->z*M_Q*(v1+v2-r)) 
       

    // m,z=0
    /*
    result *= 16.0*SQR(M_PI)
        * ( SQR(v2x) + SQR(v2y) - 0.25*(SQR(v1x) + SQR(v1y)) )
        / ( SQR(v2_m_05v1) * SQR(v2_p_05v1) );
    */

    
    if (isnan(result))
    {
        cerr << "Result is NaN at " << LINEINFO << ", loglog=" << loglog << " besse1 " << v1_bessel[0] << endl;
        return 0;
    }
    if (isinf(result))
        cerr << result << " at " << LINEINFO  << endl;
    
    
    return result;

}

#endif //ifdef USE_MPI

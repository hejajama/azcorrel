#include "xs.hpp"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <cmath>

using std::cos;
using std::sin;

/*
 * Calculate the additional term to two-body correlations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

struct Inthelper_correction
{
    double pt1,pt2,ya;
    double phi;
    AmplitudeLib* N;
    size_t calln; int monte;

    double u1,u2,r,theta1,theta2,thetar,z;
};

double Inthelperf_correction(double* vec, size_t dim, void* p);
double NormalizeAngle(double phi);

double Lx(double x, double y, AmplitudeLib* N);
size_t calls = 5e9;
double CrossSection2::CorrectionTerm(double pt1, double pt2, double ya, double phi,
        double z)

{

    /*
     * We need to integrate over u,u',r and their corresponding angles
     */
    // vec[0]= ln u1, vec[1]=ln u2, vec[2]= ln
    // vec[3]=theta_u1, vec[4]=theta_u2
    // vec[5]=theta_r
    // pt1=k, pt2=q

    N->SetOutOfRangeErrors(false);
    //double minlnr = std::log(1e-5);
    //double maxlnr = std::log(1e8);
    double minr=std::log(0.001); double maxr=std::log(10000);
    

    cerr << "Calculating correction at phi=" << phi <<", z=" << z
    << " ya=" << ya <<", pt1=" << pt1 << " pt2=" << pt2 <<" M_Q=" << M_Q << "GeV" <<
    " maxr " << std::exp(maxr) << " calls " << calls << endl;

    

    ///debug
    /*
    for (double r=0.01; r<20; r+=0.2)
    {
        for (double thetar=0; thetar<2.0*M_PI; thetar+=0.2)
        {
            
            double u1=1.6; double u2=1.2; double theta1=0.5; double theta2=1;
            double vec[6]={u1, u2, std::log(r), theta1, theta2, thetar};
            cout << r << " " << thetar << " " << Inthelperf_correction(vec, 6, &helper) << endl;
        }
    }
    exit(1);
    */
    

    


    int cube=0;
    // Parallel: divide cube into 8 parts
    /*
    #pragma omp parallel for
    for (int cube=1; cube<=8; cube++)
    {
        */
        const gsl_rng_type *T = gsl_rng_default;
        gsl_rng_env_setup();
        //gsl_rng *r = gsl_rng_alloc(T);
        
        /*double lower[6] = {minlnr, minlnr, minlnr, 0, 0, 0};
        double upper[6] = {maxlnr, maxlnr,maxlnr, 2.0*M_PI, 2.0*M_PI, 2.0*M_PI};
        */
        
        /*
        // cubes 1,3,5,7: u1 starts from log(MinR())
        if (cube%2==0)  
            lower[0]=0.5*(minlnr+maxlnr);
        else
            upper[0]=0.5*(minlnr+maxlnr);

        // for 1,2,5,6 u2: minlnr->
        if (cube==1 or cube==2 or cube==5 or cube==6)
            upper[1]=0.5*(minlnr+maxlnr);
        else
            lower[1]=0.5*(minlnr+maxlnr);

        // for 1,2,3,4 r:minlnr->
        if (cube<=4)
            upper[2] = 0.5*(minlnr+maxlnr);
        else
            lower[2] = 0.5*(minlnr+maxlnr);
    */
       
        cout << "-----------" << endl;
    double res;
    // Constants taken out from integral in order to reduse number of
    // products in integrand
    double constants = -8.0*SQR(M_PI)*SQR(z)*SQR(M_Q)/std::pow(2.0*M_PI, 6);
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            Inthelper_correction helper;
            helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
            helper.N=N; helper.z=z; helper.calln=0; helper.monte=1;
            gsl_monte_function G = {&Inthelperf_correction, 6, &helper};
    
            double result,abserr;
            double lower[6] = {minr, minr, minr, 0, 0, 0};
            double upper[6] = {maxr, maxr,maxr, 2.0*M_PI, 2.0*M_PI, 2.0*M_PI};
            gsl_rng *r = gsl_rng_alloc(T);
            gsl_monte_plain_state *s = gsl_monte_plain_alloc (6);
            gsl_monte_plain_integrate (&G, lower, upper, 6, calls, r, s, 
                                        &result, &abserr);
            gsl_monte_plain_free (s);
            gsl_rng_free(r);
            result *= constants;
            cout << "# plain result: " << result << " relerr " << std::abs(abserr/result) << endl;
       }
        #pragma omp section
        {
            Inthelper_correction helper;
            helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
            helper.N=N; helper.z=z; helper.calln=0; helper.monte=2;
            gsl_monte_function G = {&Inthelperf_correction, 6, &helper};
            
            double result,abserr;
            double lower[6] = {minr, minr, minr, 0, 0, 0};
            double upper[6] = {maxr, maxr,maxr, 2.0*M_PI, 2.0*M_PI, 2.0*M_PI};
            gsl_rng *r = gsl_rng_alloc(T);
             gsl_monte_miser_state *s = gsl_monte_miser_alloc (6);
             gsl_monte_miser_integrate (&G, lower, upper, 6, calls, r, s, 
                                        &result, &abserr);
             gsl_monte_miser_free (s);
            gsl_rng_free(r);
            result *= constants;    // integration measures
             cout <<"# miser result: " << result << " relerr " << std::abs(abserr/result) << endl;
        }
        #pragma omp section
        {
            Inthelper_correction helper;
            helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
            helper.N=N; helper.z=z; helper.calln=0; helper.monte=3;
            gsl_monte_function G = {&Inthelperf_correction, 6, &helper};
            
            double result,abserr;
            double lower[6] = {minr, minr, minr, 0, 0, 0};
            double upper[6] = {maxr, maxr,maxr, 2.0*M_PI, 2.0*M_PI, 2.0*M_PI};
            gsl_rng *r = gsl_rng_alloc(T);
            gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (6);
         
            gsl_monte_vegas_integrate (&G, lower, upper, 6, calls/5, r, s,
                                        &result, &abserr);
            result *= constants;    // integration measures
             cout << "# cube " << cube << " Vegas warmup: " << result << " relerr " << std::abs(abserr/result) << endl;
            int i=0;
            do
            {
                helper.calln=0;
                gsl_monte_vegas_integrate (&G, lower, upper, 6, calls/2, r, s,
                                        &result, &abserr);
                result *= constants;   // integration measures
                cout << "cube " << cube << " result = " << result << " relerr = " << std::abs(abserr/result)
                         << " chisq/dof = " << gsl_monte_vegas_chisq (s) << endl;
                i++;
            }
            while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5 or i<5);

            gsl_monte_vegas_free (s);
            gsl_rng_free(r);

            res=result;
        }// end omp section

        /*#pragma omp section
        {
            // own parallel monte carlo
            Inthelper_correction helper;
            helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
            helper.N=N; helper.z=z; helper.calln=0; helper.monte=4;
            gsl_rng *r = gsl_rng_alloc(T);
            double sum=0;
            double vec[6];
            #pragma omp parallel for private(vec) reduction(+:sum)
            for (size_t i=0; i<calls; i++)
            {
                vec[0]=gsl_rng_uniform(r)*2.0*M_PI;
                vec[1]=gsl_rng_uniform(r)*2.0*M_PI;
                vec[2]=gsl_rng_uniform(r)*2.0*M_PI;
                vec[3] = minr + gsl_rng_uniform(r)*(maxr-minr);
                vec[4] = minr + gsl_rng_uniform(r)*(maxr-minr);
                vec[5] = minr + gsl_rng_uniform(r)*(maxr-minr);
                double tmp = Inthelperf_correction(vec, 6, &helper);
                sum = sum+tmp;
            }
            gsl_rng_free(r);
            cout <<"# Parallel: " << sum*std::pow(2.0*M_PI*(maxr-minr),3)/calls << endl;
        }*/


        
   } // end omp sections
    //gsl_rng_free(r);
   // }


      


    

    return res;
}


double Inthelperf_correction(double* vec, size_t dim, void* p)
{
    
    Inthelper_correction* par = (Inthelper_correction*)p;
    #pragma omp critical
    par->calln++;
    if (par->calln%(calls/50)==0) cout <<"# Monte " << par->monte << " done " << (double)par->calln/calls << endl;

    // vec[0]= ln u1, vec[1]=ln u2, vec[2]= ln
    // vec[3]=theta_u1, vec[4]=theta_u2
    // vec[5]=theta_r
    // angles are within range [0, 2\pi]
    double u1 = std::exp(vec[0]); double u2=std::exp(vec[1]); double r=std::exp(vec[2]);
    //double u1=vec[0]; double u2=vec[1]; double r=vec[2];
    double theta1=vec[3];
    double theta2=vec[4]; double thetar = vec[5];

    AmplitudeLib* N = par->N;

    double pt1 = par->pt1;
    double pt2 = par->pt2;
    double phi = par->phi;
    double y = par->ya;

    // Dot products
    // pt1=k, pt2=q

    double pt1_dot_u1 = pt1*u1*cos(theta1);
    double pt1_dot_u2 = pt1*u2*cos(theta2);
    double pt1_dot_r  = pt1*r*cos(thetar);
    //double pt2_dot_u1 = pt2*u1*cos(NormalizeAngle(theta1 - phi));
    //double pt2_dot_u2 = pt2*u2*cos(NormalizeAngle(theta2 - phi));
    double pt2_dot_r  = pt2*r*cos(thetar - phi);
    double u1_dot_r   = u1*r*cos(theta1-thetar);
    double u2_dot_r   = u2*r*cos(theta2-thetar);
    double u1_dot_u2  = u1*u2*cos(theta2-theta1);

    // Norms
    // (u1-u2+r)^2
    double u1u2r_sqr = SQR(u1)+SQR(u2)+SQR(r) + 2.0*u1_dot_r - 2.0*u2_dot_r
        - 2.0*u1_dot_u2;
    double u1u2r = std::sqrt(u1u2r_sqr);
    // (u1+r)^2
    double u1r_sqr = SQR(u1)+SQR(r) + 2.0*u1_dot_r;
    double u1r = std::sqrt(u1r_sqr);
    // (u2-r)^2
    double u2r_sqr = SQR(u2)+SQR(r) - 2.0*u2_dot_r;
    double u2r = std::sqrt(u2r_sqr);

    // S(u1-u2+r)
    double s_u1u2r = N->S(u1u2r,y);
    // S(u1)
    double s_u1 = N->S(u1,y);
    // S(u2)
    double s_u2 = N->S(u2,y);
    // S(r)
    double s_r = N->S(r,y);

    // Contribution from e^(ik(u'-u)) e^(-i \Delta r),
    // k = pt1, \delta = pt1+pt2
    // Take only real part
    //double result = cos(pt1_dot_u2 - pt1_dot_u1)*cos(pt1_dot_r + pt2_dot_r)
   //     + sin(pt1_dot_u2 - pt1_dot_u1)*sin(pt1_dot_r + pt2_dot_r);

    // SHOULD go to zero with this
    double result = sin(pt1_dot_u2 - pt1_dot_u1)*cos(pt1_dot_r+pt2_dot_r)
        - cos(pt1_dot_u2 - pt1_dot_u1)*sin(pt1_dot_r+pt2_dot_r);

    //result *= N->S(u1u2r, y) *
    //    ( N->S(u1, y) * N->S(u2,y) - N->S(r, y) * N->S(u1u2r, y) );
    result *= s_u1u2r * (s_u1 * s_u2 - s_r * s_u1u2r);

    // Finally multiply by complicated F/F factor
    /*REAL f1 = Lx(u2r,y,N) - Lx(r,y,N) + Lx(u1r,y,N) - Lx(u1u2r,y,N);
    REAL f2 = Lx(u1,y,N) - Lx(r,y,N) + Lx(u2,y,N) - Lx(u1u2r,y,N);
    result *= f1/f2;
    */

    REAL minsprod = 1e-15;  /// result should not depend on this!
    REAL nom1 = N->S(u2r, y)*N->S(u1r, y);
    REAL denom1 = s_r * s_u1u2r; //N->S(r,y)*N->S(u1u2r,y);
    if (nom1 < minsprod or denom1 < minsprod) return 0;

    REAL nom2 = s_u1*s_u2; //N->S(u1, y)*N->S(u2, y);
    REAL denom2 = denom1; // optimize, same as before, N->S(r,y)*N->S(u1u2r,y);
    if (nom2 < minsprod or denom2 < minsprod or
        std::abs(nom2/denom2-1.0) < minsprod ) return 0;

    REAL f1f2 = std::log(nom1/denom1) / std::log(nom2/denom2);

    result *= f1f2;

    // Wave function product
    // \phi^*(u2) \phi(u) summed over spins and polarization
    // simple resul in massless z=0 case
    //result *= -16.0*SQR(M_PI) * u1_dot_u2 / ( SQR(u1)*SQR(u2) );
    // with masses and z!=0 a bit more complicated:
    // -8pi^2 m^2 z^2 [ K_1(mzu) K_1(mzu') u.u'/(uu') [1+(1-z)^2]
    //                - z^2 K_0(mzu)K_0(mzu') ]
    // In order to optimize, factor -8*pi^2*m^2*z^2
    // is outside the integral
    result *= (1.0+SQR(1.0-par->z))*u1_dot_u2 /(u1*u2)
                * gsl_sf_bessel_K1(par->z*M_Q*u1)
                * gsl_sf_bessel_K1(par->z*M_Q*u2)
              - SQR(par->z)*gsl_sf_bessel_K0(par->z*M_Q*u1)
                * gsl_sf_bessel_K0(par->z*M_Q*u2)  ;

    result *= u1*u2*r;   // d^2 u_1 d^2 u2 d^2 r = du_1 du_2 dr u_1 u_2 u_3 d\theta_i

    result *= u1*u2*r;  // multiply by e^(ln u1) e^(ln u2) e^(ln r) due to the
                        // change of variables (Jacobian)

    //cout << "u1 " << u1 << " u2 " << u2 << " r " << r << " result " << result << endl;

    //if (isnan(result))
   //  cerr <<" nan at " << "u1 " << u1 << " u2 " << u2 << " r " << r  << endl;

    return result;
}

const int RINTPOINTS = 4;
const int THETAINTPOINTS = 2;
const double MINLNR = std::log(1e-3);
const double MAXLNR = std::log(50);

double Inthelperf_thetar(double thetar, void* p)
{
    Inthelper_correction* par = (Inthelper_correction*)p;

    double vec[6] = {par->u1, par->u2, par->r, par->theta1, par->theta2,
        thetar };
    
    return Inthelperf_correction(vec, 6, par);
}

double Inthelperf_theta2(double theta2, void* p)
{
    Inthelper_correction* par = (Inthelper_correction*)p;
    par->theta2=theta2;

    gsl_function fun;
    fun.params = par;
    fun.function = Inthelperf_thetar;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);
    double result,abserr;

    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI,
            0, 0.3, THETAINTPOINTS, GSL_INTEG_GAUSS31, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    /*if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl;*/
    return result;
}

double Inthelperf_theta1(double theta1, void* p)
{
    Inthelper_correction* par = (Inthelper_correction*)p;
    par->theta1=theta1;

    gsl_function fun;
    fun.params = par;
    fun.function = Inthelperf_theta2;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);
    double result,abserr;

    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI,
            0, 0.3, THETAINTPOINTS, GSL_INTEG_GAUSS31, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    /*if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl;
    */
    return result;
}
int ri=0;
double Inthelperf_r(double r, void* p)
{
    ri++;
    cout << "# ri: " << ri << endl;
    Inthelper_correction* par = (Inthelper_correction*)p;
    par->r=r;

    gsl_function fun;
    fun.params = par;
    fun.function = Inthelperf_theta1;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);
    double result,abserr;

    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI,
            0, 0.3, THETAINTPOINTS, GSL_INTEG_GAUSS31, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl;
    return result;
}

double Inthelperf_u2(double u2, void* p)
{
    ri=0;
    Inthelper_correction* par = (Inthelper_correction*)p;
    par->u2=u2;

    gsl_function fun;
    fun.params = par;
    fun.function = Inthelperf_r;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);
    double result,abserr;

    int status = gsl_integration_qag(&fun, MINLNR, MAXLNR,
            0, 0.3, RINTPOINTS, GSL_INTEG_GAUSS31, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl;
    return result;
}
int u1i=0;
double Inthelperf_u1(double u1, void* p)
{
    u1i++;
    cout << "# u1: " << u1i << "/100" << endl;
    Inthelper_correction* par = (Inthelper_correction*)p;
    par->u1=u1;

    gsl_function fun;
    fun.params = par;
    fun.function = Inthelperf_u2;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);
    double result,abserr;

    int status = gsl_integration_qag(&fun, MINLNR, MAXLNR,
            0, 0.3, RINTPOINTS, GSL_INTEG_GAUSS31, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl;
    return result;
}

double CrossSection2::CorrectionTerm_nomc(double pt1, double pt2, double ya, double phi)
{
    /*
     * We need to integrate over u,u',r and their corresponding angles
     */
    // pt1=k, pt2=q

    N->SetOutOfRangeErrors(false);

    Inthelper_correction helper;
    helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
    helper.N=N;

    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_u1;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);
    double result,abserr;

    int status = gsl_integration_qag(&fun, MINLNR, MAXLNR,
            0, 0.3, RINTPOINTS, GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl;

    return result;
}

double NormalizeAngle(double phi)
{
    //return phi;
    /*
     * An angle between two vectors is allways between [0,pi]
     * phi is calculated by directly substracting the direction of one
     * vector from the direction of an another vector
     * This method forces this difference within an interval [0,pi]
     * Algorighm:
     * If angle > pi => angle = 2*pi-angle
     * if angle<0 angle = -angle
     *
     * Assumes that phi \in [-2\pi, 2\pi]
     */

     if (phi<0) phi=-phi;
     if (phi>M_PI) phi = 2.0*M_PI - phi;

     return phi;
     

}

double Lx(double x, double y, AmplitudeLib* N)
{
    REAL s = N->S(x,y);
    if (s<=0) return -999; ///TODO: result should not be sensitive to this
    REAL result = std::log(s);
    return result;
}

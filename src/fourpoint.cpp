#include "xs.hpp"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include <ctime>
#include "config.hpp"
#ifdef USE_MPI
    #include <mpi.h>

extern "C"
{
    #include "pvegas/vegas.h"
}
#endif

using std::cos;
using std::sin;

/*
 * Calculate the additional term to two-body correlations
 * d\sigma / (dpt1 dpt2 dy1 dy2 d\theta)
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

struct Inthelper_correction
{
    double pt1,pt2,ya;
    double phi;
    AmplitudeLib* N;
    size_t calln; int monte;
    CrossSection2* xs;

    double u1,u2,r,theta1,theta2,thetar,z;
};

double Inthelperf_correction(double* vec, size_t dim, void* p);
double Inthelperf_correction2(double* vec, size_t dim, void* p);
void pVegas_wrapper(double x[6], double f[1], void* par);
double NormalizeAngle(double phi);

bool cyrille=true;
bool correction=false;

/*
 * Calculates the correction term/full integral over 6-dim. hypercube
 */

double CrossSection2::CorrectionTerm(double pt1, double pt2, double ya, double phi,
        double z)

{
    std::time_t start = std::time(NULL);
    #ifdef USE_MPI
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if(id==0) 
    #endif
    cout << "# Starting MC integral, marquet: " << cyrille << ", correction: " << correction << endl;


    if (!N->InterpolatorInitialized(ya))
        cerr << "Interpolator is not initialized for y=" << ya << "!!!" << endl;
    
    /*
     * We need to integrate over u,u',r and their corresponding angles
     */
    // vec[0]= ln u1, vec[1]=ln u2, vec[2]= ln r
    // vec[3]=theta_u1, vec[4]=theta_u2
    // vec[5]=theta_r
    // pt1=k, pt2=q

    N->SetOutOfRangeErrors(false);
    double minr=std::log(0.005); double maxr=std::log(10000);
    //double minr=0; double maxr=2000;

    double (*integrand)(double*,size_t,void*) = Inthelperf_correction2;
    //cout <<"# Integrand is full" << endl;

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
    

    size_t calls = mcintpoints;

    const gsl_rng_type *T = gsl_rng_default;
    //gsl_rng_env_setup();

       
    //cout << "-----------" << endl;
    double res=0;
    // Constants taken out from integral in order to reduse number of
    // products in integrand
    double constants = 8.0*SQR(M_PI)*SQR(z)*SQR(M_Q())/std::pow(2.0*M_PI, 6);
    
    //#pragma omp parallel sections
    //{
        /*
        #pragma omp section
        {
            Inthelper_correction helper;
            helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
            helper.N=N; helper.z=z; helper.calln=0; helper.monte=1;
            gsl_monte_function G = {integrand, 6, &helper};
    
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
            cout << "# " << phi << " plain result: " << result << " relerr " << std::abs(abserr/result*constants) << endl;
       }*/
       
        //#pragma omp section
        /*
        {
            Inthelper_correction helper;
            helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
            helper.N=N; helper.z=z; helper.calln=0; helper.monte=2;
            helper.xs=this;
            gsl_monte_function G = {integrand, 6, &helper};
            
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
             cout <<"# " << phi << " miser result: " << result << " relerr " << std::abs(abserr/result*constants) << endl;
             return result;
        }
        
         */
        //#pragma omp section
        //

        Inthelper_correction helper;
            helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
            helper.N=N; helper.z=z; helper.calln=0; helper.monte=3;
            helper.xs=this;
        
        //{   
            #ifndef USE_MPI
            gsl_monte_function G = {integrand, 6, &helper};
            
            double result,abserr;
            double lower[6] = {minr, minr, minr, 0, 0, 0};
            double upper[6] = {maxr, maxr,maxr, 2.0*M_PI, 2.0*M_PI, 2.0*M_PI};
            gsl_rng *r = gsl_rng_alloc(T);
            gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (6);
         
            gsl_monte_vegas_integrate (&G, lower, upper, 6, calls/5, r, s,
                                        &result, &abserr);
            result *= constants;    // integration measures
            //#pragma omp critical
            //cout << "# " << phi << " Vegas warmup: " << result << " relerr " << std::abs(abserr/result*constants) << endl;
            int i=0;
            double prevres=result;
            do
            {
                prevres=result;
                helper.calln=0;
                std::time_t monte_start = std::time(NULL);
                gsl_monte_vegas_integrate (&G, lower, upper, 6, calls, r, s,
                                        &result, &abserr);
                result *= constants;   // integration measures
                time_t t;
                time(&t);
                #pragma omp critical
                cout << "# " << phi << " z " << z << " vegas iter " << i+1 << " result = " << result << " relerr = " << std::abs(abserr/result*constants)
                         << " change " << std::abs((result-prevres)/prevres)
                         << " chisq/dof = " << gsl_monte_vegas_chisq (s) 
                         << " time " << (std::time(NULL)-monte_start)/60.0/60.0 << " h - "
                         << std::ctime(&t);
                
                i++;

                if (i>5)
                {
                    #pragma omp critical
                    {
                        cerr <<"# !!!!!!!!!!Vegas integral at phi=" << phi
                            <<", z=" << z <<" doesn't converge! "
                            << " iterations took " << (std::time(NULL)-start)/60.0/60.0
                            <<" hours " << LINEINFO << endl;
                    }
                    //return -1;
                    return result;
                }
            }
            while(i<2 or std::abs((result-prevres)/prevres)>0.1
                or (i<3 and std::abs(gsl_monte_vegas_chisq (s)-1.0) > 0.6) );
            // so if chisq is not close to one, require at least 3 iterations

            gsl_monte_vegas_free (s);
            gsl_rng_free(r);

            res=0.5*(result+prevres);
            #pragma omp critical
            {
                cout <<"# " << phi << " z " << z << " MC " << res << " change " << std::abs((result-prevres)/prevres)
                    << " relerror " << std::abs(abserr*constants/result) << endl;
                cout <<"# " << phi << " MC integral with " << calls << " points and " << i << " iterations took " << (std::time(NULL)-start)/60.0/60.0 <<" hours " << endl;
            }
        //}// end omp section
        #endif

    
        
  // } // end omp sections




    /// MPI & PVEGAS
    #ifdef USE_MPI
    double estim[1];   // estimators for integrals 
    double std_dev[1]; // standard deviations              
    double chi2a[1];   // chi^2/n       
    int workers = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &workers);
    workers--;
    if (workers==0)
    {
        cerr << "No workers! MPI needs at least two threads! " << endl;
        return -1;
    }
    
    const int dim=6;
    double reg[2*dim];
    int functions = 1;
    reg[0]=minr; reg[1]=minr; reg[2]=minr;
    reg[dim]=maxr; reg[dim+1]=maxr; reg[dim+2]=maxr;
    reg[3]=0; reg[4]=0; reg[5]=0;
    reg[dim+3]=2.0*M_PI; reg[dim+4]=2.0*M_PI; reg[dim+5]=2.0*M_PI;
    
    // set up the grid (init = 0) with 5 iterations of 1000 samples,
    // no need to compute additional accumulators (fcns = 1),
    // no parallelization yet (wrks = 1).
    vegas(reg, dim, pVegas_wrapper,
        0, mcintpoints/100, 5, NPRN_INPUT | NPRN_RESULT,
        functions, 0, 1,
        estim, std_dev, chi2a, &helper);
      // refine the grid (init = 1) with 5 iterations of 10000 samples,
      // collect in additional accumulators (fcns = FUNCTIONS),
      // two parallel workers (wrks = 2). 
      vegas(reg, dim, pVegas_wrapper,
            1, mcintpoints/10, 5, NPRN_INPUT | NPRN_RESULT,
            functions, 0, 2,
            estim, std_dev, chi2a, &helper);
      // final sample, inheriting previous results (init = 2)
    if (id==0)
      cout <<"# Constants: " << constants << endl;
    vegas(reg, dim, pVegas_wrapper,
            2, mcintpoints, 2, NPRN_INPUT | NPRN_RESULT,
            functions, 0,  workers,
            estim, std_dev, chi2a, &helper);

    if (id==0)
      cout << "# Result: " << estim[0] << " +/- " << std_dev[0] << " time "
       << (std::time(NULL)-start)/60.0/60.0 << " h"<< endl;

    res = estim[0]*constants;
    #endif

    return res;
}










double Inthelperf_correction(double* vec, size_t dim, void* p)
{
    
    Inthelper_correction* par = (Inthelper_correction*)p;
    // vec[0]= ln u1, vec[1]=ln u2, vec[2]= ln r
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
    double z = par->z;

    // Dot products
    // pt1=k, pt2=q
    // choose pt1 along the x axis -> theta(pt1)=0
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


    double result = 0;

    // 4/6-point function
    if (cyrille)
        result = s_u1u2r*s_u1*s_u2; // lowest contribution

    // Correction beyond Marquet's paper
    if (correction)
    {
        double minsprod = 1e-30;  /// result should not depend on this!
        double nom1 = N->S(u2r, y)*N->S(u1r, y);
        double denom1 = s_r * s_u1u2r; //N->S(r,y)*N->S(u1u2r,y);

        double nom2 = s_u1*s_u2; //N->S(u1, y)*N->S(u2, y);
        double denom2 = denom1; // optimize, same as before, N->S(r,y)*N->S(u1u2r,y);
        double f1f2=0;
        if (nom2 <= 0 or denom2 <= 0 or
            std::abs(nom2/denom2-1.0) < minsprod or nom1 <= 0)
                f1f2 = 0;
        else
            f1f2 = std::log(nom1/denom1) / std::log(nom2/denom2);
        result -= f1f2*s_u1u2r*(s_u1*s_u2 - s_r*s_u1u2r);
    }
    // 3-point functions
    // -S(u)S(u+r-zu') - S(u')S(u'-r-zu)
    
    if (cyrille)
    {
        result += -s_u1*N->S(
            sqrt( SQR(u1) + SQR(r) + SQR(z*u2) + 2.0*u1_dot_r - 2.0*z*u1_dot_u2
                - 2.0*z*u2_dot_r) , y)
            - s_u2*N->S(
            sqrt( SQR(u2) + SQR(r) + SQR(z*u1) - 2.0*u2_dot_r + 2.0*z*u1_dot_r
                - 2.0*z*u1_dot_u2) , y);

        // 2-point function
        // +S(r+zu-zu')
        result += N->S( std::sqrt( SQR(r) + SQR(z)*(SQR(u1)+SQR(u2)-2.0*u1_dot_u2)
            + 2.0*z*(u1_dot_r - u2_dot_r)) , y);
    }
  
    
    // Wave function product
    // \phi^*(u2) \phi(u) summed over spins and polarization
    // simple resul in massless z=0 case
    //result *= 16.0*SQR(M_PI) * u1_dot_u2 / ( SQR(u1)*SQR(u2) );
    // with masses and z!=0 a bit more complicated:
    // 8pi^2 m^2 z^2 [ K_1(mzu) K_1(mzu') u.u'/(uu') [1+(1-z)^2]
    //                + z^2 K_0(mzu)K_0(mzu') ]
    // In order to optimize, factor 8*pi^2*m^2*z^2
    // is outside the integral
    // factor 1/k^+ is thrown away as it cancels eventually
    double M_Q = par->xs->M_Q();
    double bessel_u1[2];
    double bessel_u2[2];
    gsl_sf_bessel_Kn_array(0,1,par->z*M_Q*u1, bessel_u1);
    gsl_sf_bessel_Kn_array(0,1,par->z*M_Q*u2, bessel_u2);
    result *= (1.0+SQR(1.0-par->z))*u1_dot_u2 /(u1*u2)
                * bessel_u1[1] // gsl_sf_bessel_K1(par->z*M_Q*u1)
                * bessel_u2[1] //gsl_sf_bessel_K1(par->z*M_Q*u2)
              + SQR(par->z)* bessel_u1[0]//gsl_sf_bessel_K0(par->z*M_Q*u1)
                * bessel_u2[0]; //gsl_sf_bessel_K0(par->z*M_Q*u2)  ;

    // Contribution from e^(ik(u'-u)) e^(-i \Delta r),
    // k = pt1, \delta = pt1+pt2
    // Take only real part
    result *= cos( pt1_dot_u1 - pt1_dot_u2 + pt2_dot_r + pt1_dot_r );

    result *= u1*u2*r;   // d^2 u_1 d^2 u2 d^2 r = du_1 du_2 dr u_1 u_2 u_3 d\theta_i
    
    //double max=500;
    //if ((u1>max or u2>max or r>max) and std::abs(result)>1e-20)
    //    cout << u1 << " " << u2 << " " << r << " " << result << endl;
       
    result *= u1*u2*r;  // multiply by e^(ln u1) e^(ln u2) e^(ln r) due to the
                        // change of variables (Jacobian)

    //cout << vec[0] << " " << vec[1] << " " << vec[2] << " " << result << endl;



    return result;
}













double Inthelperf_correction2(double* vec, size_t dim, void* p)
{
    
    Inthelper_correction* par = (Inthelper_correction*)p;
    // vec[0]= ln u1, vec[1]=ln u2, vec[2]= ln r
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
    double z = par->z;

    // Dot products
    // pt1=k, pt2=q
    // choose pt1 along the x axis -> theta(pt1)=0
    double pt1_dot_u1 = pt1*u1*cos(theta1);
    double pt1_dot_u2 = pt1*u2*cos(theta2);
    double pt1_dot_r  = pt1*r*cos(thetar);
    double pt2_dot_u1 = pt2*u1*cos(theta1 - phi);
    double pt2_dot_u2 = pt2*u2*cos(theta2 - phi);
    double pt2_dot_r  = pt2*r*cos(thetar - phi);
    double u1_dot_r   = u1*r*cos(theta1-thetar);
    double u2_dot_r   = u2*r*cos(theta2-thetar);
    double u1_dot_u2  = u1*u2*cos(theta2-theta1);

    // Amplitudes
    double s_u1=N->S(u1, y);
    double s_u2=N->S(u2,y);



    double result = 0;


    if (cyrille)
    {
        result = (std::cos(-pt1_dot_r - pt2_dot_r - pt2_dot_u2 + pt2_dot_u1)
            * s_u1*s_u2 -
            std::cos(-pt1_dot_r - pt2_dot_r + pt2_dot_u1 + (1.0-z)*pt1_dot_u2-z*pt2_dot_u2)
             * s_u1
            - std::cos(-pt1_dot_r - pt2_dot_r - pt2_dot_u2 - (1.0-z)*pt1_dot_u1 + z*pt2_dot_u1)
             * s_u2
            + std::cos(-pt1_dot_r - pt2_dot_r + (1.0-z)*(pt1_dot_u2 - pt1_dot_u1)
                - z*(pt2_dot_u2 - pt2_dot_u1) )
            ) * N->S(r, y);
    }


    if (correction)
    {
        double correction=0;

        // |r - z(u-u')|
        double r_m_z_u1u2 = std::sqrt( SQR(r) + SQR(z)*(SQR(u1)+SQR(u2)-2.0*u1_dot_u2)
            - 2.0*z*(u1_dot_r - u2_dot_r) );
        double s_r_m_zu1u2 = N->S(r_m_z_u1u2, y);

        // |r + (1-z)(u-u')|
        double r_p_1mz_u1u2 = std::sqrt( SQR(r) + SQR(1.0-z)*(SQR(u1)+SQR(u2)-2.0*u1_dot_u2)
            + 2.0*(1.0-z)*(u1_dot_r - u2_dot_r));
        double s_r_p_1mz_u1u2 = N->S(r_p_1mz_u1u2, y);

        // |r - (1-z)u2 - zu1|
        double r_m_1mz_u2_m_zu1 = std::sqrt( SQR(r) + SQR(z)*SQR(u1) + SQR(1.0-z)*SQR(u2)
            - 2.0*(1.0-z)*u2_dot_r - 2.0*z*u1_dot_r + 2.0*z*(1.0-z)*u1_dot_u2 );
        double s_r_m_1mz_u2_m_zu1 = N->S(r_m_1mz_u2_m_zu1, y);

        // |r + (1-z)u + zu'|
        double r_p_1mz_u_p_zu2 = std::sqrt( SQR(r) + SQR(1.0-z)*SQR(u1)+SQR(z)*SQR(u2)
            + 2.0*(1.0-z)*u1_dot_r + 2.0*z*u2_dot_r + 2.0*z*(1.0-z)*u1_dot_u2 );
        double s_r_p_1mz_u_p_zu2 = N->S(r_p_1mz_u_p_zu2, y);

        double f1 = s_r_m_1mz_u2_m_zu1 * s_r_p_1mz_u_p_zu2 / ( s_r_m_zu1u2 * s_r_p_1mz_u1u2);

        double f2 = s_u1*s_u2 / ( s_r_m_zu1u2 * s_r_p_1mz_u1u2);

        if (std::abs(f2-1.0) < 1e-30 or f1<1e-30 or f2<1e-30 or
            s_r_m_1mz_u2_m_zu1 * s_r_p_1mz_u_p_zu2 < 1e-10 )    // r >> u,u', log/log->0
            correction=0;
        else
        {

            correction = std::log(f1)/std::log(f2) * s_r_p_1mz_u1u2
                *( s_u1*s_u2 - s_r_m_zu1u2 * s_r_p_1mz_u1u2 );

            correction *= std::cos(-pt1_dot_r - pt2_dot_r + (1.0-z)*(pt1_dot_u2 - pt1_dot_u1)
                    - z*(pt2_dot_u2 - pt2_dot_u1) );

            if (isnan(correction) or isinf(correction))
                return 0;
            
            result -=correction;
        }
    }

  
    
    // Wave function product
    // \phi^*(u2) \phi(u) summed over spins and polarization
    // simple resul in massless z=0 case
    //result *= 16.0*SQR(M_PI) * u1_dot_u2 / ( SQR(u1)*SQR(u2) );
    // with masses and z!=0 a bit more complicated:
    // 8pi^2 m^2 z^2 [ K_1(mzu) K_1(mzu') u.u'/(uu') [1+(1-z)^2]
    //                + z^2 K_0(mzu)K_0(mzu') ]
    // In order to optimize, factor 8*pi^2*m^2*z^2
    // is outside the integral
    // factor 1/k^+ is thrown away as it cancels eventually
    double M_Q = par->xs->M_Q();
    double bessel_u1[2];
    double bessel_u2[2];
    gsl_sf_bessel_Kn_array(0,1,par->z*M_Q*u1, bessel_u1);
    gsl_sf_bessel_Kn_array(0,1,par->z*M_Q*u2, bessel_u2);
    result *= (1.0+SQR(1.0-par->z))*u1_dot_u2 /(u1*u2)
                * bessel_u1[1] // gsl_sf_bessel_K1(par->z*M_Q*u1)
                * bessel_u2[1] //gsl_sf_bessel_K1(par->z*M_Q*u2)
              + SQR(par->z)* bessel_u1[0]//gsl_sf_bessel_K0(par->z*M_Q*u1)
                * bessel_u2[0]; //gsl_sf_bessel_K0(par->z*M_Q*u2)  ;


    if (isnan(result))
        cerr << "NAN!\n";

    result *= u1*u2*r;   // d^2 u_1 d^2 u2 d^2 r = du_1 du_2 dr u_1 u_2 u_3 d\theta_i
    
    //double max=500;
    //if ((u1>max or u2>max or r>max) and std::abs(result)>1e-20)
    //    cout << u1 << " " << u2 << " " << r << " " << result << endl;
       
    result *= u1*u2*r;  // multiply by e^(ln u1) e^(ln u2) e^(ln r) due to the
                        // change of variables (Jacobian)

    //cout << vec[0] << " " << vec[1] << " " << vec[2] << " " << result << endl;



    return result;
}




void pVegas_wrapper(double x[6], double f[1], void* par)
{
    f[0] = Inthelperf_correction2(x, 6, par);
}



void CrossSection2::SetMCIntPoints(size_t points)
{
    mcintpoints=points;
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



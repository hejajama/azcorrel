#include "xs.hpp"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
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
    REAL pt1,pt2,ya;
    REAL phi;
    AmplitudeLib* N;
};

REAL Inthelperf_correction(REAL* vec, size_t dim, void* p);
REAL NormalizeAngle(REAL phi);

REAL Lx(REAL x);

REAL CrossSection2::CorrectionTerm(REAL pt1, REAL pt2, REAL ya, REAL phi)

{

    /*
     * We need to integrate over u,u',r and their corresponding angles
     */
    // vec[0]=u1, vec[1]=u2, vec[2]=r, vec[3]=theta_u1, vec[4]=theta_u2
    // vec[5]=theta_r
    // pt1=k, pt2=q
    
    REAL result,abserr;
    REAL lower[6] = {N->MinR(), N->MinR(), N->MinR(), 0, 0, 0};
    REAL upper[6] = {N->MaxR(), N->MaxR(), N->MaxR(), 2.0*M_PI, 2.0*M_PI, 2.0*M_PI};

    N->SetOutOfRangeErrors(false);
    
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(T);

    Inthelper_correction helper;
    helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
    helper.N=N;

    size_t calls = 3000005;

    gsl_monte_function G = {&Inthelperf_correction, 6, &helper};

    
        gsl_monte_plain_state *s = gsl_monte_plain_alloc (6);
        gsl_monte_plain_integrate (&G, lower, upper, 6, calls, r, s, 
                                    &result, &abserr);
        gsl_monte_plain_free (s);
    

    /*
         gsl_monte_miser_state *s = gsl_monte_miser_alloc (6);
         gsl_monte_miser_integrate (&G, lower, upper, 6, calls, r, s, 
                                    &result, &abserr);
         gsl_monte_miser_free (s);
    */


    cout << "# result " << result << " relerror " << std::abs(abserr/result) << endl;    


    gsl_rng_free(r);

    result /= std::pow(2.0*M_PI, 6);

    return result;
}


REAL Inthelperf_correction(REAL* vec, size_t dim, void* p)
{
    Inthelper_correction* par = (Inthelper_correction*)p;

    // vec[0]=u1, vec[1]=u2, vec[2]=r, vec[3]=theta_u1, vec[4]=theta_u2
    // vec[5]=theta_r
    // angles are within range [0, 2\pi]
    REAL u1 = vec[0]; REAL u2=vec[1]; REAL r=vec[2]; REAL theta1=vec[4];
    REAL theta2=vec[5]; REAL thetar = vec[5];

    AmplitudeLib* N = par->N;

    REAL pt1 = par->pt1;
    REAL pt2 = par->pt2;
    REAL phi = par->phi;
    REAL y = par->ya;

    // Dot products
    // pt1=k, pt2=q

    REAL pt1_dot_u1 = pt1*u1*cos(theta1);
    REAL pt1_dot_u2 = pt1*u2*cos(theta2);
    REAL pt1_dot_r  = pt1*r*cos(thetar);
    //REAL pt2_dot_u1 = pt2*u1*cos(NormalizeAngle(theta1 - phi));
    //REAL pt2_dot_u2 = pt2*u2*cos(NormalizeAngle(theta2 - phi));
    REAL pt2_dot_r  = pt2*r*cos(NormalizeAngle(thetar - phi));
    REAL u1_dot_r   = u1*r*cos(NormalizeAngle(theta1-thetar));
    REAL u2_dot_r   = u2*r*cos(NormalizeAngle(theta2-thetar));
    REAL u1_dot_u2  = u1*u2*cos(NormalizeAngle(theta2-theta1));

    // Norms
    // (u1-u2+r)^2
    REAL u1u2r_sqr = SQR(u1)+SQR(u2)+SQR(r) + 2.0*u1_dot_r - 2.0*u2_dot_r
        - 2.0*u1_dot_u2;
    REAL u1u2r = std::sqrt(u1u2r_sqr);
    // (u1+r)^2
    REAL u1r_sqr = SQR(u1)+SQR(r) + 2.0*u1_dot_r;
    REAL u1r = std::sqrt(u1r_sqr);
    // (u2-r)^2
    REAL u2r_sqr = SQR(u2)+SQR(r) - 2.0*u2_dot_r;
    REAL u2r = std::sqrt(u2r_sqr);

    // Contribution from e^(ik(u'-u)) e^(-i \Delta r),
    // k = pt1, \delta = pt1+pt2
    // Take only real part
    REAL result = cos(pt1_dot_u2 - pt1_dot_u1)*cos(pt1_dot_r + pt2_dot_r)
        + sin(pt1_dot_u2 - pt1_dot_u1)*sin(pt1_dot_r + pt2_dot_r);

    result *= N->S(u1u2r, y) *
        ( N->S(u1, y) * N->S(u2,y) - N->S(r, y) * N->S(u1u2r, y) );

    // Finally multiply by complicated F/F factor
    result *= ( Lx(u2r) - Lx(r) + Lx(u1r) - Lx(u1u2r) )
        / ( Lx(u1) - Lx(r) + Lx(u2) - Lx(u1u2r) );

    // Wave function product
    // \phi^*(u2) \phi(u) summed over spins and polarization
    // simple resul in massless z=0 case
    result *= -16.0*SQR(M_PI) * u1_dot_u2 / ( SQR(u1)*SQR(u2) );

    return result;


}

REAL NormalizeAngle(REAL phi)
{
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

REAL Lx(REAL x)
{
    return -SQR(x) * std::log( M_E + 1.0/(LAMBDAQCD2*SQR(x)) );
}

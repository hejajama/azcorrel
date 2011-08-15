/*
 * Class to calculate cross sections related to two-body correlations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "xs.hpp"
#include <gsl/gsl_integration.h>

extern "C"
{
    #include "amplitudelib/fourier/fourier.h"
}

/*
 * Return d\sigma / (dpt1 dpt2 dy1 dy2 d\theta) lo approximation
 */
REAL CrossSection2::dSigma_lo(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL theta, REAL sqrts)
{
    REAL tmpz = z(pt1, pt2, y1, y2);
    REAL tmpxa = xa(pt1, pt2, y1, y2, sqrts);

    REAL result = std::pow(1.0-tmpz, 3.0)*(1.0+SQR(1.0-tmpz));
    REAL qs = 1.0/N->SaturationScale(std::log(0.01/tmpxa), 0.5);
    result *= SQR(qs);

    REAL delta = Delta(pt1, pt2, theta);
    //result /= SQR(pt2) * SQR(delta)
    //    * ((1.0-2.0*tmpz)*SQR(pt1)+SQR(tmpz)*SQR(delta) - 2.0*tmpz*pt1*pt2*std::cos(theta));
    result /= SQR(pt2) * SQR(delta)
        * ( SQR(1.0-tmpz)*SQR(pt1) + SQR(tmpz*pt2) - 2.0*(1.0-tmpz)*tmpz*pt1*pt2*std::cos(theta) );

    REAL tmpxh = xh(pt1, pt2, y1, y2, sqrts);
    result *= 2.0*(pdf->xq(tmpxh, delta, UVAL) + pdf->xq(tmpxh, delta, DVAL));
        // factor 2 from isospin symmetry: xf_u,p = xf_d,n

    result *= pt1*pt2;

    result *= ALPHAS*Cf/(2.0*std::pow(M_PI,3));

    return result;
}

/*
 * Ref. 0708.0231 eq (57)&(54), no k_t factorization/other approximations
 */
REAL CrossSection2::dSigma(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL phi, REAL sqrts)
{
    //return dSigma_lo(pt1, pt2, y1,y2, phi, sqrts);
    REAL result=0;
    REAL tmpz = z(pt1, pt2, y1, y2);
    REAL tmpxa = xa(pt1, pt2, y1, y2, sqrts);
    REAL ya = std::log(0.01/tmpxa);
    REAL delta = Delta(pt1,pt2,phi);
    //cout << "# phi " << phi << " delta " << delta << endl;
    REAL g = G(pt2, tmpxa);
    //cout << "## g=" << g << endl;
    // Quite uqly...
    REAL f=0;
    if (std::abs(delta - fcachek) < 0.01) f=fcacheval;
    else
    {
        f = N->S_k(delta, ya);
        fcachek=delta; fcacheval=f;
    }
    //cout <<"## f=" << f << endl;
    

    // k - z\delta = (1-z)^2 pt1^2 + z^2 pt2^2 - 2*z*(1-z)*pt1*pt2*cos \phi
    REAL kzdeltasqr = SQR(1.0-tmpz)*SQR(pt1) + SQR(tmpz*pt2) - 2.0*tmpz*(1.0-tmpz)
                                * pt1*pt2*std::cos(phi);
    
    result = SQR(g) + 1.0/kzdeltasqr
        - 2.0*g*( (1.0-tmpz)*SQR(pt1) - tmpz*pt1*pt2*std::cos(phi) ) / ( pt1*kzdeltasqr );

    result *= 2.0*(1.0+SQR(1.0-tmpz) );
    result *= (1.0-tmpz);
    REAL tmpxh = xh(pt1, pt2, y1, y2, sqrts);
    result *= 2.0*(pdf->xq(tmpxh, delta, UVAL) + pdf->xq(tmpxh, delta, DVAL));
        // factor 2 from isospin symmetry: xf_u,p = xf_d,n

    result *= f;
    result *= pt1*pt2;

    result *= ALPHAS*Cf/(4.0*SQR(M_PI));
    //cout <<"## result " << result << endl;

    return result;
}



/*
 * Total cross section, eg. dSigma integrated over theta
 */
struct Inthelper_sigma
{
    CrossSection2 *xs;
    REAL pt1,pt2,y1,y2,sqrts;
};
REAL Inthelperf_sigma(REAL theta, void* p)
{
    Inthelper_sigma *par = (Inthelper_sigma*) p;
    REAL res1 = par->xs->dSigma(par->pt1, par->pt2, par->y1, par->y2, theta, par->sqrts);
    REAL res2 = par->xs->dSigma(par->pt2, par->pt1, par->y2, par->y1, theta, par->sqrts);

    if (isnan(res1) or isnan(res2) or isinf(res1) or isinf(res2))
    {
        cerr << "nan/inf at theta=" << theta << endl;
    }
    return res1+res2;
}

REAL CrossSection2::Sigma(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL sqrts)
{
    const REAL THETAINTPOINTS = 20;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);
    Inthelper_sigma par;
    par.pt1=pt1; par.pt2=pt2; par.y1=y1; par.y2=y2; par.sqrts=sqrts;
    par.xs=this;
    gsl_function fun;
    fun.params = &par;
    fun.function = Inthelperf_sigma;

    REAL mintheta = 0.0;
    REAL maxtheta = M_PI;

    int status; REAL result, abserr;
    status = gsl_integration_qag(&fun, mintheta,
            maxtheta, 0, 0.01, THETAINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl;

    return result*2.0;  // *2.0 as we integrated over [0,pi], not [0,2pi]

}

/*
 * Funktion G_{x_A} = \int dr (1-N(r)) J_1(k*r)
 * as in ref. 0708.0231 but w.o. vector k/|k|
 */
struct G_helper { REAL y; AmplitudeLib* N; REAL kt; };
REAL G_helperf(REAL r, void* p);
REAL CrossSection2::G(REAL kt, REAL x)
{
    if (std::abs(kt - gcachek) < 0.01)
        return gcacheval;
    
    G_helper helper;
    helper.y = std::log(0.01/x);
    helper.N=N; helper.kt=kt;

    set_fpu_state();
    init_workspace_fourier(500);   // number of bessel zeroes, max 2000

    REAL result = fourier_j1(kt, G_helperf, &helper);
    gcachek=kt; gcacheval=result;
    return result;

}

REAL G_helperf(REAL r, void *p)
{
    G_helper* par = (G_helper*) p;
    REAL result=0;
    if (r< par->N->MinR()) result = 1.0;
    else if (r > par->N->MaxR()) result=0;
    else result = 1.0 - par->N->N(r, par->y);

    return result;
}


CrossSection2::CrossSection2(AmplitudeLib* N_, PDF* pdf_,FragmentationFunction* frag)
{
    N=N_; pdf=pdf_; fragfun=frag;
    gcacheval=-1;
    gcachek=-1;
    fcacheval=-1; fcachek=-1;
}


REAL CrossSection2::Delta(REAL pt1, REAL pt2, REAL theta)
{
    return std::sqrt( SQR(pt1) + SQR(pt2) + 2.0*pt1*pt2*std::cos(theta) );
}

REAL CrossSection2::z(REAL pt1, REAL pt2, REAL y1, REAL y2)
{
    return std::abs(pt1)*std::exp(y1)
        / (std::abs(pt1)*std::exp(y1) + std::abs(pt2)*std::exp(y2) );
}

REAL CrossSection2::xa(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL sqrts)
{
    return ( std::abs(pt1)*std::exp(-y1) + std::abs(pt2) * std::exp(-y2) )
        / sqrts;
}

REAL CrossSection2::xh(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL sqrts)
{
    return ( std::abs(pt1)*std::exp(y1) + std::abs(pt2) * std::exp(y2) )
        / sqrts;
}

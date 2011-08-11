/*
 * Class to calculate cross sections related to two-body correlations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "xs.hpp"
#include <gsl/gsl_integration.h>



/*
 * Return d\sigma / (dpt1 dpt2 dy1 dy2 d\theta)
 */
REAL CrossSection2::dSigma(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL theta, REAL sqrts)
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
    result *= 2.0*(pdf->xq(tmpxh, delta, USEA) + pdf->xq(tmpxh, delta, DSEA));
        // factor 2 from isospin symmetry: xf_u,p = xf_d,n

    result *= pt1*pt2;

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
    return par->xs->dSigma(par->pt1, par->pt2, par->y1, par->y2, theta, par->sqrts)
        + par->xs->dSigma(par->pt2, par->pt1, par->y2, par->y1, theta, par->sqrts);
}

REAL CrossSection2::Sigma(REAL pt1, REAL pt2, REAL y1, REAL y2, REAL sqrts)
{
    const REAL THETAINTPOINTS = 200;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);
    Inthelper_sigma par;
    par.pt1=pt1; par.pt2=pt2; par.y1=y1; par.y2=y2; par.sqrts=sqrts;
    par.xs=this;
    gsl_function fun;
    fun.params = &par;
    fun.function = Inthelperf_sigma;

    REAL mintheta = 0.0;
    REAL maxtheta = 2.0*M_PI;

    int status; REAL result, abserr;
    status = gsl_integration_qag(&fun, mintheta,
            maxtheta, 0, 0.01, THETAINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl;

    return result;


}


CrossSection2::CrossSection2(AmplitudeLib* N_, PDF* pdf_)
{
    N=N_; pdf=pdf_;
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

/*
 * Class to calculate cross sections related to two-body correlations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "xs.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>

#include <gsl/gsl_monte.h>


extern "C"
{
    #include "amplitudelib/fourier/fourier.h"
}

/*
 * Return d\sigma / (dpt1 dpt2 dy1 dy2 d\theta) lo approximation
 *
 * if multiply_pdf=true (default) then multiply by parton distribution function
 */
double CrossSection2::dSigma_lo(double pt1, double pt2, double y1, double y2, double theta,
    double sqrts, bool multiply_pdf)
{
    double tmpz = z(pt1, pt2, y1, y2);
    double tmpxa = xa(pt1, pt2, y1, y2, sqrts);

    double result = std::pow(1.0-tmpz, 3.0)*(1.0+SQR(1.0-tmpz));
    double qs = 1.0/N->SaturationScale(std::log(0.01/tmpxa), 0.5);
    result *= SQR(qs);

    double delta = Delta(pt1, pt2, theta);
    //result /= SQR(pt2) * SQR(delta)
    //    * ((1.0-2.0*tmpz)*SQR(pt1)+SQR(tmpz)*SQR(delta) - 2.0*tmpz*pt1*pt2*std::cos(theta));
    result /= SQR(pt2) * SQR(delta)
        * ( SQR(1.0-tmpz)*SQR(pt1) + SQR(tmpz*pt2) - 2.0*(1.0-tmpz)*tmpz*pt1*pt2*std::cos(theta) );

    double tmpxh = xh(pt1, pt2, y1, y2, sqrts);

    if (multiply_pdf)
        result *= 2.0*(pdf->xq(tmpxh, delta, UVAL) + pdf->xq(tmpxh, delta, DVAL));
        // factor 2 from isospin symmetry: xf_u,p = xf_d,n

    result *= pt1*pt2;

    result *= ALPHAS*Cf/(2.0*std::pow(M_PI,3));

    return result;
}

/*
 * Ref. 0708.0231 eq (57)&(54), no k_t factorization/other approximations
 *
 * If pdf=true (default), then multiply by parton distribution function
 */
double CrossSection2::dSigma(double pt1, double pt2, double y1, double y2, double phi,
    double sqrts, bool multiply_pdf)
{
    //return dSigma_lo(pt1, pt2, y1,y2, phi, sqrts);
    double result=0;
    double tmpz = z(pt1, pt2, y1, y2);
    double tmpxa = xa(pt1, pt2, y1, y2, sqrts);
    double ya = std::log(0.01/tmpxa);
    //N->InitializeInterpolation(ya);
    double delta = Delta(pt1,pt2,phi);
    //cout << "# phi " << phi << " delta " << delta << endl;
    double g=0,f=0;
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            g = G(pt2, tmpxa);
        }
        #pragma omp section
        {
            if (std::abs(delta - fcachek) < 0.01) f=fcacheval;
            else
            {
                f = N->S_k(delta, ya);
                fcachek=delta; fcacheval=f;
            }
        }
    }

    

    // (k - z\delta)^2 = (1-z)^2 pt1^2 + z^2 pt2^2 - 2*z*(1-z)*pt1*pt2*cos \phi
    double kzdeltasqr = SQR(1.0-tmpz)*SQR(pt1) + SQR(tmpz*pt2) - 2.0*tmpz*(1.0-tmpz)
                                * pt1*pt2*std::cos(phi);
    
    result = SQR(g) + 1.0/kzdeltasqr
        - 2.0*g*( (1.0-tmpz)*SQR(pt1) - tmpz*pt1*pt2*std::cos(phi) ) / ( pt1*kzdeltasqr );

    result *= 2.0*(1.0+SQR(1.0-tmpz) );
    

    result *= f;
    
    // Add correction term following from more exact calculationg of the
    // 4-point function
    double correction = CorrectionTerm(pt1,pt2, ya, phi);
    cout << "# relcorrection at pt2=" << pt2 << " = " << std::abs(correction/result) << endl;
    result+=correction;

    double tmpxh = xh(pt1, pt2, y1, y2, sqrts);

    if (multiply_pdf)
        result *= 2.0*(pdf->xq(tmpxh, delta, UVAL) + pdf->xq(tmpxh, delta, DVAL));
        // factor 2 from isospin symmetry: xf_u,p = xf_d,n

    // Multiply terms which follow from change of variable d\sigma/d^3kd^q
    // => d\simga / d^2kd^2 q dy_1 dy_2
    // and integration over x (takes away the delta function)
    // meaning (1-z)*k^+
    // However k^+ has already cancelled out and is also cancelled out in
    // CorrectionTerm()
    result *= (1.0-tmpz);

    result *= pt1*pt2;  // d^2 pt_1 d^2 pt2 => dpt_1 dpt_2

    // Overall constants
    result *= ALPHAS*Cf/(4.0*SQR(M_PI));
    return result;
}



/*
 * Total cross section, eg. dSigma integrated over theta
 */
struct Inthelper_sigma
{
    CrossSection2 *xs;
    double pt1,pt2,y1,y2,sqrts;
};
double Inthelperf_sigma(double theta, void* p)
{
    Inthelper_sigma *par = (Inthelper_sigma*) p;
    double res1 = par->xs->dSigma(par->pt1, par->pt2, par->y1, par->y2, theta, par->sqrts);
    double res2 = par->xs->dSigma(par->pt2, par->pt1, par->y2, par->y1, theta, par->sqrts);

    if (isnan(res1) or isnan(res2) or isinf(res1) or isinf(res2))
    {
        cerr << "nan/inf at theta=" << theta << endl;
    }
    return (res1+res2);
}

double CrossSection2::Sigma(double pt1, double pt2, double y1, double y2, double sqrts)
{
    const double THETAINTPOINTS = 20;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);
    Inthelper_sigma par;
    par.pt1=pt1; par.pt2=pt2; par.y1=y1; par.y2=y2; par.sqrts=sqrts;
    par.xs=this;
    gsl_function fun;
    fun.params = &par;
    fun.function = Inthelperf_sigma;

    double mintheta = 0.0;
    double maxtheta = M_PI;

    int status; double result, abserr;
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
 * Coinsidence probability
 * Ref. 1005.4065, eq. (7)
 * We integrate over a finite range of final state rapidities and transverse
 * momenta. Also take into account fragmentation => integrate over z_i from x_1
 * up to 1
 */
struct NPairHelper{
    CrossSection2* xs;
    AmplitudeLib* N;
    double sqrts;
    double phi;
    double yh1,yh2;   // Rapidities of the hadrons
};

double NPairHelperf(double* vec, size_t dim, void* p);

double CrossSection2::NPair(double phi, double sqrts)
{
    double maxpt = 10;
    const uint dim=7;
    double lower[dim] = {0, 2,     1,     0, 0 ,2.4,2.4 };
    double upper[dim] = {1, maxpt, maxpt, 1, 1 ,4,4 };


    NPairHelper helper; helper.sqrts=sqrts; helper.phi=phi;
    helper.xs=this; helper.N=N;
    //helper.yh1=3.5; helper.yh2=2;

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup ();
    gsl_monte_function montef = { &NPairHelperf, dim, &helper };   

    size_t calls = 1000;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    double result,abserr;
    
    /*gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim);
    gsl_monte_miser_integrate(&montef, lower, upper, dim, calls, r, s,
        &result, &abserr);
    gsl_monte_miser_free(s);*/
    
    gsl_monte_plain_state *s = gsl_monte_plain_alloc (dim);
    gsl_monte_plain_integrate(&montef, lower, upper, dim, calls, r, s,
        &result, &abserr);
    gsl_monte_plain_free(s);
    

    gsl_rng_free(r);

    cout << "# phi " << phi << " result " << result << " relerror "
     << std::abs(abserr/result) << endl;

    return result;
}

double NPairHelperf(double* vec, size_t dim, void* p)
{
    NPairHelper* par = (NPairHelper*)p;
    // vec[0]=x, vec[1]=p_t1, vec[2]=p_t2, vec[3]=z1,
    // vec[4]=z2  vec[5]=y1  vec[6]=y2  (vec[5] and 6 may not be integrated over)

    // If rapidity is integrated over
    double yh1 = vec[5];
    double yh2 = vec[6];
    
    //double yh1=par->yh1; double yh2=par->yh2;

    // z_1,z_2 integral
    double minx1 = vec[3]/par->sqrts*std::exp(yh1);
    double minx2 = vec[4]/par->sqrts*std::exp(yh2);
    if (vec[3]<minx1) return 0;
    if (vec[4]<minx2) return 0;

    // x integral
    if (vec[0] < minx1/vec[3]+minx2/vec[4]) return 0;
    if (minx1/vec[3] + minx2/vec[4] >=1) return 0;

    // p_2 integral
    if (vec[2]>vec[1]) return 0;

    // y1 and y2 for q->qg process
    // p^+ = p_T e^y, approx. p_T constant during the fragmentation
    // y_{i, q->qg} = ln[ p_i^+ / (z|p_i|) ]
    // y_{i, q->qg} = y_i + ln 1/z_i
    double y1 = yh1 + std::log(1.0/vec[3]);
    double y2 = yh2 + std::log(1.0/vec[4]);

    double result=0;

    FragmentationFunction* fragfun = par->xs->FragFun();

    // pdf and fragmentation scales are assumed to be p_t1 (leading hadron)
    // fragmentation to neutral pions
    result =
        par->xs->dSigma(vec[1], vec[2], y1, y2, par->phi, par->sqrts, false)
          * fragfun->Evaluate(G, PI0, vec[4], vec[1])*
          ( fragfun->Evaluate(U, PI0, vec[3], vec[1])
           +fragfun->Evaluate(D, PI0, vec[3], vec[1]) )
        + par->xs->dSigma(vec[2], vec[1], y2, y1, par->phi, par->sqrts, false)
          * fragfun->Evaluate(G, PI0, vec[3], vec[1])*
          ( fragfun->Evaluate(U, PI0, vec[4], vec[1])
           +fragfun->Evaluate(D, PI0, vec[4], vec[1]) );

    // We can multiply the whole expression by this sum of PDF's as due to the
    // isospin symmetry we have the same combination for both u and d quarks
    result *= par->xs->Pdf()->xq(vec[0], vec[1], UVAL)
        + par->xs->Pdf()->xq(vec[0], vec[1], DVAL);
    
    return result;
}


/*
 * Funktion G_{x_A} = \int dr (1-N(r)) J_1(k*r)
 * as in ref. 0708.0231 but w.o. vector k/|k|
 */
struct G_helper { double y; AmplitudeLib* N; double kt; };
double G_helperf(double r, void* p);
double CrossSection2::G(double kt, double x)
{
    if (std::abs(kt - gcachek) < 0.01)
        return gcacheval;
    
    G_helper helper;
    helper.y = std::log(0.01/x);
    helper.N=N; helper.kt=kt;

    set_fpu_state();
    init_workspace_fourier(700);   // number of bessel zeroes, max 2000

    double result = fourier_j1(kt, G_helperf, &helper);
    gcachek=kt; gcacheval=result;
    return result;

}

double G_helperf(double r, void *p)
{
    G_helper* par = (G_helper*) p;
    double result=0;
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


double CrossSection2::Delta(double pt1, double pt2, double theta)
{
    return std::sqrt( SQR(pt1) + SQR(pt2) + 2.0*pt1*pt2*std::cos(theta) );
}

double CrossSection2::z(double pt1, double pt2, double y1, double y2)
{
    return std::abs(pt1)*std::exp(y1)
        / (std::abs(pt1)*std::exp(y1) + std::abs(pt2)*std::exp(y2) );
}

double CrossSection2::xa(double pt1, double pt2, double y1, double y2, double sqrts)
{
    return ( std::abs(pt1)*std::exp(-y1) + std::abs(pt2) * std::exp(-y2) )
        / sqrts;
}

double CrossSection2::xh(double pt1, double pt2, double y1, double y2, double sqrts)
{
    return ( std::abs(pt1)*std::exp(y1) + std::abs(pt2) * std::exp(y2) )
        / sqrts;
}

PDF* CrossSection2::Pdf()
{
    return pdf;
}

FragmentationFunction* CrossSection2::FragFun()
{
    return fragfun;
}

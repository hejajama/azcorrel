/*
 * Class to calculate cross sections related to two-body correlations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "xs.hpp"
#include <tools/interpolation2d.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <fstream>
#include <string>
#include <sstream>


extern "C"
{
    #include <fourier/fourier.h>
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

    //double result = std::pow(1.0-tmpz, 3.0)*(1.0+SQR(1.0-tmpz));
    double result = std::pow(1.0-tmpz, 2.0)*(1.0+SQR(1.0-tmpz));
    double qs = 1.0/N->SaturationScale(std::log(N->X0()/tmpxa), 0.5);
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

    result *= 2.0/M_PI; // 1/(2\pi)^2 is factorized out

    /*result *= pt1*pt2;

    result *= ALPHAS*Cf/(2.0*std::pow(M_PI,3));
    */

    return result;
}

/*
 * Ref. 0708.0231 eq (57)&(54), no k_t factorization/other approximations
 *
 * If pdf=true (default), then multiply by parton distribution function
 * Returns -1 if integral doesn't converge and MC integral is used
 */
double CrossSection2::dSigma(double pt1, double pt2, double y1, double y2, double phi,
    double sqrts, bool multiply_pdf)
{
    //return dSigma_lo(pt1, pt2, y1,y2, phi, sqrts);
    double result=0;
    double tmpz = z(pt1, pt2, y1, y2);
    double tmpxa = xa(pt1, pt2, y1, y2, sqrts);
    double ya = std::log(N->X0()/tmpxa);
    //N->InitializeInterpolation(ya);
    double delta = Delta(pt1,pt2,phi);

    double g=0,f=0;
    /*
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            g = G(pt2, tmpxa, tmpz);
        }
        #pragma omp section
        {
            if (std::abs(delta - fcachek) < 0.01) f=fcacheval;
            else
            {
                f = N->S_k(delta, ya)/SQR(2.0*M_PI);
                fcachek=delta; fcacheval=f;
            }
        }
    }
    */
    g = G(pt2, tmpxa, tmpz);
    f=N->S_k(delta, ya)/SQR(2.0*M_PI);


    // (k - z\delta)^2 = (1-z)^2 pt1^2 + z^2 pt2^2 - 2*z*(1-z)*pt1*pt2*cos \phi
    double kzdeltasqr = SQR(1.0-tmpz)*SQR(pt1) + SQR(tmpz*pt2) - 2.0*tmpz*(1.0-tmpz)
                                * pt1*pt2*std::cos(phi);

    //cout << "# phi " << phi << " delta " << delta << "kzdeltasqr " << kzdeltasqr <<   endl;
    // m=0
    if (M_Q() < 1e-5)
    {
        result = SQR(g) + 1.0/kzdeltasqr + 2.0*g*( (1.0-tmpz)*pt1*pt2*std::cos(phi)
                    - tmpz*SQR(pt2) ) / ( pt2*kzdeltasqr );
        result *= 2.0*(1.0+SQR(1.0-tmpz) );
    }
    else // m!= 0
        result = 2.0*(1.0+SQR(1.0-tmpz))*(
            SQR(g) + kzdeltasqr/SQR(kzdeltasqr + SQR(M_Q()*tmpz))
                + 2.0*g*( (1.0-tmpz)*pt1*pt2*std::cos(phi) - tmpz*SQR(pt2) )
                    / (pt2* ( kzdeltasqr + SQR(M_Q()*tmpz) ) )
            )
            + 2.0*SQR( H(pt2, tmpxa, tmpz) - M_Q()*SQR(tmpz)/(kzdeltasqr + SQR(M_Q()*tmpz)) );
    
    
    

    result *= f;
    //#pragma omp critical
    //cout << "# phi=" << phi << ", result w.o. corrections = " << result << endl;

    // CorrectionTerm returns -1 if the integral doesn't converge
    
    double correction = CorrectionTerm(pt1,pt2,ya,phi,tmpz);
    //if (std::abs(correction+1.0)<0.001) return -1;
   	//result +=correction;
    result = correction;
    
    /*result = CorrectionTerm_fft(pt1, pt2, ya, phi);
    #pragma omp critical
    cout <<"# phi=" << phi <<" MC result " << result << endl;
    */

    if (!multiply_pdf)
    {
        // return Marquet's k^+|foo|^2F() with correction term, not multiplied
        // by any constants, so ~d\sigma/d^3kd^3q
        // Marquet's M w.o. S/(2\pi)^2
        return result;
    }
    
    double tmpxh = xh(pt1, pt2, y1, y2, sqrts);

    result *= 2.0*(pdf->xq(tmpxh, delta, UVAL) + pdf->xq(tmpxh, delta, DVAL));
    // factor 2 from isospin symmetry: xf_u,p = xf_d,n


    // Multiply terms which follow from change of variable d\sigma/d^3kd^q
    // => d\simga / d^2kd^2 q dy_1 dy_2
    // and integration over x (takes away the delta function)
    // meaning (1-z)*k^+
    // However k^+ has already cancelled and is also cancelled in
    // CorrectionTerm()
    //result *= (1.0-tmpz);

    //result *= pt1*pt2;  // d^2 pt_1 d^2 pt2 => dpt_1 dpt_2

    // Overall constants
    //result *= ALPHAS*Cf/(4.0*SQR(M_PI));
    
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
        cerr << "nan/inf at theta=" << theta << ". " << LINEINFO << endl;
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
/*
 * Load \Delta \phi data from files pt1_#_pt2_#_y1_#_y2_#, interpolate in pt and phi
 * Integrate over kinematical variables z_1,z_2 to take into account
 * hadronization
 *
 * Returns 0 if no error occurred, otherwise <0
 */
int CrossSection2::LoadPtData(double y1, double y2)
{
   // Load datafiles 
   // Generate interpolators

    std::stringstream tmp;
    std::string y1str, y2str;
    tmp << y1; y1str = tmp.str();
    std::stringstream tmp2; 
    tmp2 << y2; y2str = tmp2.str();

/*
    ptvals.push_back(1.0); ptvals.push_back(1.5); ptvals.push_back(2.0);
    ptvals.push_back(2.5); ptvals.push_back(3.0); ptvals.push_back(3.5);
    ptvals.push_back(4.0); ptvals.push_back(4.5); ptvals.push_back(5.0);
   
    ptstrings.push_back("1"); ptstrings.push_back("1.5"); ptstrings.push_back("2");
    ptstrings.push_back("2.5"); ptstrings.push_back("3"); ptstrings.push_back("3.5");
    ptstrings.push_back("4"); ptstrings.push_back("4.5"); ptstrings.push_back("5");
*/
    for (uint i=0; i<ptinterpolators.size(); i++)
    {
        for (uint j=0; j<ptinterpolators[i].size(); j++)
        {
            delete ptinterpolators[i][j];
        }
        ptinterpolators[i].clear();
    }
    ptinterpolators.clear();

   int points=ptvals.size();
   for (int pt1ind=0; pt1ind<points; pt1ind++)
   {
       std::vector<Interpolator*> tmpinterpolators;
       for (int pt2ind=0; pt2ind<points; pt2ind++)
       {
            std::stringstream fname;
            fname << "taulukko_pt_y_pp/pt1_" << ptstrings[pt1ind] << "_pt2_"
                << ptstrings[pt2ind] << "_y1_" << y1str << "_y2_" << y2str;
            //cout << "# Loading file " << fname.str() << endl;
            std::ifstream file(fname.str().c_str());
            if (!file.is_open())
            {
                cerr << "Can't open file " << fname.str()
                    << " " << LINEINFO << endl;
                return -1;
            }
            std::vector<double> dphi;
            std::vector<double> xs;
            while(!file.eof())
            {
                string line;
                getline(file, line);
                if (line.substr(0,1)=="#" or line.length()<5)
                    continue;   // Comment
                std::istringstream iss(line);
                string angle,val;
                iss >> angle; iss >> val;

                dphi.push_back(StrToReal(angle));
                xs.push_back(StrToReal(val));
            }

            file.close();

            if (xs.size()<2)
            {
                cerr << "Datafile without valid datapoints, filename " << fname.str()
                    << " " << LINEINFO << endl;
                continue;
            }

            // Mirror \dphi if necessary
            if (dphi[dphi.size()-1] < 1.01*M_PI)
            {
                int center_index = dphi.size()-1;
                for (int i=dphi.size()-2; i>=0; i--)
                {
                    dphi.push_back(M_PI + (M_PI-dphi[i]) );
                    xs.push_back(xs[i]);
                }
            }

            

            Interpolator *tmpinterp = new Interpolator(dphi, xs);
            tmpinterp->Initialize();
            tmpinterpolators.push_back(tmpinterp);

       }
       ptinterpolators.push_back(tmpinterpolators);
   }
   //cout << "# Data loaded" << endl; 

    return 0;
}

struct dSigma_full_helper
{
    double pt1,pt2;
    double minpt2;
    Interpolator2D *pt_interpolator;
    CrossSection2 *xs;
    double x1,x2;
    double y1, y2;
    double z1;
    double phi;
    double xh;
    bool deuteron;
    double sqrts;
};

double dSigma_full_helperf_z1(double z1, void* p);
double dSigma_full_helperf_z2(double z2, void* p);

const int PTINT_INTERVALS=5;

/*
 * Calculate dN/d^2 p_1 d^2 p_2 dy_1 dy_2, integrated over z1 and z2
 */

double CrossSection2::dSigma_full(double pt1, double pt2, double y1, double y2,
    double phi, double sqrts, bool deuteron)
{
    // Prepare 2D grid to interpolate
    std::vector< std::vector< double > > data;   //[pt1][pt2]
    for (uint pt1ind = 0; pt1ind < ptvals.size(); pt1ind++)
    {
        std::vector<double> tmpvec;
        for (uint pt2ind=0; pt2ind < ptvals.size(); pt2ind++)
        {
            tmpvec.push_back(ptinterpolators[pt1ind][pt2ind]->Evaluate(phi));
        }
        data.push_back(tmpvec);
    }
    Interpolator2D interpolator(ptvals, data);
    double x1 = pt1*std::exp(y1)/sqrts;
    double x2 = pt2*std::exp(y2)/sqrts;

    dSigma_full_helper helper;
    helper.pt1=pt1; helper.pt2=pt2;
    helper.pt_interpolator = &interpolator;
    helper.xs=this;
    helper.x1=x1; helper.x2=x2; helper.phi=phi;
    helper.xh = x1+x2;
    helper.deuteron=deuteron;

    gsl_function f;
    f.function=dSigma_full_helperf_z1;
    f.params=&helper;


    double result,abserr;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(PTINT_INTERVALS);
    int status = gsl_integration_qag(&f, x1, 1.0,
            0, 0.01, PTINT_INTERVALS, GSL_INTEG_GAUSS51, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    if (status)
        cerr << "z1int failed at " << LINEINFO <<", result " << result
        << " relerror " << std::abs(abserr/result) << " dphi " << phi << endl;

    // if we calculate dN/d^2kd^2qdy_1dy_2 instead of dN/d^3kd^3q
    result *= (1.0 - z(pt1, pt2, y1, y2));

    result *= ALPHAS * Cf / (4.0*M_PI); // NB! \int d^2b = S_T is dropped as it
    // should cancel, doesn't work anymore if we calculated something b-dependent!!!!

    return result;
}

double dSigma_full_helperf_z1(double z1, void* p)
{
    dSigma_full_helper *par = (dSigma_full_helper*)p;
    par->z1=z1;
    gsl_function f; f.function = dSigma_full_helperf_z2;
    f.params = par;

    double result,abserr;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(PTINT_INTERVALS);
    int status = gsl_integration_qag(&f, par->x2,1.0,
            0, 0.01, PTINT_INTERVALS, GSL_INTEG_GAUSS51, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    /*if (status)
        cerr << "z2int failed at " << LINEINFO <<", result " << result
        << " relerror " << std::abs(abserr/result) << " dphi " << par->phi << endl;
    */return result;
}

double dSigma_full_helperf_z2(double z2, void* p)
{    
    dSigma_full_helper *par = (dSigma_full_helper*)p;
    
    if (par->x1/par->z1 + par->x2/z2 > 1) return 0;
    if (par->pt1/par->z1 >= par->xs->MaxPt() or par->pt2/z2 >= par->xs->MaxPt()) return 0;

    double result=0;
    double qscale = std::max(par->pt1, par->pt2);
    ///TODO: why this was multiplied by z1?
    result = par->pt_interpolator->Evaluate(par->pt1/par->z1, par->pt2/z2);
    double xf_frag1 = 
          par->xh*par->xs->Pdf()->xq(par->xh, qscale, U)
            * par->xs->FragFun()->Evaluate(U, PI0, par->z1, qscale)
            * par->xs->FragFun()->Evaluate(G, PI0, z2, qscale)
          + par->xs->Pdf()->xq(par->xh, qscale, D)
            *  par->xs->FragFun()->Evaluate(D, PI0, par->z1, qscale)
            *  par->xs->FragFun()->Evaluate(G, PI0, z2, qscale);

    bool deuteron = par->deuteron;
    // u in proton = d in deuteron
    if (deuteron)
        xf_frag1 +=
         par->xh*par->xs->Pdf()->xq(par->xh, qscale, D) // u in neutron
            * par->xs->FragFun()->Evaluate(U, PI0, par->z1, qscale)
            * par->xs->FragFun()->Evaluate(G, PI0, z2, qscale)
          + par->xs->Pdf()->xq(par->xh, qscale, U) // d in neutron
            *  par->xs->FragFun()->Evaluate(D, PI0, par->z1, qscale)
            *  par->xs->FragFun()->Evaluate(G, PI0, z2, qscale);
    result *= xf_frag1;

    // exchange pt1<->pt2
    ///TODO: why this was multiplied by z2?
    
    double xf_frag2 =
        par->xh*par->xs->Pdf()->xq(par->xh, qscale, U)
            * par->xs->FragFun()->Evaluate(U, PI0, z2, qscale)
            * par->xs->FragFun()->Evaluate(G, PI0, par->z1, qscale)
          + par->xs->Pdf()->xq(par->xh, qscale, D)
            *  par->xs->FragFun()->Evaluate(D, PI0, z2, qscale)
            *  par->xs->FragFun()->Evaluate(G, PI0, par->z1, qscale);
    if (deuteron)
        xf_frag2 +=
         par->xh*par->xs->Pdf()->xq(par->xh, qscale, D) // u in neutron
            * par->xs->FragFun()->Evaluate(U, PI0, z2, qscale)
            * par->xs->FragFun()->Evaluate(G, PI0, par->z1, qscale)
          + par->xs->Pdf()->xq(par->xh, qscale, U)  // d in neutron
            *  par->xs->FragFun()->Evaluate(D, PI0, z2, qscale)
            *  par->xs->FragFun()->Evaluate(G, PI0, par->z1, qscale);

    result += par->pt_interpolator->Evaluate(par->pt2/z2, par->pt1/par->z1)*xf_frag2;

   /*cerr << "pt1/z1 " << par->pt1/par->z1 << " pt2/z2 " << par->pt2/z2 << " amp1 "
    << par->pt_interpolator->Evaluate(par->pt1/par->z1, par->pt2/z2)
        << " amp2 " << par->pt_interpolator->Evaluate(par->pt2/z2, par->pt1/par->z1) << endl;

    */return result;
}

/*
 * Integrate dSigma_full = dN/d^2 p_1 d^2 p_2 dy_1 dy_2 over
 * pt_1, pt_2, y_1, y_2 and one angle => get dN/d\phi
 */
double Inthelperf_cp_pt1(double pt1, void* p);
double Inthelperf_cp_pt2(double pt2, void* p);
const int PTINT_INTERVALS_CP=1;
const double PTINT_ACCURACY_CP=0.01;

double CrossSection2::dSigma_integrated(double minpt1, double minpt2, double miny, double maxy,
            double phi, double sqrts, bool deuteron)
{
    // Assume slow y-dependence at first....
    double y = 0.5*(miny+maxy);

    dSigma_full_helper helper;
    helper.minpt2=minpt2;
    helper.xs=this;
    helper.phi=phi;
    helper.deuteron=deuteron;
    helper.y1=y; helper.y2=y;
    helper.sqrts=sqrts;

    std::vector<double> yvals;
    yvals.push_back(miny);
    yvals.push_back(3.2);
    yvals.push_back(maxy);

    gsl_function fun;
    fun.function = Inthelperf_cp_pt1;
    fun.params=&helper;

    double result=0;

    for (int y1ind=0; y1ind<3; y1ind++)
    {
        double y2intres=0;
        for (int y2ind=0; y2ind<3; y2ind++)
        {
            cout << "# " << phi << " y1 " << yvals[y1ind] << " y2 " << yvals[y2ind] << endl;
            helper.y1=yvals[y1ind];
            helper.y2=yvals[y2ind];
            LoadPtData(helper.y1, helper.y2);
            double intresult,abserr;
            
            gsl_integration_workspace *workspace 
             = gsl_integration_workspace_alloc(PTINT_INTERVALS_CP);
            int status = gsl_integration_qag(&fun, minpt1, minpt1*2.0,
                    0, PTINT_ACCURACY_CP, PTINT_INTERVALS_CP, GSL_INTEG_GAUSS15, workspace,
                    &intresult, &abserr);
            gsl_integration_workspace_free(workspace);

            if (status)
            {
                cerr << "pt1 int failed at " << LINEINFO << " intresult: " << intresult
                    << " relerr " << std::abs(abserr/result) << " phi: " << phi << endl;
            }
            if (y2ind==1) y2intres += 4.0*intresult;
            else y2intres += intresult;
        }
        y2intres *= (maxy-miny)/6.0;
        if (y1ind==1) result += 4.0*y2intres;
        else result += y2intres;
    }
    result *= (maxy-miny)/6.0;
    
    //result *= (maxy-miny)*(maxy-miny);  // Assume y-indep.
    
    return result*2.0*M_PI; //2pi from one angular integral
}

double Inthelperf_cp_pt1(double pt1, void* p)
{
    dSigma_full_helper* par = (dSigma_full_helper*)p;
    par->pt1=pt1;

    gsl_function fun;
    fun.function = Inthelperf_cp_pt2;
    fun.params=par;

    double result,abserr;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(PTINT_INTERVALS_CP);
    int status = gsl_integration_qag(&fun, par->minpt2, pt1,
            0, PTINT_ACCURACY_CP, PTINT_INTERVALS_CP, GSL_INTEG_GAUSS15, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        cerr << "pt2 int failed at " << LINEINFO << " result: " << result
            << " relerr " << std::abs(abserr/result) << " pt1: " << pt1 << endl;
    }

    return result;
}

double Inthelperf_cp_pt2(double pt2, void* p)
{
    dSigma_full_helper* par = (dSigma_full_helper*)p;

    return par->pt1*pt2*par->xs->dSigma_full(par->pt1, pt2, par->y1, par->y2, par->phi, par->sqrts,
                        par->deuteron);
}

/*
 * Funktion G_{x_A} = \int dr S(r) J_1(k*r) (m=0)
 * Or with nonzero m,z
 * mz \int dr r S(r) K_1(mzr) J1(kr)
 * as in ref. 0708.0231 but w.o. vector k/|k|
 * Default value of z=0
 */
struct G_helper { double y; AmplitudeLib* N; double kt; double z; CrossSection2 *xs; };
double G_helperf(double r, void* p);
double CrossSection2::G(double kt, double x, double z)
{
    if (std::abs(kt - gcachek) < 0.001)
        return gcacheval;
    
    G_helper helper;
    helper.y = std::log(N->X0()/x);
    helper.N=N; helper.kt=kt; helper.z=z; helper.xs=this;

    set_fpu_state();
    init_workspace_fourier(700);   // number of bessel zeroes, max 2000

    double result = fourier_j1(kt, G_helperf, &helper);
    gcachek=kt; gcacheval=result;
    return result;

}

double G_helperf(double r, void *p)
{
    G_helper* par = (G_helper*) p;
    double M_Q = par->xs->M_Q();
    // Massless case
    if (M_Q<1e-5 or par->z<1e-5)
    {
        //cerr<<"Calculating massless case, are you sure? " << LINEINFO << endl;
        double result=0;
        if (r< par->N->MinR()) result = 1.0;
        else if (r > par->N->MaxR()) result=0;
        else result = 1.0 - par->N->N(r, par->y);

        return result;
    }

    // m,z>0
    double result=0;
    if (r< par->N->MinR()) result = r*gsl_sf_bessel_K1(r*M_Q*par->z);
    else if (r>par->N->MaxR()) result=0;
    else result = r*gsl_sf_bessel_K1(r*M_Q*par->z)*par->N->S(r, par->y);

    return M_Q*par->z*result;
}


/*
 * Funktion H_{x_A} = m*z^2*\int dr r S(r) K_0(mzr)J_0(kr)
 * as in ref. 0708.0231 
 */
double H_helperf(double r, void* p);
double CrossSection2::H(double kt, double x, double z)
{
    //TODO: Cache?
    if (M_Q() < 1e-5) return 0;
    G_helper helper;
    helper.y = std::log(N->X0()/x);
    helper.N=N; helper.kt=kt; helper.z=z; helper.xs=this;

    set_fpu_state();
    init_workspace_fourier(700);   // number of bessel zeroes, max 2000

    double result = fourier_j0(kt, H_helperf, &helper);
    return M_Q()*SQR(z)*result;

}

double H_helperf(double r, void *p)
{
    G_helper* par = (G_helper*) p;
    // Massless case
    if (par->xs->M_Q()<1e-5 or par->z<1e-5)
    {
        cerr<<"Calculating massless case, are you sure? " << LINEINFO << endl;

        return 0;
    }

    // m,z>0
    double M_Q=par->xs->M_Q();
    double result=0;
    if (r< par->N->MinR()) result = r*gsl_sf_bessel_K0(r*M_Q*par->z);
    else if (r>par->N->MaxR()) result=0;
    else result = r*gsl_sf_bessel_K0(r*M_Q*par->z)*par->N->S(r, par->y);

    return result;
}

CrossSection2::CrossSection2(AmplitudeLib* N_, PDF* pdf_,FragmentationFunction* frag)
{
    N=N_; pdf=pdf_; fragfun=frag;
    gcacheval=-1;
    gcachek=-1;
    fcacheval=-1; fcachek=-1;
    transform=NULL;
    mcintpoints=1e7;
    gsl_rng_env_setup();
    m_q=0.14;

    for (double pt = 1; pt<=9.5; pt+=0.5)
    {
        ptvals.push_back(pt);
        std::stringstream s; s << pt;
        ptstrings.push_back(s.str());
    }
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

CrossSection2::~CrossSection2()
{
    if (transform!=NULL)
        delete[] transform;
    for (uint i=0; i<ptinterpolators.size(); i++)
    {
        for (uint j=0; j<ptinterpolators[i].size(); j++)
        {
            delete ptinterpolators[i][j];
        }
    }
}

void CrossSection2::SetM_Q(double mq_)
{
    m_q=mq_;
}

double CrossSection2::M_Q()
{
    return m_q;
}

double CrossSection2::MaxPt()
{
    return ptvals[ptvals.size()-1];
}

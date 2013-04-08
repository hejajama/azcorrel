/*
 * Class to calculate cross sections related to two-body correlations
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2012
 */

#include "xs.hpp"
#include "config.hpp"
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
#include <algorithm>

#ifdef USE_MPI
    #include <mpi.h>
#endif


extern "C"
{
    #include <fourier/fourier.h>
}

const bool STAR=true;	// false: phenix, true: STAR kinematics

/*
 * Return d\sigma / (dpt1 dpt2 dy1 dy2 d\theta) lo approximation
 *
 * if multiply_pdf=true (default) then multiply by parton distribution function
 */
double CrossSection2::dSigma_lo(double pt1, double pt2, double y1, double y2, double theta,
    double sqrts, bool multiply_pdf)
{
    double tmpz = Z(pt1, pt2, y1, y2);
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
 * Returns -99999999 if integral doesn't converge and MC integral is used
 *
 * pt1=k: gluon
 * pt2=q: quark
 */
double CrossSection2::dSigma(double pt1, double pt2, double y1, double y2, double phi,
    double sqrts)
{
	#ifdef USE_MPI
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	#endif
    //return dSigma_lo(pt1, pt2, y1,y2, phi, sqrts);
    double result=0;
    double tmpz = Z(pt1, pt2, y1, y2);
    double tmpxa = xa(pt1, pt2, y1, y2, sqrts);
    if (xh(pt1, pt2, y1, y2, sqrts)>1)	// Kinematically forbidden, dont waste time
	{
		return 0;
	}
    double ya = std::log(N->X0()/tmpxa);
    if (ya<0)
        cerr << "Negative rapidity interval at " << LINEINFO << endl;
    
    if (!N->InterpolatorInitialized(ya))
        N->InitializeInterpolation(ya);
    double delta = Delta(pt1,pt2,phi);

    double g=0,f=0,h=0;
 
    f=N->S_k(delta, ya)/SQR(2.0*M_PI);

    if (f<0)
    {
        cerr << "f=" << f << " at phi=" << phi << endl;
    }
    g = G(pt2, tmpxa, tmpz); h=H(pt2, tmpxa, tmpz);

    // \kappa^2 = (k - z\delta)^2 = (1-z)^2 pt1^2 + z^2 pt2^2 - 2*z*(1-z)*pt1*pt2*cos \phi
    double kzdeltasqr = SQR(1.0-tmpz)*SQR(pt1) + SQR(tmpz*pt2) - 2.0*tmpz*(1.0-tmpz)
                                * pt1*pt2*std::cos(phi);
    double kappa = std::sqrt(kzdeltasqr);

    //cout << "# phi " << phi << " delta " << delta << "kzdeltasqr " << kzdeltasqr <<   endl;
    // m=0
    if (M_Q() < 1e-5 )
    {
        result = SQR(g) + 1.0/kzdeltasqr + 2.0*g*( (1.0-tmpz)*pt1*pt2*std::cos(phi)
                    - tmpz*SQR(pt2) ) / ( pt2*kzdeltasqr );
        result *= 2.0*(1.0+SQR(1.0-tmpz) );
    }
    else // m!= 0
    {
        result = (1.0+SQR(1.0-tmpz))*SQR(g) + SQR(h)
            + 1.0/SQR( kzdeltasqr + SQR(M_Q()*tmpz) )
                * ( kzdeltasqr * (1.0+SQR(1.0-tmpz)) + SQR(M_Q()*tmpz*tmpz) )
            - 1.0/(kzdeltasqr + SQR(M_Q()*tmpz)) *
                ( 2.0*(1.0+SQR(1.0-tmpz))*(tmpz*SQR(pt2)-(1.0-tmpz)*pt1*pt2*std::cos(phi))*g/pt2
                 + 2.0*M_Q()*SQR(tmpz)*h );
        result *= 2.0;
    }

    result *= f;
    
    //#pragma omp critical
    //cout << "# phi=" << phi << ", result w.o. corrections = " << result << endl;
    
        double nc_suppress=0;		// pt1  pt2  pt1_dot_pt2
    if (FiniteNc())
    {
		nc_suppress += 2.0*WavefSqr_k(pt1, kappa, (1.0-tmpz)*SQR(pt1) - tmpz*pt1*pt2*std::cos(phi), tmpz); // 2\psi(k)\psi(\kappa)^*
		nc_suppress -= WavefSqr_k(kappa, kappa, kappa*kappa, tmpz);	// -\psi(\kappa)\psi^*(\kappa)
		//nc_suppress -= WavefSqr_k(pt1, pt1, pt1*pt1, tmpz);			// -\psi(k)\psi^*(k), numerically as contains DPS contribution
		nc_suppress *= f;
		
		// In finite-Nc case we subtract a bit less than the Marque's result,
		// compensate it here (subtract, as we add Marquet's result later 
		// = add too much.
		if (M_Q()<1e-5)
			nc_suppress -=  2.0*(1.0+SQR(1.0-tmpz) ) * SQR(g)  * f;
		else
			nc_suppress -=  2.0*((1.0+SQR(1.0-tmpz))*SQR(g) + SQR(h) ) * f;
		nc_suppress *= 1.0/SQR(Nc);
	}
	//nc_suppress=0;///DEBUG

    // CorrectionTerm returns -1 if the integral doesn't converge
    double correction=0;
    #ifdef USE_MPI
		if (id==0)
		    cout << "# uncorrected result " << result << " nc-suppressed " << nc_suppress << endl;
	#endif
	if (channel==Q_QG)
	    correction = CorrectionTerm(pt1,pt2,ya,phi,tmpz);
	else if (channel==G_QQ)
	{
		result = GluonQuarkQuark(pt1,pt2,y1,y2,phi,sqrts);
	} 
	else if (channel==G_GG)
	{
		result = GluonGluonGluon(pt1, pt2, y1, y2, phi, sqrts);
		return result;

	}
    if (correction < -9999)
    {
		// Integral didn't converge
		if (mcintpoints > mcintpoints_orig * 3) // If # of intpoints is already increased twise
		{
			#ifdef USE_MPI
			if (id==0)
			#endif
			{
			cout <<"# ERROR! Integral didn't converge even though we increased the number of mcintpoints. Quitting..." << endl;
			}
			return -99999999;	
		}
		mcintpoints*=3;
		#ifdef USE_MPI
		if (id==0)
		#endif
		{
			
			cout <<"# Integral didn't converge, increasing num. of intpoints by factor 3 to " << mcintpoints << endl;
		}
		return dSigma(pt1, pt2, y1, y2, phi, sqrts);
	}
    // Nc-supprsesed terms analytically

    correction +=nc_suppress;  
    if (channel != Q_QG) cerr << "Shoudn't be here...\n";	///DEBUG
	else result += correction;

    // return Marquet's k^+|foo|^2F() with correction term, not multiplied
    // by any constants, so ~d\sigma/d^3kd^3q
    // Marquet's M w.o. S/(2\pi)^2
    return result;
   
	// The following terms are taken into account when calculating
	// convolution with pdf and FF

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
    
    //return result;
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
 * Use saved ptinterpolators and ptinterpolators_rev (for z and 1-z)
 * to initialize 2D interpolators ptinterpolator2d and ptinterpolator2d_rev
 * respectively which return amplitude at given phi as a function of pt1 and
 * pt2
 */
void CrossSection2::Prepare2DInterpolators(double phi)
{
    if (ptinterpolator2d != NULL)
        delete ptinterpolator2d;
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
    ptinterpolator2d = new Interpolator2D(ptvals, data);


    if (ptinterpolator2d_rev != NULL)
        delete ptinterpolator2d_rev;
    data.clear();   //[pt1][pt2]
    for (uint pt1ind = 0; pt1ind < ptvals.size(); pt1ind++)
    {
        std::vector<double> tmpvec;
        for (uint pt2ind=0; pt2ind < ptvals.size(); pt2ind++)
        {
            tmpvec.push_back(ptinterpolators_rev[pt1ind][pt2ind]->Evaluate(phi));
        }
        data.push_back(tmpvec);
    }
    ptinterpolator2d_rev = new Interpolator2D(ptvals, data);

    if (!apply_corrections)
        return;
    
    if (ptinterpolator2d_correction != NULL)
        delete ptinterpolator2d_correction;
    data.clear();   //[pt1][pt2]
    for (uint pt1ind = 0; pt1ind < ptvals.size(); pt1ind++)
    {    
        std::vector<double> tmpvec;
        for (uint pt2ind=0; pt2ind < ptvals.size(); pt2ind++)
        {
            tmpvec.push_back(ptinterpolators_correction[pt1ind][pt2ind]->Evaluate(phi));
        }
        data.push_back(tmpvec);
    }
    ptinterpolator2d_correction = new Interpolator2D(ptvals, data);


    if (ptinterpolator2d_rev_correction != NULL)
        delete ptinterpolator2d_rev_correction;
    data.clear();   //[pt1][pt2]
    for (uint pt1ind = 0; pt1ind < ptvals.size(); pt1ind++)
    {
        std::vector<double> tmpvec;
        for (uint pt2ind=0; pt2ind < ptvals.size(); pt2ind++)
        {
            tmpvec.push_back(ptinterpolators_rev_correction[pt1ind][pt2ind]->Evaluate(phi));
        }
        data.push_back(tmpvec);
    }
    ptinterpolator2d_rev_correction = new Interpolator2D(ptvals, data);
}

struct dSigma_full_helper
{
    double pt1,pt2;
    double minpt1, minpt2;
    double miny,maxy;
    Interpolator2D *pt_interpolator, *pt_interpolator_rev;
    Interpolator2D *pt_interpolator_cor, *pt_interpolator_rev_cor;
    CrossSection2 *xs;
    double x1,x2;
    double y1, y2;
    double z1;
    double phi;
    bool deuteron;
    double sqrts;
    Channel channel;
};

double dSigma_full_helperf_z1(double z1, void* p);
double dSigma_full_helperf_z2(double z2, void* p);

const int PTINT_INTERVALS=3;
const int ZINT_INTERVALS=3;

/*
 * Calculate dN/d^2 p_1 d^2 p_2 dy_1 dy_2, integrated over z1 and z2
 */

double CrossSection2::dSigma_full(double pt1, double pt2, double y1, double y2,
    double phi, double sqrts, bool deuteron)
{
    double x1 = pt1*std::exp(y1)/sqrts;
    double x2 = pt2*std::exp(y2)/sqrts;

    if (x1>1 or x2>1)
        return 0;     // kinematically forbidden
    
    dSigma_full_helper helper;
    helper.pt1=pt1; helper.pt2=pt2;
    helper.pt_interpolator = ptinterpolator2d;
    helper.pt_interpolator_rev = ptinterpolator2d_rev;
    helper.pt_interpolator_cor = ptinterpolator2d_correction;
    helper.pt_interpolator_rev_cor = ptinterpolator2d_rev_correction;
    helper.channel=channel;
    helper.xs=this;
    helper.x1=x1; helper.x2=x2; helper.phi=phi;
    helper.y1=y1; helper.y2=y2; helper.sqrts=sqrts;
    
    helper.deuteron=deuteron;

    gsl_function f;
    f.function=dSigma_full_helperf_z1;
    f.params=&helper;


    double result,abserr;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(ZINT_INTERVALS);
    int status = gsl_integration_qag(&f, x1, 1.0,
            0, 0.01, ZINT_INTERVALS, GSL_INTEG_GAUSS21, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    /*if (status)
        cerr << "z1int failed at " << LINEINFO <<", result " << result
        << " relerror " << std::abs(abserr/result) << " dphi " << phi << endl;
    */

	double alphas=Alpha_s(SQR(std::max(pt1, pt2))); 
	//cout << "maxpt " << std::max(pt1,pt2) << " alphas " << alphas << endl;
	//alphas=0.2;
    result *= alphas / (4.0*SQR(M_PI)); // NB! \int d^2b = S_T is dropped as it
    // should cancel, wouldn't work anymore if we calculated something b-dependent!!!!
    
    if (channel==Q_QG)
		result *= Nc/2.0;
	if (channel==G_GG)
		result *= 8.0 * Nc/2.0;
    // Quark channel color factor is Cf*NC/(2Cf) = Nc/2
    // Nc/2 as we have Cf * Nc/(2Cf) = Nc/2 in large and finite-nc
    

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
     = gsl_integration_workspace_alloc(ZINT_INTERVALS);
    double minz2 = par->x2*z1/(z1+0.00001-par->x1); //minz1=x1, add small eps to force finiteness
    if (minz2>1) return 0;
    /*int status = gsl_integration_qag(&f, par->x2,1.0,
            0, 0.01, ZINT_INTERVALS, GSL_INTEG_GAUSS21, workspace,
            &result, &abserr);*/
    int status = gsl_integration_qag(&f, minz2,1.0,
            0, 0.01, ZINT_INTERVALS, GSL_INTEG_GAUSS21, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    /*if (status)
        cerr << "z2int failed at " << LINEINFO <<", result " << result
        << " relerror " << std::abs(abserr/result) << " dphi " << par->phi << endl;
        */
    return result/SQR(z1);
}

double dSigma_full_helperf_z2(double z2, void* p)
{    
	
	
    dSigma_full_helper *par = (dSigma_full_helper*)p;
    
    if (par->x1/par->z1 + par->x2/z2 > 1)
        return 0; 
    if (par->pt1/par->z1 >= par->xs->MaxPt() or par->pt2/z2 >= par->xs->MaxPt())
        return 0;


    double result=0;
    double qscale = std::max(par->pt1, par->pt2);
    bool deuteron = par->deuteron;
    Channel channel = par->channel;
    if (channel != Q_QG and deuteron)
		cerr<< "Gluon channel does not support deuteron...." << endl;

    double xh = par->xs->xh(par->pt1/par->z1, par->pt2/z2, par->y1, par->y2, par->sqrts);
    /*double xa = par->xs->xa(par->pt1/par->z1, par->pt2/z2, par->y1, par->y2, par->sqrts);
    par->xs->GetN()->InitializeInterpolation(std::log(par->xs->GetN()->X0()/xa));
    */

    double xf_u = par->xs->Pdf()->xq(xh, qscale, U);
    double xf_d = par->xs->Pdf()->xq(xh, qscale, D);
    double xf_g = par->xs->Pdf()->xq(xh, qscale, G);

    double frag_u_pi0_z1 = par->xs->FragFun()->Evaluate(U, PI0, par->z1, qscale);
    double frag_d_pi0_z1 = par->xs->FragFun()->Evaluate(D, PI0, par->z1, qscale);
    double frag_s_pi0_z1 = par->xs->FragFun()->Evaluate(S, PI0, par->z1, qscale);
    double frag_g_pi0_z1 = par->xs->FragFun()->Evaluate(G, PI0, par->z1, qscale);
    double frag_u_pi0_z2 = par->xs->FragFun()->Evaluate(U, PI0, z2, qscale);
    double frag_d_pi0_z2 = par->xs->FragFun()->Evaluate(D, PI0, z2, qscale);
    double frag_s_pi0_z2 = par->xs->FragFun()->Evaluate(S, PI0, z2, qscale);
    double frag_g_pi0_z2 = par->xs->FragFun()->Evaluate(G, PI0, z2, qscale);

    // Interpolators give parton level cross section which is not multiplied by (1-z)
    // Evaluate(gluon momentum, quark momentum)
    double zfac =  1.0 - par->xs->Z(par->pt1/par->z1, par->pt2/z2, par->y1, par->y2);
    // in g->qq channel this factor is z(1-z) instead of (1-z)
    if (channel==G_QQ or channel==G_GG) zfac *= par->xs->Z(par->pt1/par->z1, par->pt2/z2, par->y1, par->y2);
    
    result =
        (par->pt_interpolator->Evaluate(par->pt1/par->z1, par->pt2/z2)  )
        * zfac;
    
		
      

    double xf_frag1 = 
          xf_u
            * frag_u_pi0_z2
            * frag_g_pi0_z1
          + xf_d
            *  frag_d_pi0_z2
            *  frag_g_pi0_z1;

	// g->uu,dd,ss and quarks fragment, assume that u->pi0 and ubar->pi0 are same
	if (channel==G_QQ)
		xf_frag1 = xf_g *
			( frag_u_pi0_z1 * frag_u_pi0_z2 + frag_d_pi0_z1 *frag_d_pi0_z2 
				+ frag_s_pi0_z1 * frag_s_pi0_z2 );
	else if (channel==G_GG)
	{
		xf_frag1 = xf_g * frag_g_pi0_z1 * frag_g_pi0_z2;
	}
    
    // u in proton = d in deuteron
    if (deuteron)
    {
        xf_frag1 +=
         xf_d // u in neutron
            * frag_u_pi0_z2
            * frag_g_pi0_z1
          + xf_u // d in neutron
            *  frag_d_pi0_z2
            *  frag_g_pi0_z1;
    }
    result *= xf_frag1;


    // exchange pt1<->pt2

    double xf_frag2 =
        xf_u
            * frag_u_pi0_z1
            * frag_g_pi0_z2
          + xf_d
            *  frag_d_pi0_z1
            *  frag_g_pi0_z2;
    if (channel==G_QQ)
		xf_frag2 = xf_g *
			( frag_u_pi0_z1 * frag_u_pi0_z2 + frag_d_pi0_z1 *frag_d_pi0_z2 
				+ frag_s_pi0_z1 * frag_s_pi0_z2 );
	else if (channel == G_GG)
		xf_frag2=0;	// in g->gg channel both situations are already taken into account
		
    if (deuteron)
    {
        xf_frag2 +=
         xf_d // u in neutron
            * frag_u_pi0_z1
            * frag_g_pi0_z2
          + xf_u  // d in neutron
            *  frag_d_pi0_z1
            *  frag_g_pi0_z2;
    }
    
    zfac = (1.0 - par->xs->Z(par->pt2/z2, par->pt1/par->z1, par->y2, par->y1));
    if (channel==G_QQ) zfac *=   par->xs->Z(par->pt2/z2, par->pt1/par->z1, par->y2, par->y1);
    
    result += (par->pt_interpolator_rev->Evaluate(par->pt2/z2, par->pt1/par->z1)  )
        *  zfac
        * xf_frag2;
    
    
    return result/SQR(z2);
}

/*
 * Integrate dSigma_full = dN/d^2 p_1 d^2 p_2 dy_1 dy_2 over
 * pt_1, pt_2, y_1, y_2 and one angle => get dN/d\phi
 * 
 * TODO: pt and y ranges are hardcoded, parameters are ignored!!
 */
double Inthelperf_cp_pt1(double pt1, void* p);
double Inthelperf_cp_pt2(double pt2, void* p);
double Inthelperf_y1(double pt2, void* p);
double Inthelperf_y2(double pt2, void* p);
const int PTINT_INTERVALS_CP=1;
const int YINT_INTERVALS=1;
const double PTINT_ACCURACY_CP=0.01;

double CrossSection2::dSigma_integrated(double minpt1, double minpt2, double miny, double maxy,
            double phi, double sqrts, bool deuteron)
{

    dSigma_full_helper helper;
    helper.minpt2=minpt2;
    helper.minpt1=minpt1;
    helper.xs=this;
    helper.phi=phi;
    helper.deuteron=deuteron;
    helper.sqrts=sqrts;

    std::vector<double> yvals;
    if (STAR)
    {
		yvals.push_back(2.4);
		yvals.push_back(2.8);
		yvals.push_back(3.2);
		yvals.push_back(3.6);
		yvals.push_back(4);
	}
	else
	{
		yvals.push_back(3);
		//yvals.push_back(3.2);
		yvals.push_back(3.4);
		//yvals.push_back(3.6);
		yvals.push_back(3.8);
	}
    
    
    miny=yvals[0];
    maxy=yvals[yvals.size()-1];
    helper.miny=miny; helper.maxy=maxy;
    

    double maxpt = MaxPt();
    
    
    gsl_function fun;
    fun.function = Inthelperf_cp_pt1;
    fun.params=&helper;

    double result=0;
    
    for (uint y1ind=0; y1ind<yvals.size(); y1ind++)
    {
        double y2intres=0;
        for (uint y2ind=0; y2ind<yvals.size(); y2ind++)
        {
            double y1 = yvals[y1ind]; double y2=yvals[y2ind];
            
            helper.y1=y1;
            helper.y2=y2;
            

            // Load data and prepare interpolation in (pt1,pt2) plane
            LoadPtData(helper.y1, helper.y2);
            Prepare2DInterpolators(phi);
            double intresult,abserr;
            double minpt=0, maxpt=0;
            if (STAR)
            {
				minpt=2; maxpt=4;
			} else {
				minpt=1.6; maxpt=2;
			}
                       
            gsl_integration_workspace *workspace 
             = gsl_integration_workspace_alloc(PTINT_INTERVALS_CP);
            int status = gsl_integration_qag(&fun, minpt, maxpt,
                    0, PTINT_ACCURACY_CP, PTINT_INTERVALS_CP, GSL_INTEG_GAUSS15, workspace,
                    &intresult, &abserr);
            gsl_integration_workspace_free(workspace);
            cout << "# " << phi << " y1 " << y1 << " y2 " << y2 << " pt1range " << minpt << " - " << maxpt << " res " << intresult << "\n";
            cout.flush();
            if (status)
            {
                cerr << "pt1 int failed at " << LINEINFO << " intresult: " << intresult
                    << " relerr " << std::abs(abserr/intresult) << " phi: " << phi << "\n";
            }

			if (yvals.size()==3)
			{
				if (y2ind==1) y2intres += 4.0*intresult;
				else y2intres += intresult;
			}
			else if (yvals.size()==4)
			{
				if (y2ind==1 or y2ind==2) y2intres += 3.0*intresult;
				else y2intres += intresult;
			}
			else if (yvals.size()==5)
			{
				if (y2ind==1 or y2ind == 3) y2intres += 4.0*intresult;
				else if (y2ind==2) y2intres += 2.0*intresult;
				else y2intres += intresult;
			}
			else if (yvals.size()==1)
				y2intres += intresult;
			else
				cerr << "Unknown amount of yvals at "<< LINEINFO << endl;
        }
        if (yvals.size()==3)
        {
			y2intres *= (maxy-miny)/6.0;
			if (y1ind==1) result += 4.0*y2intres;
			else result += y2intres;
        }
        else if (yvals.size()==4)
        {
			y2intres *= (maxy-miny)/8.0;
			if (y1ind==1 or y1ind==2) result += 3.0*y2intres;
			else result += y2intres;
		}
        else if (yvals.size()==5)
        {
			y2intres *= (maxy-miny)/12.0;
			if (y1ind==1 or y1ind==3) result += 4.0*y2intres;
			else if (y1ind==2) result += 2.0*y2intres;
			else result += y2intres;
		}
		else if (yvals.size()==1)
			result += y2intres;
		else
			cerr <<  "Unknown amount of yvals at "<< LINEINFO << endl;
    }
    if (yvals.size()==3)
		result *= (maxy-miny)/6.0;
	else if (yvals.size()==5)
		result *= (maxy-miny)/12.0;
	else if (yvals.size()==4)
		result *= (maxy-miny)/8.0;
    
    return result*2.0*M_PI; //2pi from one angular integral
}

double Inthelperf_y1(double y1, void* p)
{
    dSigma_full_helper* par = (dSigma_full_helper*)p;
    par->y1=y1;

    gsl_function fun;
    fun.function=Inthelperf_y2;
    fun.params=par;
    
    double result,abserr;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(YINT_INTERVALS);
    int status = gsl_integration_qag(&fun, par->miny, par->maxy,
            0, 0.01, YINT_INTERVALS, GSL_INTEG_GAUSS15, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);

    return result;
}

double Inthelperf_y2(double y2, void* p)
{
    dSigma_full_helper* par = (dSigma_full_helper*)p;
    par->y2=y2;

    gsl_function fun;
    fun.function=Inthelperf_cp_pt1;
    fun.params=par;
    
    
    double result,abserr;
            
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(PTINT_INTERVALS_CP);
    int status = gsl_integration_qag(&fun, par->minpt1, par->xs->MaxPt(),
            0, PTINT_ACCURACY_CP, PTINT_INTERVALS_CP, GSL_INTEG_GAUSS15, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        cerr << "pt1 int failed at " << LINEINFO << " intresult: " << result
            << " relerr " << std::abs(abserr/result) << " phi: " << par->phi << endl;
    }


    return result;
}



double Inthelperf_cp_pt1(double pt1, void* p)
{
    dSigma_full_helper* par = (dSigma_full_helper*)p;
    par->pt1=pt1;

    gsl_function fun;
    fun.function = Inthelperf_cp_pt2;
    fun.params=par;
    
    double minpt=0, maxpt=0;
    if (STAR)
    {
		minpt=1; maxpt=pt1;
	} else{
		minpt=0.5; maxpt=0.75;
	}

    double result,abserr;
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(PTINT_INTERVALS_CP);
    int status = gsl_integration_qag(&fun, minpt, maxpt,
            0, PTINT_ACCURACY_CP, PTINT_INTERVALS_CP, GSL_INTEG_GAUSS15, workspace,
            &result, &abserr);
    gsl_integration_workspace_free(workspace);
/*
    if (status)
    {
        cerr << "pt2 int failed at " << LINEINFO << " result: " << result
            << " relerr " << std::abs(abserr/result) << " pt1: " << pt1 << endl;
    }
*/
    return result;
}

double Inthelperf_cp_pt2(double pt2, void* p)
{
    dSigma_full_helper* par = (dSigma_full_helper*)p;

    return par->pt1*pt2*par->xs->dSigma_full(par->pt1, pt2, par->y1, par->y2, par->phi, par->sqrts,
                        par->deuteron);
}

/*
 * Funktion G_{x_A} = \int dr S(r) J_1(k*r) (m=0 or if gluon channel)
 * Or with nonzero m,z
 * mz \int dr r S(r) K_1(mzr) J1(kr)
 * as in ref. 0708.0231 but w.o. vector k/|k|
 * Default value of z=0
 */
struct G_helper { double y; AmplitudeLib* N; double kt; double z; CrossSection2 *xs; Channel channel; };
double G_helperf(double r, void* p);
double G_helperf_smallpt(double r, void* p);
double CrossSection2::G(double kt, double x, double z)
{

    
    G_helper helper;
    helper.y = std::log(N->X0()/x); 
    helper.channel=channel;
    helper.N=N; helper.kt=kt; helper.z=z; helper.xs=this;
    
    if (kt < 1e-3)
    {
		gsl_function fun; fun.function=G_helperf_smallpt; fun.params=&helper;
		gsl_integration_workspace *workspace 
		 = gsl_integration_workspace_alloc(10);
		 double gslresult,abserr;
		int status = gsl_integration_qag(&fun, N->MinR(), N->MaxR(),
				0, 0.001, 10, GSL_INTEG_GAUSS51, workspace,
				&gslresult, &abserr);
		gsl_integration_workspace_free(workspace);
    
		return gslresult;
		
	}
	
	if ( std::abs(kt - gcachek) < 0.00001 and std::abs(x/gcachex-1.0) < 0.0001 and std::abs(z-gcachez) < 0.0001)
        return gcacheval;

    set_fpu_state();
    init_workspace_fourier(FOURIER_ZEROS);   // number of bessel zeroes, max 2000
    double result = fourier_j1(std::abs(kt), G_helperf, &helper);
    gcachek=kt; gcachex = x; gcachez = z;
    gcacheval=result;
    

    
    return result;
    
    

}

double G_helperf(double r, void *p)
{
    G_helper* par = (G_helper*) p;
    double M_Q = par->xs->M_Q();
    // Massless case
    if (M_Q<1e-5 or par->z<1e-5 or par->channel==G_QQ or par->channel==G_GG)
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

double G_helperf_smallpt(double r, void* p)
{
	G_helper* par = (G_helper*) p;
	return par->kt*r/2.0*G_helperf(r, p);
}

/*
 * Funktion H_{x_A} = m*z^2*\int dr r S(r) K_0(mzr)J_0(kr)
 * as in ref. 0708.0231 
 */
double H_helperf(double r, void* p);
double CrossSection2::H(double kt, double x, double z)
{
    if (std::abs(kt - gcachek) < 0.001 and std::abs(x/hcachex-1.0) < 0.001 and std::abs(z-hcachez) < 0.001)
        return hcacheval;
        
    if (M_Q() < 1e-5) return 0;
    G_helper helper;
    helper.y = std::log(N->X0()/x);
    helper.N=N; helper.kt=kt; helper.z=z; helper.xs=this;

    set_fpu_state();
    init_workspace_fourier(700);   // number of bessel zeroes, max 2000

    double result = fourier_j0(std::abs(kt), H_helperf, &helper);

    hcacheval=M_Q()*SQR(z)*result;
    hcachek=kt; hcachex=x; hcachez=z;
    return hcacheval;

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

/*
 * Wave function product summed over spins and helicities
 * Momentum space
 * zfrac: momentum fraction z
 */
double CrossSection2::WavefSqr_k(double k1, double k2, double k1_dot_k2, double zfrac)
{
	return 2.0*(k1_dot_k2*(1.0+SQR(1.0-zfrac)) + SQR(M_Q()*zfrac*zfrac))
		/ ( (SQR(k1)+SQR(M_Q()*zfrac))*(SQR(k2)+SQR(M_Q()*zfrac)) );
}

CrossSection2::CrossSection2(AmplitudeLib* N_, PDF* pdf_,FragmentationFunction* frag)
{
    N=N_; pdf=pdf_; fragfun=frag;
    gcacheval=-1;
    gcachek=-1; gcachex=-1; gcachez=-1;
    hcacheval=-1; hcachek=-1; hcachex=-1; hcachez=-1;
    transform=NULL;
    mcintpoints=1e7;
    gsl_rng_env_setup();
    m_q=0.14;

    ptinterpolator2d = NULL;
    ptinterpolator2d_rev = NULL;
    ptinterpolator2d_correction = NULL;
    ptinterpolator2d_rev_correction = NULL;

	///paperi maxpt 7.5
	double minpt=0.5, maxpt=8, ptstep=0.5;
    for (double pt=minpt; pt<=maxpt; pt+=ptstep)
    {
        ptvals.push_back(pt);
        std::stringstream s; s << pt;
        ptstrings.push_back(s.str());
    }

    apply_corrections = false;
    finite_nc=false;
    
        
        

    //fileprefix = "final_result/ircutoff_lambdaqcd/mv1_qs072/largenc/";
    fileprefix = "gluon/q_vs_g/lhc/ggg/";
 
    fileprefix_cor = "NOT USED";
    
	postfix = "";//"_largenc";
	#ifdef USE_MPI
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	if (id==0)
	#endif
		    cout << "# Loading data from " << fileprefix << " postfix " << postfix << " minpt " << minpt << " maxpt " << maxpt << " ptstep " << ptstep << endl;	
}


double CrossSection2::Delta(double pt1, double pt2, double theta)
{
    return std::sqrt( SQR(pt1) + SQR(pt2) + 2.0*pt1*pt2*std::cos(theta) );
}

double CrossSection2::Z(double pt1, double pt2, double y1, double y2)
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
            delete ptinterpolators_rev[i][j];
        }
    }

    for (uint i=0; i<ptinterpolators_correction.size(); i++)
    {
        for (uint j=0; j<ptinterpolators_correction[i].size(); j++)
        {
            delete ptinterpolators_correction[i][j];
            delete ptinterpolators_rev_correction[i][j];
        }
    }

    if (ptinterpolator2d != NULL)
        delete ptinterpolator2d;
    if (ptinterpolator2d_rev != NULL)
        delete ptinterpolator2d_rev;
    if (ptinterpolator2d_correction != NULL)
        delete ptinterpolator2d_correction;
    if (ptinterpolator2d_rev_correction != NULL)
        delete ptinterpolator2d_rev_correction;
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
    return ptvals[ptvals.size()-1]*0.999;
}

Interpolator2D * CrossSection2::Ptinterpolator2d()
{
    return ptinterpolator2d;
}
Interpolator2D * CrossSection2::Ptinterpolator2d_rev()
{
    return ptinterpolator2d_rev;
}

AmplitudeLib* CrossSection2::GetN()
{
    return N;
}

void CrossSection2::SetFiniteNc(bool fnc)
{
	finite_nc=fnc;
}

bool CrossSection2::FiniteNc()
{
	return finite_nc;
}

 // The following piece of code is so ugly it is hidden at the bottom of this
 // huge file. Sry.
/*
 * Load \Delta \phi data from files pt1_#_pt2_#_y1_#_y2_#, interpolate in pt and phi
 *
 * Returns 0 if no error occurred, otherwise <0
 */
int CrossSection2::LoadPtData(double y1, double y2)
{
	//cout << "Loading: y1=" << y1 <<", y2=" << y2 << endl;
   // Load datafiles 
   // Generate interpolators

    std::stringstream tmp;
    std::string y1str, y2str;
    tmp << y1; y1str = tmp.str();
    std::stringstream tmp2; 
    tmp2 << y2; y2str = tmp2.str();   

    for (uint i=0; i<ptinterpolators.size(); i++)
    {
        for (uint j=0; j<ptinterpolators[i].size(); j++)
        {
            delete ptinterpolators[i][j];
            delete ptinterpolators_rev[i][j];
        }
        ptinterpolators[i].clear();
        ptinterpolators_rev[i].clear();
    }
    ptinterpolators.clear();
    ptinterpolators_rev.clear();
    for (uint i=0; i<ptinterpolators_correction.size(); i++)
    {
        for (uint j=0; j<ptinterpolators_correction[i].size(); j++)
        {
            delete ptinterpolators_correction[i][j];
            delete ptinterpolators_rev_correction[i][j];
        }
        ptinterpolators_correction[i].clear();
        ptinterpolators_rev_correction[i].clear();
    }
    ptinterpolators_correction.clear();
    ptinterpolators_rev_correction.clear();

   int points=ptvals.size();
   for (int pt1ind=0; pt1ind<points; pt1ind++)
   {
       std::vector<Interpolator*> tmpinterpolators;
       std::vector<Interpolator*> tmpinterpolators_rev;
       std::vector<Interpolator*> tmpinterpolators_cor;
       std::vector<Interpolator*> tmpinterpolators_rev_cor;
       for (int pt2ind=0; pt2ind<points; pt2ind++)
       {
		   ///TODO: Hardcoded rhic kinematics
		   double x = xh(ptvals[pt1ind], ptvals[pt2ind], y1, y2, 5020);
		   double x_rev = xh(ptvals[pt1ind], ptvals[pt2ind], y2, y1, 5020);
            std::stringstream fname, fname_rev, fname_cor, fname_rev_cor;
            /*fname << fileprefix << "pt1_" << ptstrings[pt1ind] << "_pt2_"
                << ptstrings[pt2ind] << "_y1_" << y1str << "_y2_" << y2str << postfix;
            fname_rev << fileprefix << "pt1_" << ptstrings[pt1ind] << "_pt2_"
                << ptstrings[pt2ind] << "_y1_" << y2str << "_y2_" << y1str << postfix;
            */
               
            fname << fileprefix << "y1_" << y1str <<"_y2_" << y2str << "/pt1_"
                << ptstrings[pt1ind] <<"_pt2_" << ptstrings[pt2ind];
            fname_rev << fileprefix << "y1_" << y2str <<"_y2_" << y1str << "/pt1_"
                << ptstrings[pt1ind] <<"_pt2_" << ptstrings[pt2ind];
              
		

            //cout << "# Loading files " << fname.str() << " and " << fname_rev.str() << endl;
            
            std::ifstream file(fname.str().c_str());
            std::ifstream file_rev(fname_rev.str().c_str());
            
            
            std::vector<double> dphi;
            std::vector<double> dphi_rev;
            std::vector<double> xs;
            std::vector<double> xs_rev;

            if (!file.is_open()  )
            {
				if (x<1)
					cerr << "Can't open " << fname.str()  << " " << LINEINFO << ", assuming=0" << " (x_h=" << x <<")" << endl;
                dphi.push_back(0.5); dphi.push_back(1.4); dphi.push_back(2.5); dphi.push_back(3.141592);
                xs.push_back(0); xs.push_back(0); xs.push_back(0);xs.push_back(0);
                //return -1;
            }
            while(!file.eof() and file.is_open())
            {
                string line;
                getline(file, line);
                if (line.substr(0,1)!="#" and line.length()>3)  // not comment/empty
                {
                    std::istringstream iss(line);
                    string angle,val;
                    iss >> angle; iss >> val;
                    xs.push_back(StrToReal(val));
                    dphi.push_back(StrToReal(angle));
                }
            }
            if (dphi.size()<3)
            {
				if (x<1)
					cerr << "No enough data points in " << fname.str() << " (x_h=" << x << ", assuming=0): " << LINEINFO << endl;
				dphi.clear(); xs.clear();
				dphi.push_back(0.5); dphi.push_back(1.4); dphi.push_back(2.5); dphi.push_back(3.141592);
                xs.push_back(0); xs.push_back(0); xs.push_back(0);xs.push_back(0);
			}
            if (dphi[dphi.size()-1]<3.1)
            {
				if (x<1)
					cerr << "Max dphi in file " << fname.str() << " is " << dphi[dphi.size()-1] << " at " << LINEINFO  << " (x_h=" << x << ", assuming=0)" << endl;
				dphi.clear(); dphi.push_back(0.5); dphi.push_back(1); dphi.push_back(1.5); dphi.push_back(3.141592);
				xs.clear(); xs.push_back(0);xs.push_back(0);xs.push_back(0);xs.push_back(0);
			}
			if (dphi[0]>1)
				cerr << "Min dphi in file " << fname.str() << " is " << dphi[0] << " (x_h=" << x << ") " << LINEINFO << endl;
            
            if (!file_rev.is_open()  )
            {
				if (x_rev<1)
					cerr << "Can't open " << fname_rev.str()  << " " << LINEINFO << ", assuming=0" << " (x_h=" << x_rev << "): " << endl;
                dphi_rev.push_back(0.5); dphi_rev.push_back(1.4); dphi_rev.push_back(2.5); dphi_rev.push_back(3.141592);
                xs_rev.push_back(0); xs_rev.push_back(0); xs_rev.push_back(0);xs_rev.push_back(0);
                //return -1;
            }
            while (!file_rev.eof() and file_rev.is_open())
            {
                string line;
                getline(file_rev, line);
                if (line.substr(0,1)!="#" and line.length()>3)  // not comment/empty
                {
                    std::istringstream iss(line);
                    string angle,val;
                    iss >> angle; iss >> val;
                    xs_rev.push_back(StrToReal(val));
                    dphi_rev.push_back(StrToReal(angle));
                }                
            }
            if (dphi_rev.size()<3)
            {
				//if (x_rev<1)
					cerr << "Not enough data points in " << fname_rev.str() << ": " << LINEINFO << " (x_h=" << x_rev << ", assuming=0): " << endl;
				dphi_rev.clear(); xs_rev.clear();
				dphi_rev.push_back(0.5);  dphi_rev.push_back(1); dphi_rev.push_back(1.5); dphi_rev.push_back(3.141592);
				xs_rev.clear(); xs_rev.push_back(0);xs_rev.push_back(0);xs_rev.push_back(0);xs_rev.push_back(0);
			}
            if (dphi_rev[dphi_rev.size()-1]<3.1)
            {
				//if (x_rev<1)
					cerr << "Max dphi in file " << fname_rev.str() << " is " << dphi_rev[dphi_rev.size()-1] <<  " (x_h=" << x_rev << ", assuming=0): " << endl;
				dphi_rev.clear(); 		dphi_rev.push_back(0.5);  dphi_rev.push_back(1); dphi_rev.push_back(1.5); dphi_rev.push_back(3.141592);
				xs_rev.clear(); xs_rev.push_back(0);xs_rev.push_back(0);xs_rev.push_back(0);xs_rev.push_back(0);
			}
			if (dphi[0]>1)
				cerr << "Min dphi in file " << fname.str() << " is " << dphi[0] << " " << LINEINFO << endl;


            file.close();
            file_rev.close();

            // Mirror \dphi if necessary
            if (dphi[dphi.size()-1] < 1.01*M_PI)
            {
                for (int i=dphi.size()-2; i>=0; i--)
                {
                    dphi.push_back(M_PI + (M_PI-dphi[i]) );
                    xs.push_back(xs[i]);
                }
            }
            // Mirror \dphi if necessary
            if (dphi_rev[dphi_rev.size()-1] < 1.01*M_PI)
            {
                for (int i=dphi_rev.size()-2; i>=0; i--)
                {
                    dphi_rev.push_back(M_PI + (M_PI-dphi_rev[i]) );
                    xs_rev.push_back(xs_rev[i]);
                }
            }
            
			//cerr << fname.str() << endl;
            Interpolator *tmpinterp = new Interpolator(dphi, xs);
            tmpinterp->Initialize();
            tmpinterpolators.push_back(tmpinterp);
			//cerr << fname_rev.str() << endl;
            Interpolator *tmpinterp_rev = new Interpolator(dphi_rev, xs_rev);
            tmpinterp_rev->Initialize();
            tmpinterpolators_rev.push_back(tmpinterp_rev);


        }
        ptinterpolators.push_back(tmpinterpolators);
        ptinterpolators_rev.push_back(tmpinterpolators_rev);

   }
   //cout << "# Data loaded" << endl;

  

    return 0;
}

void CrossSection2::SetChannel(Channel c_)
{
	channel=c_;
}

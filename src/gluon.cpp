/*
 * Computes g->gg scattering
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */
#include "xs.hpp"
#include "config.hpp"
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
extern "C"
{
    #include <fourier/fourier.h>
}


const double IR_CUTOFF=LAMBDAQCD;
void pVegas_gluon_wrapper(double x[6], double f[1], void* par);

const double GGG_MAXR = 80;

double CrossSection2::GluonCorrectionTerm(double pt1, double pt2, double ya, double phi,
        double z)

{
    std::time_t start = std::time(NULL);
    #ifdef USE_MPI
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if(id==0) 
    #endif
    cout << "# Starting GLUON MC integral, "
        << " finite-Nc " << FiniteNc() << " ir-cutoff " << IR_CUTOFF << " GeV " << endl;
    


    if (!N->InterpolatorInitialized(ya))
        cerr << "Interpolator is not initialized for y=" << ya << "!!! " << LINEINFO << endl;
    
    /*
     * We need to integrate over u,u',r and their corresponding angles
     */
    // vec[0]= ln u1, vec[1]=ln u2, vec[2]= ln r
    // vec[3]=theta_u1, vec[4]=theta_u2
    // vec[5]=theta_r
    // pt1=k, pt2=q

    N->SetOutOfRangeErrors(false);
    double minr=std::log(0.0001); double maxr=std::log(GGG_MAXR);
    //double minr=0; double maxr=100;

       
    double res=0;
    // Constants taken out from integral in order to reduse number of
    // products in integrand
    double constants=0;
    if (channel != G_GG)
    {
		cerr << "Wrong channel " << channel <<" and we are computing g->gg " << LINEINFO << endl;
		return 0;
	}


	constants=1.0/std::pow(2.0*M_PI, 6) * SQR(2.0*M_PI) * 4.0
			* ( z/(1.0-z) + (1.0-z)/z + z*(1.0-z) );


        Inthelper_correction helper;
            helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
            helper.N=N; helper.z=z; helper.calln=0; helper.monte=3;
            helper.xs=this; helper.finite_nc=FiniteNc();
            


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

    time_t t;
    time(&t);
    
    // set up the grid (init = 0) with 5 iterations of 1000 samples,
    // no need to compute additional accumulators (fcns = 1),
    if (id==0)
    {
        cout <<"# Constants: " << constants << " workers " << workers << endl;
        cout << "# First round at " << std::ctime(&t) ;
    }
    vegas(reg, dim, pVegas_gluon_wrapper,
        0, mcintpoints/100, 5,  0,
        functions, 0, workers,
        estim, std_dev, chi2a, &helper);

      // refine the grid (init = 1) 
      // collect in additional accumulators (fcns = FUNCTIONS),
      time(&t);
    if (id==0)
        cout << "# 2nd round at " << std::ctime(&t);  
    vegas(reg, dim, pVegas_gluon_wrapper,
            1, mcintpoints/10, 5, NPRN_RESULT,
            functions, 0, workers,
            estim, std_dev, chi2a, &helper);
      // final sample, inheriting previous results (init = 2)
    time(&t);
    if (id==0)
        cout << "# 3rd round at " << std::ctime(&t) ;  
    vegas(reg, dim, pVegas_gluon_wrapper,
            2, mcintpoints, 2,   NPRN_RESULT,
            functions, 0,  workers,
            estim, std_dev, chi2a, &helper);

    if (id==0)
    {
      cout << "# phi=" << phi << " result: " << estim[0] << " +/- " << std_dev[0]
        << " effect " << estim[0]*constants << " +/- " << std_dev[0]*constants << " relerr " << std::abs(std_dev[0]/estim[0]) << " time "
       << (std::time(NULL)-start)/60.0/60.0 << " h"<< endl;
	}
    if (std::abs(std_dev[0]/estim[0]) > 0.1)
    {
		if (id==0)
		{
            //cout << "#!!!!!!!!!!!!!!! INTEGRAL DOESN'T CONVERGE!?!?!?" << endl;
            cout <<"# INTEGRAL DOESN'T CONVERGE! phi=" << phi
                << " result: " << estim[0] << " +/- " << std_dev[0]
                << " relerr " << std::abs(std_dev[0]/estim[0]) << " pt1 " << pt1 << " pt2 "
                << pt2 << " z " << z << endl;//" increasing mcintpoints from " << mcintpoints << " to " << mcintpoints*3 << endl;
		}
		return -99999999;
    }
	//else if (std::abs(std_dev[0]/estim[0])<0.01) // Can reduce number of mcintpoints
		//mcintpoints /= 2.0;
    

    
    
    res = estim[0]*constants;
    #endif
    return res;
}

double Inthelperf_gluon_correction(double* vec, size_t dim, void* p)
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
    
    bool finite_nc = par->finite_nc;
    if (finite_nc) cerr << "GLUON CHANNEL DOES NOT SUPPORT FINITE NC" << endl;
    
    if (isnan(u1) or isnan(u2))
    {
        cerr << "u1/u2 nan, WTF??? " << LINEINFO << endl;
    }

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
    double s_r=N->S(r,y);
    
  
    // r - u1
    double r_m_u1 = std::sqrt(SQR(r) + SQR(u1) - 2.0*u1_dot_r);
    double s_r_m_u1 = N->S(r_m_u1, y);
    // r + u1
    double r_p_u1 = std::sqrt(SQR(r) + SQR(u1) + 2.0*u1_dot_r);
    double s_r_p_u1 = N->S(r_p_u1, y);
    // r + u'
    double r_p_u2 = std::sqrt(SQR(r) + SQR(u2) + 2.0*u2_dot_r);
    double s_r_p_u2 = N->S(r_p_u2, y);
    // r - u'
    double r_m_u2 = std::sqrt(SQR(r) + SQR(u2) - 2.0*u2_dot_r);
    double s_r_m_u2 = N->S(r_m_u2, y);
    // r - u + u' = r + (u'-u)
    double r_m_u1_p_u2 = std::sqrt( SQR(r) + (SQR(u1) + SQR(u2) - 2.0*u1_dot_u2)
                    + 2.0*(u2_dot_r - u1_dot_r) );
	double s_r_m_u1_p_u2 = N->S(r_m_u1_p_u2, y);
	
		
		// Quadrupole, large Nc
		double f1 = s_r_m_u1 * s_r_p_u2 / (s_r_m_u1_p_u2 * s_r);
		double f2 = s_u1 * s_u2 / (s_r_m_u1_p_u2 * s_r);
		double loglog=0;     
		loglog = std::log(f1) / std::log(f2);
		// If u1,u2>>r loglog->1
		// (but this is a DPS contribution and removed later)
		if (s_u1==0 and s_u2==0 and s_r>0 )
			loglog=1;
		// All other NaN limits vanish
			
		if (isnan(loglog) )
			loglog=0;
		if (isinf(loglog))
			loglog=0;
					
		double q = s_u1*s_u2 - loglog * (s_u1 * s_u2 - s_r * s_r_m_u1_p_u2);
		
		// e^(-i(k+q).r) e^(iq.(u-u')) S(r)S(r-u+u')Q(r,u,u') - DPS
		double result = s_r * s_r_m_u1_p_u2 * q ; //
		result -= s_r*s_r*s_u1*s_u2;		// Subtract this from here and add it to GluonGluonGluon() where it can be computed analytically
		
		result -=  s_r*s_r*s_r_m_u1_p_u2*s_r_m_u1_p_u2 
					* step(u1-1.0/IR_CUTOFF)* step(u2-1.0/IR_CUTOFF);	// \theta(u1-1/CUTOFF)\tehta(u2-1/CUTOFF)
		
		result *= std::cos(-pt1_dot_r  - pt2_dot_r + pt2_dot_u1 - pt2_dot_u2);
		/*
		///// DEBUG: FULL result, not just "correction"
		result += - s_r*s_r_m_u2 * s_u2 * std::cos( -pt1_dot_r - pt2_dot_r + pt1_dot_u2 - (1.0-z)*pt1_dot_u1 + z*pt2_dot_u1 )
			- s_r*s_r_p_u1 * s_u1 * std::cos( -pt1_dot_r - pt2_dot_r - pt1_dot_u1 + (1.0-z)*pt1_dot_u2 - z*pt2_dot_u2 )
			+ s_r*s_r * std::cos( -pt1_dot_r - pt2_dot_r + (1.0-z)*pt1_dot_u2 - z*pt2_dot_u2 - (1.0-z)*pt1_dot_u1 + z*pt2_dot_u1 );
		*/
		// wavef
		result *= u1_dot_u2 / (SQR(u1*u2) );
		
		result *= SQR(u1*u2*r);	// d^2 u1 d^2 u2 d^2 r = du1 du2 dr u1 u2 r, and then jacobian u->ln u
		//result *= u1*u2*r;
		
		return result;
		
}

void pVegas_gluon_wrapper(double x[6], double f[1], void* par)
{
    f[0] = Inthelperf_gluon_correction(x, 6, par);
}
void pVegas_gluon_uncorrected_wrapper(double x[4], double f[2], void* par);


/*
 * Calculates the uncorrected part of the gluon channel
 * Integrals are 4 or less dimensional
 * 
 * No other prefactors included than (2\pi)^(-6)
 */

struct Inthelper_gluon
{
	double pt1, pt2, ya, phi;
	double u,r,theta_u;
	bool x;	// true if x-coord, false if y
	double z;
	AmplitudeLib* N;
};
void Inthelperf_gluon_uncorrected(double vec[4], double f[2], void* p);

double Inthelperf_gluon_uncorrected_u(double u, void* p);
double Inthelperf_gluon_uncorrected_r(double r, void* p);
double Inthelperf_gluon_uncorrected_theta_u(double theta_u, void* p);
double Inthelperf_gluon_uncorrected_theta_r(double theta_r, void* p);
const int GLUONINTPOINTS = 4;

double CrossSection2::GluonGluonGluon(double pt1, double pt2, double y1, double y2, double phi, double sqrts)
{
	#ifdef USE_MPI
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	#endif

    double tmpz = Z(pt1, pt2, y1, y2);
    double tmpxa = xa(pt1, pt2, y1, y2, sqrts);

    double ya = std::log(N->X0()/tmpxa);
    if (ya<0)
        cerr << "Negative rapidity interval at " << LINEINFO << endl;
    
    if (!N->InterpolatorInitialized(ya))
        N->InitializeInterpolation(ya);
    double delta = Delta(pt1,pt2,phi);
    
    // \kappa^2 = (k - z\delta)^2 = (1-z)^2 pt1^2 + z^2 pt2^2 - 2*z*(1-z)*pt1*pt2*cos \phi
    double kzdeltasqr = SQR(1.0-tmpz)*SQR(pt1) + SQR(tmpz*pt2) - 2.0*tmpz*(1.0-tmpz)
                                * pt1*pt2*std::cos(phi);
    double kappa = std::sqrt(kzdeltasqr);

    double g=0,f2=0;
    
    ////DEBUG
    //return tmpz*(1.0-tmpz)*Nc/2.0*GluonCorrectionTerm(pt1, pt2, ya, phi, tmpz);
 
    f2=N->S_k(delta, ya, false, 2.0)/SQR(2.0*M_PI);	// 2d ft of S(r)^2 (adjoint rep=false) / (2pi)^2
    //g = G(delta, tmpxa);	// \int dr S(r) J_1(k*r) = 1/2pi \int d^2 r e^(iq.r)S(r)
    g = G(pt2, tmpxa);
    
    double wf_z = 4.0*SQR(2.0*M_PI)* ( tmpz/(1.0-tmpz) + (1.0-tmpz)/tmpz + tmpz*(1.0-tmpz) );
    
    double result = g * g * f2  * 1.0/(SQR(2.0*M_PI) );		// / (2pi)^2 as we have (2pi)^2 in wf_z
    if (id==0) cout << "# kappa: " << kappa << " f2 " << f2 << " g " << g << endl;
    if (id==0) cout << "# 1: " << result << endl;
    result += 1.0/SQR(kappa) * f2 / std::pow(2.0*M_PI,2);
    if (id==0) cout << "# 2: " << 1.0/SQR(kappa) * f2 / std::pow(2.0*M_PI,2) << endl;
    
  
    //cout << "wo interferece: " << result << endl;
    
    
    /////////////////////////
    // then subtract interference-like terms, a bit more tricky...
	
	Inthelper_gluon par;
	par.pt1=pt1; par.pt2=pt2; par.ya=ya; par.N=N; par.phi=phi;
	par.z=tmpz;
	std::time_t start = std::time(NULL);
	
	/*
	gsl_function fun;
    fun.params = &par;
    fun.function = Inthelperf_gluon_uncorrected_u;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(GLUONINTPOINTS);
    double intresult,abserr;

    int status = gsl_integration_qag(&fun, MINLNR, MAXLNR,
            0, 0.1, GLUONINTPOINTS, GSL_INTEG_GAUSS51, workspace, &intresult, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl; 
        
     cout << "xres " << intresult << " pm " << abserr << endl;  
   */
     
	const int dim=4;
    double reg[2*dim];
    int functions = 2;
    N->SetOutOfRangeErrors(false);
    double minr=std::log(0.0001); double maxr=std::log(GGG_MAXR);
    reg[0]=minr; reg[1]=minr;
    reg[dim]=maxr; reg[dim+1]=maxr; 
    reg[2]=0; reg[3]=0; 
    reg[dim+2]=2.0*M_PI; reg[dim+3]=2.0*M_PI;
    
    double estim[2];   // estimators for integrals 
    double std_dev[2]; // standard deviations              
    double chi2a[2];   // chi^2/n       
    int workers = 1;
    
    MPI_Comm_size(MPI_COMM_WORLD, &workers);
    workers--;
    if (workers==0)
    {
        cerr << "No workers! MPI needs at least two threads! " << endl;
        return -1;
    }
    
    time_t t;
    time(&t);
    
    // set up the grid (init = 0) with 5 iterations of 1000 samples,
    // no need to compute additional accumulators (fcns = 1),
    vegas(reg, dim, pVegas_gluon_uncorrected_wrapper,
        0, mcintpoints/100, 5,  0,
        functions, 0, workers,
        estim, std_dev, chi2a, &par);

      // refine the grid (init = 1) 
      // collect in additional accumulators (fcns = FUNCTIONS),
      time(&t);
    if (id==0)
        cout << "# 2nd round at " << std::ctime(&t);  
    vegas(reg, dim, pVegas_gluon_uncorrected_wrapper,
            1, mcintpoints/10, 5, NPRN_RESULT,
            functions, 0, workers,
            estim, std_dev, chi2a, &par);
      // final sample, inheriting previous results (init = 2)
    time(&t);
    if (id==0)
        cout << "# 3rd round at " << std::ctime(&t) ;  
    vegas(reg, dim,pVegas_gluon_uncorrected_wrapper,
            2, mcintpoints, 2,   NPRN_RESULT,
            functions, 0,  workers,
            estim, std_dev, chi2a, &par);

    if (id==0)
    {
      cout << "# phi=" << phi << " result: " << estim[0] << " +/- " << std_dev[0]
        << " relerr " << std::abs(std_dev[0]/estim[0]) << " time "
       << (std::time(NULL)-start)/60.0/60.0 << " h"<< endl;
	}
    if (std::abs(std_dev[0]/estim[0]) > 0.1)
    {
		if (id==0)
		{
            //cout << "#!!!!!!!!!!!!!!! INTEGRAL DOESN'T CONVERGE!?!?!?" << endl;
            cout <<"# INTEGRAL DOESN'T CONVERGE! phi=" << phi
                << " result: " << estim[0] << " +/- " << std_dev[0]
                << endl;
		}
		return -99999999;
    }
	
	/*result -= -2.0/std::pow(2.0*M_PI, 5) * ( 
		(1.0-tmpz)*pt1*(std::cos(phi)*estim[0] + std::sin(phi)*estim[1])
		+ tmpz*pt2*estim[0] );*/
	double interf = estim[0];
	result += interf;

	if (id==0)
	cout << "# interference: " << interf << endl;
	/*cout << "#interferece contrib: "<< 2.0/std::pow(2.0*M_PI, 5) * ( 
		(1.0-tmpz)*pt1*(std::cos(phi)*estim[0] + std::sin(phi)*estim[1])
		+ tmpz*pt2*estim[0] ) << endl;*/
	
	result *= wf_z;
	
	double correction=0;
	correction = GluonCorrectionTerm(pt1, pt2, ya, phi, tmpz);
	if (id==0)
	cout << "# uncorrected: " << result << " correction " << correction << endl;
	result +=correction;
	    
    // Normalization factor that is not in CrossSection2::dSigma_full
    // where result is multiplied by z(1-z)S_T\alpha_s/(4pi^2)
    result *= tmpz*(1.0-tmpz)*Nc/2.0;
	return result;
	
}

double Inthelperf_gluon_uncorrected_u(double lnu, void* p)
{
	Inthelper_gluon* par = (Inthelper_gluon*)p;
	par->u = std::exp(lnu);
	gsl_function fun;
    fun.params = par;
    fun.function = Inthelperf_gluon_uncorrected_r;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(GLUONINTPOINTS);
    double result,abserr;

    int status = gsl_integration_qag(&fun, MINLNR, MAXLNR,
            0, 0.1, GLUONINTPOINTS, GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl; 
    
    return result;
	
}

double Inthelperf_gluon_uncorrected_r(double lnr, void* p)
{
	Inthelper_gluon* par = (Inthelper_gluon*)p;
	par->r = std::exp(lnr);
	gsl_function fun;
    fun.params = par;
    fun.function = Inthelperf_gluon_uncorrected_theta_u;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(GLUONINTPOINTS);
    double result,abserr;

    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI,
            0, 0.1, GLUONINTPOINTS, GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl; 

    
    return result;
	
}

double Inthelperf_gluon_uncorrected_theta_u(double theta_u, void* p)
{
	Inthelper_gluon* par = (Inthelper_gluon*)p;
	par->theta_u = theta_u;
	gsl_function fun;
    fun.params = par;
    fun.function = Inthelperf_gluon_uncorrected_theta_r;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(GLUONINTPOINTS);
    double result,abserr;

    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI,
            0, 0.1, GLUONINTPOINTS, GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl; 
    
    return result;
	
}

double Inthelperf_gluon_uncorrected_theta_r(double theta_r, void* p)
{
	Inthelper_gluon* par = (Inthelper_gluon*)p;
	double theta_u = par->theta_u;
	double u = par->u; double r = par->r;
	double pt1 = par->pt1; double pt2 = par->pt2;
	double phi = par->phi;
	double z = par->z;
	
	// all angles are measured from pt1 to positive (counterclockwise) direction
	double pt1_dot_r = pt1*r*std::cos(theta_r);
	double pt2_dot_r = pt2*r*std::cos(theta_r - phi);
	double pt1_dot_u = pt1*u*std::cos(theta_u);
	double pt2_dot_u = pt2*u*std::cos(theta_u - phi);
	double kappa_sqr = SQR(1.0-z)*SQR(pt1) + SQR(z*pt2) - 2.0*(1.0-z)*z*pt1*pt2*std::cos(phi);
	double kappa = std::sqrt(kappa_sqr);
	
	double kappa_dot_u = (1.0-z)*pt1_dot_u - z*pt2_dot_u;
	
	double r_m_u = std::sqrt( SQR(r) + SQR(u) - 2.0*r*u*std::cos(theta_r - theta_u));
	
	double sss = par->N->S(r, par->ya) * par->N->S(u, par->ya) * par->N->S(r_m_u, par->ya);
	
	double result = -2.0 * 2.0*M_PI / std::pow(2.0*M_PI, 6) * kappa_dot_u / (SQR(kappa*u)) 
					* std::sin( -pt1_dot_r - pt2_dot_r + pt1_dot_u) * sss;
	
	result *= SQR(u*r);	// jacobian, d^2 u d^2 r and integrate over log u, log r
	return result;	
	
	//wanha:
	/*
	double dotprod = std::sin( -(pt1*r*std::cos(phi-theta_r) + pt2*r*std::cos(theta_r)) - pt2*u*std::cos(theta_u) );
	
	double result = 0;
	if (par->x)
		result = sss*dotprod * std::cos(theta_u) / u;
	else
		result = sss*dotprod * std::sin(theta_u) / u;
	//f[1] = dotprodcos * std::sin(theta_u) / u;
	
	result *= SQR(u*r);	// jacobian, d^2 u d^2 r and integrate over log u, log r
	*/
	return result;
	
}

void pVegas_gluon_uncorrected_wrapper(double x[4], double f[2], void* par)
{
	Inthelper_gluon* p = (Inthelper_gluon*) par;
	p->u = std::exp(x[0]);
	p->r = std::exp(x[1]);
	p->theta_u = x[2];
	p->x=true;
	f[0] = Inthelperf_gluon_uncorrected_theta_r(x[3], p);
	f[1]=0;
	//p->x=false;
	//f[1]= Inthelperf_gluon_uncorrected_theta_r(x[3], p);
}

	
/*****************************************************************
 * g -> qq channel
 *****************************************************************/
struct inthelper_gqq
{
	AmplitudeLib* N;
	double delta;
	double phi;
	double pt1;
	double k,q;
	double z;
	double ya;
};
double inthelperf_gqq_theta(double theta, void* p);
double inthelperf_gqq_p(double p1, void* par);
// pt1=quark, pt2=antiquark, but symmetric in exchange!

double CrossSection2::GluonQuarkQuark(double pt1, double pt2, double y1, double y2, double phi, double sqrts)
{
	double tmpz = Z(pt1, pt2, y1, y2);
    double tmpxa = xa(pt1, pt2, y1, y2, sqrts);
    // \kappa^2 = (k - z\delta)^2 = (1-z)^2 pt1^2 + z^2 pt2^2 - 2*z*(1-z)*pt1*pt2*cos \phi
    double kappasqr = SQR(1.0-tmpz)*SQR(pt1) + SQR(tmpz*pt2) - 2.0*tmpz*(1.0-tmpz)
                                * pt1*pt2*std::cos(phi);
    double kappa = std::sqrt(kappasqr);
    double delta = Delta(pt1,pt2,phi);
    
    double minpt=LAMBDAQCD;
    
    cout <<"# Computing g->qq channel, z=" << tmpz <<" kappa=" << kappa << " pt1 " << pt1 << " pt2 " << pt2  <<" delta " << delta << " ya " << std::log(N->X0() / tmpxa) << " dps removal cutoff: " << minpt << endl;
    
    inthelper_gqq helper;
    helper.N=N;
    helper.delta=delta; helper.phi=phi;
    helper.k=pt1; helper.q=pt2; helper.z=tmpz;
    helper.ya = std::log(N->X0() / tmpxa);
    gsl_function fun;
    fun.function=inthelperf_gqq_p;
    fun.params=&helper;
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(GLUONINTPOINTS);
    double result,abserr;

    int status = gsl_integration_qag(&fun, minpt, 10.0*std::max(pt1,pt2),
            0, 0.03, GLUONINTPOINTS, GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl; 

    
    return result;
    
	
}

double inthelperf_gqq_p(double p1, void* par_)
{
	inthelper_gqq *par = (inthelper_gqq*)par_;
	par->pt1=p1;
	
	gsl_function fun;
	fun.params=par;
	fun.function=inthelperf_gqq_theta;
	
	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(GLUONINTPOINTS);
    double result,abserr;
    


    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI,
            0, 0.03, GLUONINTPOINTS, GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl; 
    
    return result;
}

// Geometry: angles are measured from k
// angle between k and q is phi
// angle between k and p is theta

double inthelperf_gqq_theta(double theta, void* par_)
{
	
	inthelper_gqq *par = (inthelper_gqq*)par_;
	
	double delta = par->delta;
	double phi = par->phi;
	double k=par->k; double q = par->q;
	double p1 = par->pt1; 
	double z = par->z;
	double ya = par->ya;
	
	
	// p1 + q
	double p1_p_q = std::sqrt(
		p1*p1 + q*q + 2.0*p1*q*std::cos(phi-theta)
	);
	
	// p1-k
	double p1_m_k = std::sqrt(
		p1*p1 + k*k - 2.0*p1*k*std::cos(theta)
	);
	
	// kappa = (1-z)k-zq
	double kappa = std::sqrt( 
		SQR(1.0-z)*k*k + z*z*q*q - 2.0*z*(1.0-z)*k*q*std::cos(phi)
	);
	
	// kappa - p1 = (1-z)k - zq - p1
	double kappa_m_p1 = std::sqrt(
		kappa*kappa + p1*p1
		 - 2.0*( (1.0-z)*p1*k*std::cos(theta) - z*p1*q*std::cos(phi-theta) )
		);
	
	double result = par->N->S_k(p1_p_q, ya) * par->N->S_k(p1_m_k, ya) / std::pow(2.0*M_PI,4);
	result *= (SQR(z) + SQR(1.0-z)) * SQR(kappa_m_p1) / ( SQR(p1)*SQR(kappa) );
	return p1*result;
	
	
	
	/*
	// \Delta - p1 = k+q-p1
	double delta_m_p1 = std::sqrt( 
		k*k + q*q + p1*p1 - 2.0*k*p1*std::cos(theta) - 2.0*q*p1*std::cos(theta-phi)
				+ 2.0*k*q*std::cos(phi) 
		);

	// kappa = (1-z)k - zq
	double kappa = std::sqrt(
		SQR(1.0-z)*k*k + SQR(z)*q*q - 2.0*z*(1.0-z)*k*q*std::cos(phi)
		);
		
		

	// \kappa + q - p1 = (1-z)(k+q) - p1
	double kappa_p_q_m_p1 = std::sqrt(
		SQR(1.0-z)*(k*k + q*q + 2.0*k*q*std::cos(phi))
		+ p1*p1
		- 2.0*(1.0-z)*(p1*k*std::cos(theta) + p1*q*std::cos(theta-phi))
		);
	
	double p1_m_q = std::sqrt(
		p1*p1 + q*q - 2.0*p1*q*std::cos(theta-phi)
		);
	
	double result = par->N->S_k(p1, ya) * par->N->S_k(delta_m_p1, ya) / std::pow(2.0*M_PI, 4);
	result *= (SQR(z) + SQR(1.0-z)) * SQR(kappa_p_q_m_p1) / (SQR(kappa) * SQR(p1_m_q) );
	
	return p1*result;
	*/
}



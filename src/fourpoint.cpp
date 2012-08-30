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
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2012
 */

struct Inthelper_correction
{
    double pt1,pt2,ya;
    double phi;
    AmplitudeLib* N;
    size_t calln; int monte;
    CrossSection2* xs;
    bool finite_nc;

    double u1,u2,r,theta1,theta2,thetar,z;
};

double Inthelperf_correction(double* vec, size_t dim, void* p);
void pVegas_wrapper(double x[6], double f[1], void* par);
double NormalizeAngle(double phi);

bool cyrille=false; 
bool correction=true;
double IR_CUTOFF = LAMBDAQCD;

/*
 * Calculates the correction term/full integral over 6-dim. hypercube
 */

double CrossSection2::CorrectionTerm(double pt1, double pt2, double ya, double phi,
        double z)

{
    std::time_t start = std::time(NULL);
    //IR_CUTOFF = std::max(pt1, pt2);
    #ifdef USE_MPI
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if(id==0) 
    #endif
    cout << "# Starting MC integral, marquet: " << cyrille << ", correction: " << correction 
        << " finite-Nc " << FiniteNc() << " ir-cutoff " << IR_CUTOFF << " GeV " << endl;
    


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
    double minr=std::log(0.0001); double maxr=std::log(40);
    //double minr=0; double maxr=100;

       
    //cout << "-----------" << endl;
    double res=0;
    // Constants taken out from integral in order to reduse number of
    // products in integrand
    double constants;
    if (M_Q()>1e-4)
        constants = 8.0*SQR(M_PI)*SQR(z)*SQR(M_Q())/std::pow(2.0*M_PI, 6);
    else
        constants = 8.0*SQR(M_PI) / std::pow(2.0*M_PI, 6);

        Inthelper_correction helper;
            helper.pt1=pt1; helper.pt2=pt2; helper.phi=phi; helper.ya=ya;
            helper.N=N; helper.z=z; helper.calln=0; helper.monte=3;
            helper.xs=this; helper.finite_nc=FiniteNc();
             
        #ifndef USE_MPI
            double (*integrand)(double*,size_t,void*) = Inthelperf_correction;
            const gsl_rng_type *T = gsl_rng_default;
            //gsl_rng_env_setup();
            gsl_monte_function G = {integrand, 6, &helper};
            
            double result,abserr;
            double lower[6] = {minr, minr, minr, 0, 0, 0};
            double upper[6] = {maxr, maxr,maxr, 2.0*M_PI, 2.0*M_PI, 2.0*M_PI};
            gsl_rng *r = gsl_rng_alloc(T);
            gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (6);
         
            gsl_monte_vegas_integrate (&G, lower, upper, 6, mcintpoints/5, r, s,
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
                gsl_monte_vegas_integrate (&G, lower, upper, 6, mcintpoints, r, s,
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
                cout <<"# " << phi << " MC integral with " << mcintpoints << " points and " << i << " iterations took " << (std::time(NULL)-start)/60.0/60.0 <<" hours " << endl;
            }
       #endif   // USE_MPI



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
    vegas(reg, dim, pVegas_wrapper,
        0, mcintpoints/100, 5,  0,
        functions, 0, workers,
        estim, std_dev, chi2a, &helper);

      // refine the grid (init = 1) 
      // collect in additional accumulators (fcns = FUNCTIONS),
      time(&t);
    if (id==0)
        cout << "# 2nd round at " << std::ctime(&t);  
    vegas(reg, dim, pVegas_wrapper,
            1, mcintpoints/10, 5, NPRN_RESULT,
            functions, 0, workers,
            estim, std_dev, chi2a, &helper);
      // final sample, inheriting previous results (init = 2)
    time(&t);
    if (id==0)
        cout << "# 3rd round at " << std::ctime(&t) ;  
    vegas(reg, dim, pVegas_wrapper,
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
		return -100;
    }
	else if (std::abs(std_dev[0]/estim[0])<0.01) // Can reduce number of mcintpoints
		mcintpoints /= 2.0;
    

    
    
    res = estim[0]*constants;
    #endif

    return res;
}

// smooth step function
double step(double x) { return 1.0/(1.0+std::exp(-8.0*x)); }


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
    
    bool finite_nc = par->finite_nc;
    
    if (isnan(u1) or isnan(u2))
    {
        cerr << "u1/u2 nan, WTF???" << endl;
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
    



    double result = 0;

	
    if (cyrille)
    {
		//cerr << "Do we really want to calculate this?\n";
		// Stupid naive Large-Nc 6-point-function
		if (!finite_nc)
			result = std::cos(-pt1_dot_r - pt2_dot_r - pt2_dot_u2 + pt2_dot_u1)
                * s_u1*s_u2*s_r;
        
        
        result += (
            - std::cos(-pt1_dot_r - pt2_dot_r + pt2_dot_u1 + (1.0-z)*pt1_dot_u2-z*pt2_dot_u2)
             * s_u1
            - std::cos(-pt1_dot_r - pt2_dot_r - pt2_dot_u2 - (1.0-z)*pt1_dot_u1 + z*pt2_dot_u1)
             * s_u2
            + std::cos(-pt1_dot_r - pt2_dot_r + (1.0-z)*(pt1_dot_u2 - pt1_dot_u1)
                - z*(pt2_dot_u2 - pt2_dot_u1) )
            ) * s_r;
            
    }
    
    double r_m_u1=0, s_r_m_u1=0, r_p_u2=0, s_r_p_u2=0, r_m_u1_p_u2=0, s_r_m_u1_p_u2=0;
    if (finite_nc or correction)
    {      
        r_m_u1 = std::sqrt(SQR(r) + SQR(u1) - 2.0*u1_dot_r);
        s_r_m_u1 = N->S(r_m_u1, y);
        // r + u'
        r_p_u2 = std::sqrt(SQR(r) + SQR(u2) + 2.0*u2_dot_r);
        s_r_p_u2 = N->S(r_p_u2, y);
        // r - u + u' = r + (u'-u)
        r_m_u1_p_u2 = std::sqrt( SQR(r) + (SQR(u1) + SQR(u2) - 2.0*u1_dot_u2)
                    + 2.0*(u2_dot_r - u1_dot_r) );
        s_r_m_u1_p_u2 = N->S(r_m_u1_p_u2, y);
        
    }
    
    
    if (correction)
    {
		
		// separate 1/Nc^2 terms which contain DPS contribution, non-DPS
		// terms are computed analytically in xs.cpp
		double nc_suppress=0;
		if (finite_nc)
		{
				nc_suppress = s_r * cos(-pt1_dot_r - pt2_dot_r + pt1_dot_u2 - pt1_dot_u1)
					* (- 1.0 +  step(u1-1.0/IR_CUTOFF)*step(u2-1.0/IR_CUTOFF));
			nc_suppress *= 1.0/SQR(Nc);
		}
		
		double s6=0,quadrupole=0;
		// Quadrupole
		// Finite-Nc
		if (finite_nc)
		{
			// F(b,x;x',b')
			double f_b1x1x2b2 = 1.0/Cf * std::log( s_r_m_u1 * s_r_p_u2  / (s_r_m_u1_p_u2 * s_r ) );
			// F(b,x';x,b')
			double f_b1x2x1b2 = 1.0/Cf * std::log( s_u1 * s_u2 / (s_r_m_u1_p_u2 * s_r ) );
			// F(b,b';x',x)
			double f_b1b2x2x1 = 1.0/Cf * std::log( s_r_m_u1 * s_r_p_u2 / (s_u1 * s_u2) );
			
			// \Delta in Ref.
			double f_delta = SQR( f_b1x2x1b2 ) + 4.0/SQR(Nc) * f_b1x1x2b2 * f_b1b2x2x1 ;
			double sqrt_fdelta = std::sqrt(f_delta);
			//if (isnan(f_delta) or isnan(f_b1x1x2b2) or isnan(f_b1x2x1b2))
			//	cout << "sqrtdelta nan: u1 " << u1 << " u2 " << u2 << " r " << r << " r-u1 " << r_m_u1 << " r+u2 " << r_p_u2 << " r-u1+u2 " << r_m_u1_p_u2 << endl;
			quadrupole = s_u1 * s_u2 *
				( ( (sqrt_fdelta + f_b1x2x1b2 )/2.0 - f_b1x1x2b2)/sqrt_fdelta
					* std::exp(Nc/4.0 * sqrt_fdelta)
				+((sqrt_fdelta - f_b1x2x1b2 )/2.0 + f_b1x1x2b2)/sqrt_fdelta
					* std::exp(-Nc/4.0 * sqrt_fdelta) 
				)
				* std::exp( -Nc/4.0*f_b1x2x1b2 + 1.0/(2.0*Nc)*f_b1x1x2b2);
			
			if (isinf(quadrupole))
				cerr << "Inf quadrupole, r " << r << " u1 " << u1 << " u2 " << u2 << endl;
			if (isnan(quadrupole)) 
			{ 
				// DPS limit, this limit is shown analytically in notes
				// u1,u2>>r: Q -> S(x-x')S(b-b') also in finite-nc 
				if (s_u1==0 and s_u2==0 and s_r>0 )
					quadrupole = s_r*s_r_m_u1_p_u2;
				else   // All other nan limits vanish, right?
					quadrupole = 0;	
			}
			s6 = s_r * quadrupole;
		}
		else
		{
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
				
			s6 = s_r*s_u1*s_u2 - s_r * loglog * (s_u1 * s_u2 - s_r * s_r_m_u1_p_u2); 		
		}
		
		//result = s6*cos( -pt1_dot_r - pt2_dot_r - pt2_dot_u2 + pt2_dot_u1 ) + nc_suppress;
		result = s6;
		
		// Remove non-Nc-supressed DPS contribution
		//if (u1 > 1.0/cutoff and u2>1.0/cutoff)
		if (!isnan(quadrupole) and !isinf(quadrupole))
		{
			result -= s_r*s_r*s_r_m_u1_p_u2 
				* step(u1-1.0/IR_CUTOFF)* step(u2-1.0/IR_CUTOFF);	// \theta(u1-1/CUTOFF)\tehta(u2-1/CUTOFF)			
		}
		// Subtract this here in order to force the integrand to vanish in UV,
		// contribution calculated/compensated analytically in xs.cpp
		double subtract_coeff = 1.0;
		if (finite_nc)
			subtract_coeff = 1.0 - 1.0/SQR(Nc);
		result -= subtract_coeff * s_r*s_u1*s_u2;
		result*=std::cos(-pt1_dot_r - pt2_dot_r - pt2_dot_u2 + pt2_dot_u1);
		
		// add nc-suppressed contribution
		result += nc_suppress;
    
    }

    /*
    if (correction and !finite_nc)
    {
		double f1 = s_r_m_u1 * s_r_p_u2 / (s_r_m_u1_p_u2 * s_r);
        double f2 = s_u1 * s_u2 / (s_r_m_u1_p_u2 * s_r);
        double loglog=0;     
        loglog = std::log(f1) / std::log(f2);
        // If u1,u2>>r loglog->1
        if (s_u1==0 and s_u2==0 and s_r>0 )
            loglog=1;
        // All other NaN limits vanish
            
        if (isnan(loglog) )
            loglog=0;
        if (isinf(loglog))
        {
           //cout << "inf at u1 " << u1 << " u2 " << u2 << " r " << r << " f1 " << f1 << " f2 " << f2 <<
           // " s " << s_r*(s_u1 * s_u2 - s_r * s_r_m_u1_p_u2) << " dpfd " << s_r* s_r * s_r_m_u1_p_u2 << endl;
            
            loglog=0;
        }
            
        double tmpresult =  s_r * loglog * (s_u1 * s_u2 - s_r * s_r_m_u1_p_u2); 
            
        // Remove DPS, in large-\Nc it is S(b-b')S(x-x')^2
        if (u1 > 1.0/cutoff and u2>1.0/cutoff)
            tmpresult += s_r * s_r * s_r_m_u1_p_u2;
         
        tmpresult *= cos( -pt1_dot_r - pt2_dot_r - pt2_dot_u2 + pt2_dot_u1 );
        result -=tmpresult;
    }*/

    
    
    // Wave function product
    // \phi^*(u2) \phi(u) summed over spins and polarization
    // simple resul in massless z=0 case
    //result *= 8.0*SQR(M_PI) * u1_dot_u2 / ( SQR(u1)*SQR(u2) );
    // with masses and z!=0 a bit more complicated:
    // 8pi^2 m^2 z^2 [ K_1(mzu) K_1(mzu') u.u'/(uu') [1+(1-z)^2]
    //                + z^2 K_0(mzu)K_0(mzu') ]
    // In order to optimize, factor 8*pi^2*m^2*z^2
    // is outside the integral
    // factor 1/k^+ is thrown away as it cancels eventually
    double M_Q = par->xs->M_Q();
    double bessel_u1[2]; double bessel_u2[2];
    if (M_Q > 1e-4)
    {
        gsl_sf_bessel_Kn_array(0,1,par->z*M_Q*u1, bessel_u1);
        gsl_sf_bessel_Kn_array(0,1,par->z*M_Q*u2, bessel_u2);
        double wf =  (1.0+SQR(1.0-par->z))*u1_dot_u2 /(u1*u2)
                    * bessel_u1[1] // gsl_sf_bessel_K1(par->z*M_Q*u1)
                    * bessel_u2[1] //gsl_sf_bessel_K1(par->z*M_Q*u2)
                  + SQR(par->z)* bessel_u1[0]//gsl_sf_bessel_K0(par->z*M_Q*u1)
                    * bessel_u2[0]; //gsl_sf_bessel_K0(par->z*M_Q*u2)  ;
        
        
        result *= wf;
        if (isnan(result)){
            cout <<"nan, u1bessel: " << par->z*M_Q*u1 << " u2 " << par->z*M_Q*u2 << " u1 " 
            << u1 << " u2 " << u2 << " z " << par->z << " m " << M_Q 
            << " r " << r << " theta1 " << theta1 << " " << LINEINFO << endl;
            //exit(1);
        }
    }
    else
    {
        double wf = u1_dot_u2 / (SQR(u1)*SQR(u2)) * (1.0+SQR(1.0-par->z));    // 8\pi^2 is outside the integral
        result *= wf;
    }



    if (isnan(result)){
        cerr << "NAN! " << LINEINFO << "\n"; exit(1); }

    result *= u1*u2*r;   // d^2 u_1 d^2 u2 d^2 r = du_1 du_2 dr u_1 u_2 u_3 d\theta_i


       
    result *= u1*u2*r;  // multiply by e^(ln u1) e^(ln u2) e^(ln r) due to the
                        // change of variables (Jacobian)

    //cout << vec[0] << " " << vec[1] << " " << vec[2] << " " << result << endl;



    return result;
}




void pVegas_wrapper(double x[6], double f[1], void* par)
{
    f[0] = Inthelperf_correction(x, 6, par);
}



void CrossSection2::SetMCIntPoints(unsigned long long points)
{
    mcintpoints=points; mcintpoints_orig = points;
}



const int RINTPOINTS = 4;
const int THETAINTPOINTS = 2;
const double MINLNR = std::log(1e-3);
const double MAXLNR = std::log(500);

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

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl;
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

    if (status)
        cerr << "Integral failed! Result " << result << " relerr " << std::abs(abserr/result)
        << " at " << LINEINFO << endl;
    
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



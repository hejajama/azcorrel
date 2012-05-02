/*
 * NAME
 *   pvegas_mpi.c
 *   Parallel version of G.P.Lepage's VEGAS-algorithm.
 *
 * SYNOPSIS
 *   void vegas(double regn[], int ndim, void (*fxn)(double x[], double f[]),
 *              int init, unsigned long long ncall, int itmx, int nprn,
 *              int fcns, int pdim, int wrks,
 *              double tgral[], double sd[], double chi2a[]);
 *
 *     regn[]: array specifying the region to be integrated, 2*ndim entries
 *     ndim: dimensionality of space
 *     (*fxn)(x[],f[]): pointer to function to be evaluated (must be MT-safe!)
 *     init: initialization level (start with 0, then 1, later 2)
 *     ncall: number of samples points per iteration
 *     itmx: number of iterations in one call
 *     nprn: bit field, see constants NPRN_* below
 *     fcns: actual number of integrands (<=FNMX),
 *           if > 1 additional function accumulators are used
 *     pdim: dimension of parallel space, 0==autosense, higher=="manual tuning"
 *     wrks: number of parallel working units
 *     tgral[]: pointer to estimate of result (maybe array)
 *     sd[]: pointer to estimate of standard deviation (maybe array)
 *     chi2a[]: pointer to chi-squared over ndim (maybe array)
 *
 *   The calling program must make sure that the MPI-environment has
 *   been set up properly. (This cannot be done from within a
 *   subprogram.) Please call MPI_Init(&argc,&argv) explicitly before
 *   calling vegas().
 *
 * DESCRIPTION
 *   pvegas is a parallel farmer-worker implementation of the popular
 *   VEGAS-algorithm. It splits up some dimensions (the first
 *   ndim_par ones) into separate chunks and then lets each worker
 *   evaluate one chunk of these dimensions (and all the remaining
 *   dimensions). The random-numbers from gfsr are parallelized as
 *   as well.
 *
 *   This is the version for MPI. It uses the same technique as the
 *   multi-threaded version (pvegas.c) for decomposition. It should
 *   compile on all MPI-based systems. Please feel free to contact
 *   <Richard.Kreckel@Uni-Mainz.DE> if you encounter any
 *   implementation-specific problems.
 * 
 *   No external random number generator (RNG) needs to be supplied. A
 *   shift register-generator (SR) is implemented which is initialized
 *   with the system-time. You can force reproducible results by
 *   defining REPRO to be some positive integer different from zero
 *   and sticking to that value. All versions of vegas provided by the
 *   author should return the same numerical values if the same values
 *   for REPRO are used.
 *
 *   Note that the RNG is guaranteed to work properly if your
 *   architecture adheres to the LP64-model, i.e. ints are 32 bit
 *   long. Another, more crucial, assumption is that chars are 8 bit
 *   long. If you are on some strange hardware you are well advised to
 *   check the random numbers manually by consulting the supplied
 *   sample program vegastest.c!
 *
 *   This version may differ considerably from other implementations
 *   found on the net with regards to the arguments passed to it. The
 *   design goal was to have a uniform interface across all versions
 *   of vegas supplied by the author (nvegas.c, pvegas.c,
 *   pvegas_mpi.c) and to make optimal use of parallel facilities.
 *   Consult vegas.h and the samples to make proper use of the
 *   interface.
 *
 * AUTHOR
 * Richard Kreckel, ThEP, Univ. Mainz, December 1997 - October 2000 */

/*
 * Changed by H.M. 01/2012-04/2012: 
 * - user can pass a void pointer to parameter struct to the integrand
 * - number of calls is an unsigned long long int variable, not double
 */

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "vegas.h"

#define TINY 1.0e-68         /* small, since we are in double-precision      */
#define MXWORK 511           /* no more than MXWORK+1 CPUs, please           */
#define NDIM_PAR 0           /* 0 = auto-determination of parallel-space     */
#define REPRO 0              /* 0 = default, others used for comparison      */

static int mds;              /* ==1: statified smpl. (see README)  (mds)     */
int nd;                      /* slices in grid (c.f. NDMX)         (nd)      */
int ng;                      /*                                    (ng)      */
unsigned long long int npg;                     /* number of calls within bin         (npg), H.M: added unsigned long long   */
int ndim_par;                /* dimensionality of parallel space             */
int gndim;                   /* global copy of ndim                          */
double dxg;                  /*                                    (dxg)     */
double xnd;                  /*                                    (xnd)     */
double xJac;                 /* Jacobian of integration            (xjac)    */
typedef struct {
  double ti;                 /* sum for f over bins                (ti)      */
  double tsi;                /* sum for variances over bins        (tsi)     */
} binAccu;                   /* accumulator over bins / hypercubes...        */
binAccu Ab[FNMX];       /* ...one for each integrand                    */
double d[NDMX][MXDIM];       /*                                    (d[][])   */
double di[NDMX][MXDIM];      /* delta i                            (di[][])  */
double dx[MXDIM];            /* width of integration region        (dx[])    */
double gregn[MXDIM];         /* global copy of regn                          */
static double xi[MXDIM][NDMX];  /*                                 (xi[][])  */
int kgl[MXDIM-1];            /* first few dimensions are stored here         */
void (*p_fxn)(double [], double [], void* par);  /* copy of fxn              */
unsigned int gfsr_m[SR_P];   /* status-field for GFSR-generators             */
int gfsr_k;                  /* pointer into status-field                    */
static int gfsr_not_initialized = 1;  /* flag: need to initialize GFSR-field */
unsigned int rdum;           /* linear congruential counter in kw_rand()     */
double gfsr_norm;            /* will be set such that gfsr is normalized     */
int functions;               /* copy of (*ctl).fcns                          */
int p_rank;                  /* stores result of MPI_Comm_rank()             */
int p_size;                  /* stores result of MPI_Comm_size()             */
double result_buf[(2*FNMX+(2*MXDIM*NDMX))];  /* buffer to send results  */
int kgl_buf[MXDIM-1];        /* buffer for sending the chunk of hypercubes   */

/*
 * This routine is used by vegas. It rebins a vector of densities xi into
 * new bins defined by a vector r.
 */
void rebin(double rc, int nd, double r[], double xin[], double xi[])
{
  int i;
  int k = 0;
  double dr = 0.0;
  double xn = 0.0;
  double xo = 0.0;
  
  for (i=0; i<nd-1; i++) {
    while (rc > dr)
      dr += r[k++];
    if (k > 1) xo = xi[k-2];
    xn = xi[k-1];
    dr -= rc;
    xin[i] = xn-(xn-xo)*dr/r[k-1];
  }
  for (i=0; i<nd-1; i++) xi[i] = xin[i];
  xi[nd-1] = 1.0;
}

/*
 * gfsr produces the random numbers once the starting values have been
 * initialized using gfsr_init.
 */
double gfsr_rand(unsigned int w[], int *k)
{
  int j;
  
  (*k)++;
  if (*k >= SR_P) *k = 0;
  j = *k + SR_Q;
  if (j >= SR_P) j -= SR_P;
  w[*k] = w[*k] ^ w[j];
  return((double)w[*k] * gfsr_norm);
}

/*
 * Simple linear congruential generator used to initialize SR.
 * The multiplier and increment were provided by KOMA, Univ.Mainz.
 * (We may be abusing them a little, but I have checked that 
 * it has maximal periodicity which is what we want.)
 */
unsigned int kw_rand(void)
{
  rdum = (1812433253*rdum + 314159265);
  return (rdum);
}

/*
 * gfsr_init initializes the sequences using values from kw_rand.
 */
void gfsr_init(long seed)
{
  int i;
 
  if (p_rank==0) printf("# Initializing SR-sequences with seed %ld\n",seed);
  gfsr_norm = (double)(1/(pow(2.0,(long double)(8*sizeof(int))) - 1));
  rdum = (unsigned int)seed;
#if (REPRO == 0)
  for (i=0; i<SR_P*p_rank; i++) {
    kw_rand();                     /* advance to our entry-point */
  }
#endif
  for (i=0; i<SR_P; i++)
    gfsr_m[i] = kw_rand();           /* initialize starting values */
  gfsr_k = -1;                       /*         initialize pointer */
  gfsr_not_initialized = 0;          /* only effective if p_rank=0 */
}

/*
 * This function sends the results collected by one worker back to
 * the master. It first assembles a result_buf in order to send
 * everything in one single operation.
 */
void p_vegasrewrite(binAccu r_Ab[FNMX],
                    double r_d[NDMX][MXDIM], double r_di[NDMX][MXDIM])
{
  int i, j;
  
  /* assemble the send-buffer */
  for (j=0; j<functions; j++) {
    result_buf[j] = r_Ab[j].ti;
  } 
  for (j=0; j<functions; j++) {
    result_buf[j+functions] = r_Ab[j].tsi;
  } 
  for (j=0; j<gndim; j++) {
    for (i=0; i<nd; i++) 
      result_buf[2*functions + j*nd + i] = r_d[i][j];
    for (i=0; i<nd; i++) 
      result_buf[2*functions + gndim*nd + j*nd + i] = r_di[i][j];
  }
  MPI_Send(result_buf, (2*functions + 2*(gndim*nd)), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
}

/*
 * The implementation of a VEGAS-worker with its local variables
 * and its code.
 */
void* p_vegasloop(void *par)
{
  int i, j, k, h, l;
  int ia[MXDIM];             /*                                    (ia[])    */
  int kg[MXDIM];             /*                                    (kg[])    */
  typedef struct {
    double f2;               /* f squared                          (f2)      */
    double fb;               /* sum for f within bin               (fb)      */
    double f2b;              /* sum for f2 within bin              (f2b)     */
    unsigned long long npg;       /* number of calls within bin f != 0            */
  } pointAccu;               /* accumulator over points x within bins...     */
  pointAccu Ax[FNMX];   /* ...one for each integrand                    */
  double f[FNMX];       /* array passed into fxn for evaluation at x    */
  double x[MXDIM];           /* evaluation point                   (x[])     */
  double rc;                 /*                                    (rc)      */
  double wgt;                /* weight                             (wgt)     */
  double xn;                 /*                                    (xn)      */
  double xo;                 /*                                    (xo)      */
  MPI_Status recv_status;
#if (REPRO != 0)
  int kgl_rep[MXDIM-1];      /* copy of kgl, used only if if REPRO != 0      */
  unsigned long i_rep;       /* number of RNs to skip                        */
  int f_rep = 0;             /* a flag: have we iterated yet?                */
#endif
  
  for (i=0; i<gndim; i++) kg[i] = 1;
  
  /* What we need to get from master before we start looping over orthogonal 
   * space:
   * int kgl[0..ndim_par-1]
   * What we need to send to the master, who in turn has to ADD them to his 
   * values:
   * Ab[0..functions-1], d[0..nd-1][0..gndim-1], di[0..nd-1][0..gndim-1]
   */
  for (;;) {   /* Here the outer iteration begins, the one which is split up */
    MPI_Send(&p_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Recv(&kgl_buf, ndim_par, MPI_INT, 0, 0, MPI_COMM_WORLD, &recv_status);
#if (REPRO != 0)
    for (i=0; i<ndim_par; i++) kgl_rep[i] = kgl[i];
#endif
    for (i=0; i<ndim_par; i++) kgl[i] = kgl_buf[i];
#if (REPRO != 0)
    i_rep = 0;
    for (j=0; j<ndim_par; j++) {
      for (k=1,i=0; i<gndim-j-1; i++) k *= ng;
      i = (kgl[j]-kgl_rep[j]);
      if (j==ndim_par-1 && f_rep) i--;
      i_rep += k*npg*i*gndim;
    }
    if (kgl[0] != 0)
      while (i_rep-- > 0)
        gfsr_rand(gfsr_m,&gfsr_k);
#endif
    if (kgl[0]==0 &&  /* we are done already! */
        ndim_par>0)   /* may not make much sense, but makes things foolproof */
      {
        p_vegasrewrite(Ab,d,di);
#if (REPRO != 0)
	if (f_rep) i_rep = 0; else {
	  for (k=1,i=0; i<gndim-ndim_par; i++) k *= ng;
	  i_rep = k*npg*gndim;
	}
        for (j=0; j<ndim_par; j++) {
          for (k=1,i=0; i<gndim-j-1; i++) k *= ng;
          i = (ng-kgl_rep[j]);
          i_rep += k*i*npg*gndim;
        }
        while (i_rep-- > 0)
          gfsr_rand(gfsr_m,&gfsr_k);
#endif
        return;
      }
    
    for (i=0; i<ndim_par; i++) kg[i] = kgl[i];
    for (;;) {  /* ...and here the inner one, which is untouched */
      for (j=0; j<functions; j++) {
        Ax[j].fb = 0.0;
        Ax[j].f2b = 0.0;
        Ax[j].npg = 0;
      }
      for (k=0; k<npg; k++) {
        wgt = xJac;
        for (j=0; j<gndim; j++) {
          xn = (kg[j]-gfsr_rand(gfsr_m,&gfsr_k))*dxg+1.0;
          ia[j] = ((int)xn<NDMX) ? ((int)xn) : (NDMX);
          ia[j] = (ia[j]>1) ? (ia[j]) : (1);
          if (ia[j] > 1) {
            xo = xi[j][ia[j]-1]-xi[j][ia[j]-2];
            rc = xi[j][ia[j]-2]+(xn-ia[j])*xo;
          } else {
            xo = xi[j][ia[j]-1];
            rc = (xn-ia[j])*xo;
          }
          x[j] = gregn[j]+rc*dx[j];
          wgt *= xo*xnd;
        }
        (*p_fxn)(x,f, par);
        for (j=0; j<functions; j++) {
          if (f[j] != 0.0) ++Ax[j].npg;
          f[j] *= wgt;
          Ax[j].f2 = f[j]*f[j];
          Ax[j].fb += f[j];
          Ax[j].f2b += Ax[j].f2;
        }
        for (j=0; j<gndim; j++) {
          di[ia[j]-1][j] += f[0];
          if (mds >= 0) d[ia[j]-1][j] += Ax[0].f2;
        }
      }
      for (j=0; j<functions; j++) {
        Ax[j].f2b = sqrt(Ax[j].f2b*Ax[j].npg);
        Ax[j].f2b = (Ax[j].f2b-Ax[j].fb)*(Ax[j].f2b+Ax[j].fb);
        if (Ax[j].f2b <= 0.0) Ax[j].f2b = TINY;
        Ab[j].ti += Ax[j].fb;
        Ab[j].tsi += Ax[j].f2b;
      }
      if (mds < 0) {
        for (j=0; j<gndim; j++) d[ia[j]-1][j] += Ax[0].f2b;
      }
      for (h=gndim; h>ndim_par; h--) {
        kg[h-1] %= ng;
        if (++kg[h-1] != 1) break;
      }
      if (h < ndim_par+1) break;
    }           /* End of the inner iteration */
#if (REPRO != 0)
    f_rep = 1;
#endif
  }             /* End of the outer iteration */
}

/*
 * The routine pvegas to be called by the user. Parameters are just like
 * in vegas, except that workers specifies the degree of parallelization.
 */
void vegas(double regn[], int ndim, void (*fxn)(double x[], double f[], void* p),
           int init, unsigned long long ncall, int itmx, int nprn,
           int fcns, int pdim, int wrks,
           double tgral[], double sd[], double chi2a[], void* par)
{
  static int ndo;            /*                                    (ndo)     */
  int it;                    /* iteration counter                  (it)      */
  static int ittot;          /* iteration counter across init>1              */
  unsigned long long calls;              /* real total number of calls to fxn  (calls)   */
  double dv2g;               /*                                    (dv2g)    */
  double xin[NDMX];          /* aux. variable for rebinning        (xin[])   */
  typedef struct {
    double Wgt;              /* weight                             (wgt)     */
    double sWgt;             /* cumulative sum for weights         (swgt)    */
    double sChi;             /* cumulative sum for chi^2           (schi)    */
    double sInt;             /* cumulative sum for integral        (si)      */
  } iterAccu;                /* accumulator for stuff at end of iteration... */
  static iterAccu Ai[FNMX];  /* ...one for each integrand                    */
  double dt[MXDIM];          /*                                    (dt[])    */
  double r[NDMX];            /*                                    (r[])     */
  unsigned long long i, j, k;               /* counters                           (i, j, k) */
  int whodunit;              /* who returned the result to master?           */
  double rc;                 /*                                    (rc)      */
  double xn;                 /*                                    (xn)      */
  double xo;                 /*                                    (xo)      */
  double wgt;                /* weight                             (wgt)     */
  int wmax;                  /* for computing limit, if wrks too large       */
  MPI_Status recv_status;
#if (REPRO != 0)
  unsigned long i_rep;       /* number of RNs to skip                        */
#endif

  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p_size);
  wrks = (wrks<1) ? 1 : ((MXWORK<wrks) ? MXWORK : wrks);
  if (wrks>=p_size) wrks = p_size-1;
  gndim = ndim; 
  if (NDIM_PAR == 0)
    ndim_par = ndim/2;
  else
    ndim_par = NDIM_PAR;
  functions = (fcns<FNMX) ? (fcns) : (FNMX);
  
  for (i=0; i<ndim_par; i++) kgl[i] = 1;
  p_fxn = fxn;
  for (j=0; j<ndim; j++) gregn[j] = regn[j];
  
#if (REPRO != 0)
  if (gfsr_not_initialized) gfsr_init(REPRO);
#else
  if (gfsr_not_initialized) gfsr_init((long)time(NULL));
#endif
  if (init <= 0) {    /* entry for cold start        */
    mds = ndo = 1;    /* Careful with that Axe, Eugene!    */
                      /* mds=0 will trash parallelization. */
    for (j=0; j<ndim; j++) xi[j][0] = 1.0;
  }
  if (init <= 1) {    /* inherit the previous grid   */
    for (j=0; j<functions; j++) {
      Ai[j].sInt = 0.0;
      Ai[j].sWgt = 0.0;
      Ai[j].sChi = 0.0;
    }
    ittot = 1;
  }
  if (init <= 2) {    /* inherit grid and results    */
    nd = NDMX;
    ng = 1;
    if (mds) {
      ng = (int)pow(ncall/2.0+0.25,1.0/ndim);
      mds = 1;
      if ((2*ng-NDMX) >= 0) {
        mds = -1;
        npg = ng/NDMX+1;
        nd = ng/npg;
        ng = npg*nd;
      }
    }
    wmax = 1; 
    for (i=0; i<ndim_par; i++) wmax *= ng;
    if (wrks>wmax) wrks = wmax;
    for (k=1,i=0; i<ndim; i++) k *= ng;
    npg = (ncall/k>2) ? (ncall/k) : (2);
    calls = npg * k;
    dxg = 1.0/ng;

    for (dv2g=1,i=0; i<ndim; i++) dv2g *= dxg;    
    dv2g = calls*calls*dv2g*dv2g/npg/npg/(npg-1.0);
    xnd = nd;
    dxg *= xnd;
    xJac = 1.0/calls;
    for (j=0; j<ndim; j++) {
      dx[j] = regn[j+ndim]-regn[j];
      xJac *= dx[j];
    }
    if (nd != ndo) {
      for (i=0; i<(nd>ndo?nd:ndo); i++) r[i] = 1.0;
      for (j=0; j<ndim; j++) rebin(ndo/xnd,nd,r,xin,xi[j]);
      ndo = nd;
    }
    if (!p_rank && nprn & NPRN_INPUT) {
      printf("# %s:  ndim= %3d  ncall= %llu   %3d+1/%3d CPU(s)\n",
             "#  Input parameters for vegas",ndim,calls,wrks,p_size);
      printf("# % 28s  ittot=%5d  itmx=%5d    %5d^%1d hypercubes\n"," ",ittot,itmx,ng,ndim_par);
      printf("# %28s  nprn=0x%04x  ALPH=%5.2f\n"," ",nprn,ALPH);
      printf("# %28s  mds=%3d  nd=%4d%15s npg=%d\n"," ",mds,nd," ",npg);
      for (j=0; j<ndim; j++) {
        printf("# %30s xl[%2d]= %11.4g xu[%2d]= %11.4g\n",
               " ",j,regn[j],j,regn[j+ndim]);
      }
    }
  }
  
  for (it=ittot; it<=itmx+ittot-1; it++) {
    kgl[0] = 1;  /* kgl[0] is also used to determine if we are done yet */
    for (i=0; i<ndim_par; i++) 
      if (kgl[i]!=1) {
        kgl[i] = 1;
      }
    for (j=0; j<functions; j++) {
      Ab[j].ti = Ab[j].tsi = 0.0;
    }
    for (j=0; j<gndim; j++) {
      for (i=0; i<nd; i++) d[i][j] = di[i][j] = 0.0;
    }
    
    if (p_rank == 0) { /* ====master================================ */
      do {
        MPI_Recv(&whodunit, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recv_status);
        for (j=0; j<ndim_par; j++) kgl_buf[j] = kgl[j];
        MPI_Send(kgl_buf, ndim_par, MPI_INT, whodunit, 0, MPI_COMM_WORLD);
        for (j=ndim_par; j>=1; j--) {
          kgl[j-1] %= ng;
          if (++kgl[j-1] != 1) break;
        }
        if (j < 1) kgl[0] = 0;  /* we are done! */
      } while (kgl[0] != 0);
      /* Send a final kgl (zeroed) to everyone */
      for (j=1; j<=wrks; j++) {
        MPI_Send(&kgl, ndim_par, MPI_INT, j, 0, MPI_COMM_WORLD);
      }
      /* Now wait for all the results to come in */
      for (k=0; k<wrks; k++) {
	MPI_Recv(&result_buf, (2*functions + 2*(gndim*nd)), MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &recv_status);
	for (j=0; j<functions; j++)        /* disassemble result_buf */
	  Ab[j].ti += result_buf[j];
	for (j=0; j<functions; j++) 
	  Ab[j].tsi += result_buf[j+functions];
        for (j=0; j<gndim; j++) {
          for (i=0; i<nd; i++) 
	    d[i][j] += result_buf[2*functions + j*nd + i];
          for (i=0; i<nd; i++) 
	    di[i][j] += result_buf[2*functions + gndim*nd + j*nd + i];
        } 
      }
      for (j=0; j<functions; j++)     /* assemble the send-buffer */
	result_buf[j] = Ab[j].ti;
      for (j=0; j<functions; j++) 
	result_buf[j+functions] = Ab[j].tsi;
      for (j=0; j<gndim; j++) {
        for (i=0; i<nd; i++) 
	  result_buf[2*functions + j*nd + i] = d[i][j];
        for (i=0; i<nd; i++) 
	  result_buf[2*functions + gndim*nd + j*nd + i] = di[i][j];
      }
      for (j=1; j<p_size; j++) {
        MPI_Send(&result_buf, (2*functions + 2*(ndim*nd)), MPI_DOUBLE, j, 2, MPI_COMM_WORLD);
      }
    }
    else {  /* ====slave============================================ */
      if (p_rank<=wrks) p_vegasloop(par);
#if (REPRO != 0)       /* Now advance the remaining RNGs in case the */
      if (p_rank>wrks) {           /* next call runs on more CPUs */
        for (i_rep=0;i_rep<calls*ndim;i_rep++) 
          gfsr_rand(gfsr_m,&gfsr_k);
      }
#endif
      MPI_Recv(&result_buf, (2*functions + 2*(ndim*nd)), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &recv_status);
      for (j=0; j<functions; j++)        /* disassemble result_buf */
	Ab[j].ti = result_buf[j];
      for (j=0; j<functions; j++) 
	Ab[j].tsi = result_buf[j+functions];
      for (j=0; j<gndim; j++) {
        for (i=0; i<nd; i++) 
	  d[i][j] = result_buf[2*functions + j*nd + i];
        for (i=0; i<nd; i++) 
	  di[i][j] = result_buf[2*functions + gndim*nd + j*nd + i];
      }
      
    }
    for (j=0; j<functions; j++) {
      Ab[j].tsi *= dv2g;
      Ai[j].Wgt = 1.0/Ab[j].tsi;
      Ai[j].sInt += Ai[j].Wgt*Ab[j].ti;
      Ai[j].sChi += Ai[j].Wgt*Ab[j].ti*Ab[j].ti;
      Ai[j].sWgt += Ai[j].Wgt;
      tgral[j] = Ai[j].sInt/Ai[j].sWgt;
      chi2a[j] = (Ai[j].sChi-Ai[j].sInt*tgral[j])/(it-0.9999);
      if (chi2a[j] < 0.0) chi2a[j] = 0.0;
      sd[j] = sqrt(1.0/Ai[j].sWgt);
      Ab[j].tsi = sqrt(Ab[j].tsi);
    }
    wgt = Ai[0].Wgt;
    /* Changed by H.M: prints also relerr */
    if (!p_rank && nprn & NPRN_RESULT) {
      printf("# %s %3d : integral = %14.7g +/-  %9.2g  (relerr %9.2g)\n",
             " iteration no.",it,Ab[0].ti,Ab[0].tsi, fabs(Ab[0].tsi/Ab[0].ti));
      printf("# %s integral =%14.7g+/-%9.2g  chi^2/IT n = %9.2g\n",
             " all iterations:  ",tgral[0],sd[0],chi2a[0]);
    }
    if (!p_rank && nprn & NPRN_SECRES) {
      for (i=1; i<functions; i++) {
        printf("#   %4d%s%14.7g+/-%9.2g  chi^2/IT n = %9.2g\n",
               i,".additional integral= ",tgral[i],sd[i],chi2a[i]);
      }
    }
    if (!p_rank && nprn & (NPRN_GRID | NPRN_GRID_2 | NPRN_GRID_4 | NPRN_GRID_8)) {
      for (j=0; j<ndim; j++) {
        printf("# data for axis  %2d\n",j);
        printf("# %6s%13s%11s%13s%11s%13s\n", 
               "X","delta i","X","delta i","X","delta i");
        for (i=0; i<nd; i += 3) {
          for (k=0; k<3 && i+k<nd; k++) {
            printf(" %8.5f%12.4g    ",xi[j][i+k],di[i+k][j]);
          }
          printf("\n");
          if (nprn & NPRN_GRID_8) k = 3*(8-1);
          if (nprn & NPRN_GRID_4) k = 3*(4-1);
          if (nprn & NPRN_GRID_2) k = 3*(2-1);
          if (nprn & NPRN_GRID) k = 3*(1-1);
          i += k;
        }
      }
    }
    for (j=0; j<ndim; j++) {
      xo = d[0][j];
      xn = d[1][j];
      d[0][j] = (xo+xn)/2.0;
      dt[j] = d[0][j];
      for (i=1; i<nd-1; i++) {
        rc = xo+xn;
        xo = xn;
        xn = d[i+1][j];
        d[i][j] = (rc+xn)/3.0;
        dt[j] += d[i][j];
      }
      d[nd-1][j] = (xo+xn)/2.0;
      dt[j] += d[nd-1][j];
    }
    for (j=0; j<ndim; j++) {
      rc = 0.0;
      for (i=0; i<nd; i++) {
        if (d[i][j] < TINY) d[i][j] = TINY;
        r[i] = pow((1.0-d[i][j]/dt[j])/
		   (log(dt[j])-log(d[i][j])),ALPH);
        rc += r[i];
      }
      rebin(rc/xnd,nd,r,xin,xi[j]);
    }
  }
  ittot += itmx;

  return;
}

#undef SR_P
#undef SR_Q
#undef NDIM_PAR
#undef REPRO
#undef MXWORK
#undef TINY

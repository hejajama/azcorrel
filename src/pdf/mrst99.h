// mrst99 class header file
// Comments regarding the C++ version to: jeppe.andersen@durham.ac.uk
// Original comments from the FORTRAN version:
//  C****************************************************************C
//  C								     C
//  C     This is a package for the new **corrected** MRST parton    C
//  C     distributions. The format is similar to the previous       C
//  C     (1998) MRST series.                                        C
//  C								     C
//  C     NOTE: 7 new sets are added here, corresponding to shifting C
//  C     the small x HERA data up and down by 2.5%, and by varying  C
//  C     the charm and strange distributions, and by forcing a      C
//  C     larger d/u ratio at large x.                               C
//  C								     C
//  C     As before, x times the parton distribution is returned,    C
//  C     q is the scale in GeV, MSbar factorization is assumed,     C
//  C     and Lambda(MSbar,nf=4) is given below for each set.        C
//  C								     C
//  C     NAMING SCHEME:                                             C
//  C						                     C
//  C  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
//  C  ----  ---    -------             --------  -------   ------   C
//  C								     C
//  C  1     COR01  central gluon, a_s    300      0.1175   0.00524  C
//  C  2     COR02  higher gluon          300      0.1175   0.00497  C
//  C  3     COR03  lower gluon           300      0.1175   0.00398  C
//  C  4     COR04  lower a_s             229      0.1125   0.00585  C
//  C  5     COR05  higher a_s            383      0.1225   0.00384  C
//  C  6     COR06  quarks up             303.3    0.1178   0.00497  C
//  C  7     COR07  quarks down           290.3    0.1171   0.00593  C
//  C  8     COR08  strange up            300      0.1175   0.00524  C
//  C  9     COR09  strange down          300      0.1175   0.00524  C
//  C  10    C0R10  charm up              300      0.1175   0.00525  C
//  C  11    COR11  charm down            300      0.1175   0.00524  C
//  C  12    COR12  larger d/u            300      0.1175   0.00515  C
//  C						                     C
//  C      The corresponding grid files are called cor01.dat etc.    C
//  C							  	     C
//  C      The reference is:                                         C
//  C      A.D. Martin, R.G. Roberts, W.J. Stirling, R.S Thorne      C
//  C      Univ. Durham preprint DTP/99/64, hep-ph/9907231 (1999)    C
//  C                                                                C
//  C      Comments to : W.J.Stirling@durham.ac.uk                   C
//  C                                                                C
//  C								     C
//  C****************************************************************C

#ifndef _MRST_H_INCLUDED_
const int np=8;
const int nx=49;
const int nq=37;
const int ntenth=23;
const double xmin=1E-5;
const double xmax=1.0;
const double qsqmin=1.25;
const double qsqmax=1E7;

struct s_partoncontent {
  double upv,dnv,usea,dsea,str,chm,bot,glu;
};

class c_mrst99function {
 private:
  double xx[nx+1];
  double qq[nq+1];
  double n0[np+1];
  double f[np+1][nx+1][nq+2];
//    double g[np+1];
 public:
  int initialise(int mode);
  struct s_partoncontent update(double x, double qsq);
};

class c_mrst {
 private:
  // Note that we start with the first function as array element zero.
  class c_mrst99function function[12]; 
 public:
  double * table[8];
  struct s_partoncontent cont;
  void mrst99(double x,double q,int mode);
  // The constructor (initialises the functions)
  c_mrst();
};
#define _MRST_H_INCLUDED_
#endif

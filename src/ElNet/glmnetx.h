#pragma once

#include <armadillo>
#include <vector>

extern "C" {
  /*
    c dense predictor matrix:
    c
    c call elnet(ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
    c            intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
    c
    c x(no,ni) = predictor data matrix flat file (overwritten)
    c
    c
    c other inputs:
    c
    c   ka = algorithm flag
    c      ka=1 => covariance updating algorithm
    c      ka=2 => naive algorithm
    c   parm = penalty member index (0 <= parm <= 1)
    c        = 0.0 => ridge
    c        = 1.0 => lasso
    c   no = number of observations
    c   ni = number of predictor variables
    c   y(no) = response vector (overwritten)
    c   w(no)= observation weights (overwritten)
    c   jd(jd(1)+1) = predictor variable deletion flag
    c      jd(1) = 0  => use all variables
    c      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
    c   vp(ni) = relative penalties for each predictor variable
    c      vp(j) = 0 => jth variable unpenalized
    c   cl(2,ni) = interval constraints on coefficient values (overwritten)
    c      cl(1,j) = lower bound for jth coefficient value (<= 0.0)
    c      cl(2,j) = upper bound for jth coefficient value (>= 0.0)
    c   ne = maximum number of variables allowed to enter largest model
    c        (stopping criterion)
    c   nx = maximum number of variables allowed to enter all models
    c        along path (memory allocation, nx > ne).
    c   nlam = (maximum) number of lamda values
    c   flmin = user control of lamda values (>=0)
    c      flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
    c      flmin >= 1.0 => use supplied lamda values (see below)
    c   ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
    c   thr = convergence threshold for each lamda solution.
    c      iterations stop when the maximum reduction in the criterion value
    c      as a result of each parameter update over a single pass
    c      is less than thr times the null criterion value.
    c      (suggested value, thr=1.0e-5)
    c   isd = predictor variable standarization flag:
    c      isd = 0 => regression on original predictor variables
    c      isd = 1 => regression on standardized predictor variables
    c      Note: output solutions always reference original
    c            variables locations and scales.
    c   intr = intercept flag
    c      intr = 0/1 => don't/do include intercept in model
    c   maxit = maximum allowed number of passes over the data for all lambda
    c      values (suggested values, maxit = 100000)
    c
    c output:
    c
    c   lmu = actual number of lamda values (solutions)
    c   a0(lmu) = intercept values for each solution
    c   ca(nx,lmu) = compressed coefficient values for each solution
    c   ia(nx) = pointers to compressed coefficients
    c   nin(lmu) = number of compressed coefficients for each solution
    c   rsq(lmu) = R**2 values for each solution
    c   alm(lmu) = lamda values corresponding to each solution
    c   nlp = actual number of passes over the data for all lamda values
    c   jerr = error flag:
    c      jerr  = 0 => no error
    c      jerr > 0 => fatal error - no output returned
    c         jerr < 7777 => memory allocation error
    c         jerr = 7777 => all used predictors have zero variance
    c         jerr = 10000 => maxval(vp) <= 0.0
    C      jerr < 0 => non fatal error - partial output:
    c         Solutions for larger lamdas (1:(k-1)) returned.
    c         jerr = -k => convergence for kth lamda value not reached
    c            after maxit (see above) iterations.
    c         jerr = -10000-k => number of non zero coefficients along path
    c            exceeds nx (see above) at kth lamda value.
  */
  void elnet_(int* ka, double* parm, int* no, int* ni, double* x, double* y,
              double* w, int* jd, double* vp, double* cl, int* ne, int* nx,
              int* nlam, double* flmin, double* ulam, double* thr, int* isd,
              int* intr, int* maxit, int* lmu, double* a0, double* ca, int* ia,
              int* nin, double* rsq, double* alm, int* nlp, int* jerr);
}

struct lambda_interp_t {
  arma::uvec left;
  arma::uvec right;
  arma::vec frac;
};

struct glm_cv_t {
  arma::vec cvm;
  arma::vec cvsd;
  double lambda_min;
  double lambda_sd;
};

class glm {
 public:
  // Base data
  arma::mat X;
  arma::vec Y;
  // Input to FORTRAN code
  int nlam; 
  int nobs;
  int nvars;
  int ne;
  int nx;
  int jd;
  double thresh;
  int isd;
  int intr;
  int ka;
  int lmu;
  int maxit;
  double flmin;
  double alpha;
  arma::vec w;
  arma::vec vp;
  arma::vec cl;
  arma::vec ulam;
  // Output of FORTRAN code
  arma::vec a0;
  arma::vec ca;
  arma::vec rsq;
  arma::vec alm;
  arma::ivec ia;
  arma::ivec nin;
  int nlp = 0;
  int jerr = 0;
  // Beta data
  arma::mat beta;
  arma::uvec sort_order;
  bool beta_success;
  // Ctor
  glm(arma::mat X, arma::vec Y, int _nlam = 6, double _flmin = 0.3,
      double _alpha = 1, arma::vec _ulam = arma::vec(1, arma::fill::zeros));
  // Calculate and return beta
  void calculate_beta();
  // Predict a new set of data
  arma::mat predict(arma::mat& newx);
  // Predict a new set of data given new vector of lambdas s
  arma::mat predict(arma::mat& newx, arma::vec& s);
  // k fold cross validation
  glm_cv_t k_fold_cv(arma::uword k);
};

glm glmnet(arma::mat&, arma::vec&, int nsteps, double fmin);

lambda_interp_t lambda_interp(arma::vec lambda, arma::vec s);

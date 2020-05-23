//
// Seidr - Create and operate on gene crowd networks
// Copyright (C) 2016-2019 Bastian Schiffthaler <b.schiffthaler@gmail.com>
//
// This file is part of Seidr.
//
// Seidr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Seidr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Seidr.  If not, see <http://www.gnu.org/licenses/>.
//

// Seidr
#include <common.h>
#include <glmnet2.h>
// External
#include <armadillo>
#include <boost/numeric/conversion/cast.hpp>

extern "C"
{
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
  void elnet_(int64_t* ka,
              double* parm,
              int64_t* no,
              int64_t* ni,
              double* x,
              double* y,
              double* w,
              int64_t* jd,
              double* vp,
              double* cl,
              int64_t* ne,
              int64_t* nx,
              int64_t* nlam,
              double* flmin,
              double* ulam,
              double* thr,
              int64_t* isd,
              int64_t* intr,
              int64_t* maxit,
              int64_t* lmu,
              double* a0,
              double* ca,
              int64_t* ia,
              int64_t* nin,
              double* rsq,
              double* alm,
              int64_t* nlp,
              int64_t* jerr);
}

glm
glmnet(arma::mat X, arma::vec Y, int64_t nsteps, double fmin)
{

  // Setup one more lambda malue than requested
  int64_t nlam = nsteps + 1;
  // Setup all variables passed to FORTRAN
  int64_t nobs = X.n_rows;
  int64_t nvars = X.n_cols;
  int64_t ne = nvars + 1;
  int64_t nx = nvars;
  int64_t jd = 0;
  double thresh = 1e-7; // NOLINT
  int64_t isd = 1;
  int64_t intr = 1;
  int64_t ka = 1;
  int64_t maxit = 1e5; // NOLINT
  double flmin = fmin; // NOLINT
  double alpha = 1;
  double* yp = Y.memptr();
  double* xp = X.memptr();

  std::unique_ptr<double> wp(new double[nobs]());
  for (int64_t i = 0; i < nobs; i++) {
    wp.get()[i] = 1.0;
  }

  std::unique_ptr<double> vpp(new double[nvars]());
  for (int64_t i = 0; i < nvars; i++) {
    vpp.get()[i] = 1.0;
  }

  std::unique_ptr<double> clp(new double[2 * nvars]());
  for (int64_t i = 0; i < 2 * nvars; i++) {
    clp.get()[i] = (i % 2 == 0 ? -9.9e35 : 9.9e35); // NOLINT
  }

  std::unique_ptr<double> ulamp(new double[1]()); // ignored

  // Allocate memory to fill (by FORTRAN code)
  int64_t lmu = 0;

  std::unique_ptr<double> a0(new double[nlam]());
  std::unique_ptr<double> ca(new double[nx * nlam]());
  std::unique_ptr<int64_t> ia(new int64_t[nx]());
  std::unique_ptr<int64_t> nin(new int64_t[nlam]());
  std::unique_ptr<double> rsq(new double[nlam]());
  std::unique_ptr<double> alm(new double[nlam]());

  int64_t nlp = 0;
  int64_t jerr = 0;

  // Run ElNet optimization
  elnet_(&ka,
         &alpha,
         &nobs,
         &nvars,
         xp,
         yp,
         wp.get(),
         &jd,
         vpp.get(),
         clp.get(),
         &ne,
         &nx,
         &nlam,
         &flmin,
         ulamp.get(),
         &thresh,
         &isd,
         &intr,
         &maxit,
         &lmu,
         a0.get(),
         ca.get(),
         ia.get(),
         nin.get(),
         rsq.get(),
         alm.get(),
         &nlp,
         &jerr);

  // Create safe unsigned variables to compare to
  auto u_lmu = boost::numeric_cast<arma::uword>(lmu);
  auto u_nx = boost::numeric_cast<arma::uword>(nx);
  auto u_nlam =  boost::numeric_cast<arma::uword>(nlam);
  auto u_nvars = boost::numeric_cast<arma::uword>(nvars);

  // Maximum coefficients that will enter the model at
  // any lambda
  int64_t ninmax = 0;
  for (arma::uword i = 0; i < u_lmu; i++) {
    if (nin.get()[i] > ninmax) {
      ninmax = nin.get()[i];
    }
  }
  auto u_ninmax = boost::numeric_cast<arma::uword>(ninmax);

  // Setup return struct
  arma::mat b(u_nvars, u_nlam);
  b.zeros();
  glm ret{ b, arma::uvec(1), false }; // returned in case of error

  if (u_ninmax > 0 && u_lmu > 0) {
    // Zero out noise arising from precision issues
    for (arma::uword i = 0; i < u_nx * u_lmu; i++) {
      double absval = ca.get()[i] < 0 ? (ca.get()[i] * -1) : ca.get()[i];
      if (absval < 1e-7) { // NOLINT
        ca.get()[i] = 0;
      }
    }
    // Create matrix in [feature,nlambda] filled by
    // model coefficients
    arma::mat cx(u_nx, u_lmu, arma::fill::zeros);
    for (arma::uword i = 0; i < u_lmu; i++) {
      for (arma::uword j = 0; j < u_nx; j++) {
        cx(j, i) = ca.get()[(i * u_nx) + j];
      }
    }
    // Subset coefficient matrix to only use values
    // actually set by GLMNet
    cx = cx.rows(0, u_ninmax - 1);

    // For selection convenience, create an arma::uvec
    // from the row indices of coefficients and for the
    // lambda columns excluding the first
    arma::uvec arma_ia(u_ninmax, arma::fill::zeros);
    for (arma::uword i = 0; i < u_ninmax; i++) {
      arma_ia(i) = ia.get()[i] - 1;
    }
    arma::uvec tcols(cx.n_cols - 1, arma::fill::zeros);
    for (arma::uword i = 0; i < cx.n_cols - 1; i++) {
      tcols(i) = i;
    }

    // Create final full coefficient matrix
    arma::mat bx(u_nx, u_nlam - 1, arma::fill::zeros);
    bx.submat(arma_ia, tcols) = cx.cols(1, cx.n_cols - 1);
    b = bx;
    ret.beta = b;
    ret.sort_order = arma_ia;
    if (u_lmu == u_nlam) {
      ret.success = true;
    }
  }

  return ret;
}
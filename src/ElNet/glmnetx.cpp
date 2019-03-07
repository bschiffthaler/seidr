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
#include <glmnetx.h>
#include <common.h>
// External
#include <armadillo>
#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <memory>

using boost::numeric_cast;

glm glmnet(arma::mat& X, arma::vec& Y, int nsteps, double fmin)
{
  glm ret(X, Y, nsteps, fmin);
  return ret;
}

glm::glm(arma::mat _X, arma::vec _Y, int _nlam, double _flmin,
         double _alpha, arma::vec _ulam) :
  lmu(0),
  beta_success(false)
{
  X = std::move(_X);
  Y = std::move(_Y);
  // Setup one more lambda malue than requested
  nlam = _nlam;
  // Setup all variables passed to FORTRAN
  nobs = X.n_rows;
  nvars = X.n_cols;
  ne = nvars + 1;
  nx = nvars;
  jd = 0;
  thresh = 1e-7;
  isd = 1;
  intr = 1;
  ka = 1;
  maxit = 100000;
  flmin = _flmin;
  alpha = _alpha;

  // Need to copy as data is overwritten during elnet_()
  arma::mat x = X;
  arma::vec y = Y;

  double * xptr = x.memptr();
  double * yptr = y.memptr();

  auto wptr = new double[nobs]();
  for (int i = 0; i < nobs; i++)
  {
    wptr[i] = 1.0;
  }

  auto vpptr = new double[nvars]();
  for (int i = 0; i < nvars; i++)
  {
    vpptr[i] = 1.0;
  }

  auto clptr = new double[2 * nvars]();
  for (int i = 0; i < 2 * nvars; i++)
  {
    clptr[i] = (i % 2 == 0 ? -9.9e35 : 9.9e35);
  }

  auto ulamptr = new double[1];

  if (flmin >= 1)
  {
    ulamptr = new double[_ulam.n_elem];
    for (arma::uword i = 0; i < _ulam.n_elem; i++)
    {
      ulamptr[i] = _ulam(i);
    }
    nlam = numeric_cast<int>(_ulam.n_elem);
  }

  auto a0ptr = new double[nlam]();
  auto captr = new double[nx * nlam]();
  auto iaptr = new int[nx]();
  auto ninptr = new int[nlam]();
  auto rsqptr = new double[nlam]();
  auto almptr = new double[nlam]();

  // Run ElNet optimization
  elnet_(&ka, &alpha, &nobs, &nvars, xptr, yptr, wptr,
         &jd, vpptr, clptr, &ne, &nx, &nlam, &flmin, ulamptr, &thresh, &isd,
         &intr, &maxit, &lmu, a0ptr, captr, iaptr,
         ninptr, rsqptr, almptr, &nlp, &jerr);

  // Create arma vectors from aux memory while at the same time
  // truncating data to reflect actual number of lambdas (solutions)
  a0 = arma::vec(a0ptr, lmu);
  ca = arma::vec(captr, nx * lmu);
  ia = arma::ivec(iaptr, nx);
  nin = arma::ivec(ninptr, lmu);
  rsq = arma::vec(rsqptr, lmu);
  alm = arma::vec(almptr, lmu);

  // Fix lambda return
  if (nlam > 1 && flmin < 1)
  {
    arma::mat lalm = arma::log(alm);
    alm(0) = exp(2 * lalm(1) - lalm(2));
  }

  // Clean up
  delete[] wptr;
  delete[] vpptr;
  delete[] clptr;
  delete[] a0ptr;
  delete[] captr;
  delete[] iaptr;
  delete[] ninptr;
  delete[] rsqptr;
  delete[] almptr;
  delete[] ulamptr;

}

void glm::calculate_beta() {
  // Setup return struct
  beta = arma::mat(nvars, nlam, arma::fill::zeros);
  sort_order = arma::uvec(1);
  beta_success = false;
  // Create safe unsigned variables to compare to
  arma::uword u_lmu = 0, u_nx = 0, u_ninmax = 0, u_nlam = 0;

  try
  {
    u_lmu = boost::numeric_cast<arma::uword>(lmu);
    u_nx = boost::numeric_cast<arma::uword>(nx);
    u_nlam = boost::numeric_cast<arma::uword>(nlam);
    // Maximum coefficients that will enter the model at
    // any lambda
    int ninmax = arma::max(nin);
    u_ninmax = boost::numeric_cast<arma::uword>(ninmax);
  }
  catch (boost::numeric::bad_numeric_cast& e)
  {
    throw std::runtime_error(e.what());
  }

  if (u_ninmax > 0 && u_lmu > 0)
  {
    // Zero out noise arising from precision issues
    // ca( arma::find(ca < 1e-7) ).fill(0);

    // Create matrix in [feature,nlambda] filled by
    // model coefficients
    arma::mat cx(ca);
    cx.reshape(u_nx, u_lmu);
    // Subset coefficient matrix to only use values
    // actually set by GLMNet
    cx = cx.rows(0, u_ninmax - 1);
    // For selection convenience, create an arma::uvec
    // from the row indices of coefficients and for the
    // lambda columns excluding the first
    arma::uvec arma_ia(u_ninmax, arma::fill::zeros);
    for (arma::uword i = 0; i < u_ninmax; i++)
    {
      arma_ia(i) = ia(i) - 1;
    }
    // Create final full coefficient matrix
    // arma::mat bx(u_nx, u_nlam - 1, arma::fill::zeros);
    // bx.submat(arma_ia, tcols) = cx.cols(1, cx.n_cols - 1);
    beta.set_size(nx, lmu);
    beta.rows(arma_ia) = cx;
    sort_order = arma_ia;
    if (u_lmu == u_nlam)
    {
      beta_success = true;
    }
  }
}

arma::mat glm::predict(arma::mat& newx)
{
  calculate_beta();

  // Create new beta matrix with intercept as first row
  arma::mat nbeta(beta.n_rows + 1, beta.n_cols, arma::fill::zeros);
  nbeta.row(0) = a0.t();
  nbeta.rows(1, nbeta.n_rows - 1) = beta;

  // Create new data matrix with intercept = 1 as first column
  arma::mat icep_newx(newx.n_rows, newx.n_cols + 1,  arma::fill::zeros);
  icep_newx.col(0).fill(1);
  icep_newx.cols(1, icep_newx.n_cols - 1) = newx;

  return icep_newx * nbeta;
}

arma::mat glm::predict(arma::mat& newx, arma::vec& s)
{
  calculate_beta();

  // Create new beta matrix with intercept as first row
  arma::mat nbeta(beta.n_rows + 1, beta.n_cols, arma::fill::zeros);
  nbeta.row(0) = a0.t();
  nbeta.rows(1, nbeta.n_rows - 1) = beta;

  // Do linear interpolation for beta coefficients
  lambda_interp_t ip = lambda_interp(alm, s);

  arma::mat dl = arma::diagmat(ip.frac, 0);

  arma::vec fr = ip.frac;
  fr.transform([](double x) {return 1 - x;});
  arma::mat dr = arma::diagmat(fr, 0);

  nbeta = nbeta.cols(ip.left) * dl + nbeta.cols(ip.right) * dr;

  // Create new data matrix with intercept = 1 as first column
  arma::mat icep_newx(newx.n_rows, newx.n_cols + 1,  arma::fill::zeros);
  icep_newx.col(0).fill(1);
  icep_newx.cols(1, icep_newx.n_cols - 1) = newx;

  return icep_newx * nbeta;
}

lambda_interp_t lambda_interp(arma::vec lambda, arma::vec s)
{
  lambda_interp_t ret;
  if (lambda.size() == 1)
  {
    arma::uvec left(s.size(), arma::fill::ones);
    ret.left = left;
    arma::uvec right(s.size(), arma::fill::ones);
    ret.right = right;
    arma::vec frac(s.size(), arma::fill::ones);
    ret.frac = frac;
  }
  else
  {
    double maxlam = arma::max(lambda);
    double minlam = arma::min(lambda);

    s(arma::find(s > maxlam)).fill(maxlam);
    s(arma::find(s < minlam)).fill(minlam);
    arma::uword k = lambda.size();
    arma::vec sfrac = s;

    double x0 = lambda(0);
    double x1 = lambda(k - 1);

    sfrac.transform([&](double v) {
      return (x0 - v) / (x0 - x1);
    });

    lambda.transform([&](double v) {
      return (x0 - v) / (x0 - x1);
    });

    arma::vec y(lambda.size());
    arma::vec yy;

    for (arma::uword i = 0; i < lambda.size(); i++)
      y(i) = i;

    interp1(lambda, y, sfrac, yy);

    arma::uvec l = arma::conv_to<arma::uvec>::from(arma::floor(yy));
    arma::uvec r = arma::conv_to<arma::uvec>::from(arma::ceil(yy));

    for (arma::uword i = 0; i < sfrac.size(); i++)
    {
      if (l(i) == r(i))
      {
        sfrac(i) = 1;
      }
      else
      {
        sfrac(i) =
          ( sfrac(i) - lambda(r(i)) ) /
          ( lambda(l(i)) - lambda(r(i)) );
      }

    }

    ret.left = l;
    ret.right = r;
    ret.frac = sfrac;
  }

  return ret;
}

glm_cv_t glm::k_fold_cv(arma::uword k)
{
  // Create sampling vector
  arma::uvec foldid(X.n_rows);
  arma::uword i = 0, j = 0;
  while (i < X.n_rows)
  {
    foldid(i) = j;
    i++;
    j = (j == k - 1) ? 0 : (j + 1);
  }
  foldid = arma::shuffle(foldid);

  // Populate list of models
  arma::mat predmat(Y.n_elem, alm.n_elem, arma::fill::zeros);
  std::vector<glm> outlist;

  // Train CV models
  arma::vec minlams(k, arma::fill::zeros);
  for (arma::uword i = 0; i < k; i++)
  {
    arma::uvec nwhich = arma::find(foldid != i);
    arma::vec y_sub = Y.elem(nwhich);
    arma::mat x_sub = X.rows(nwhich);
    glm g(x_sub, y_sub, nlam, flmin, alpha);
    outlist.push_back(g);
    minlams(i) = arma::min(g.alm);
  }

  double mlam = arma::max(minlams);

  arma::uvec which_lam = arma::find(alm >= mlam);
  arma::uvec nlams(k, arma::fill::zeros);

  arma::vec s = alm(which_lam);

  arma::uvec snlam(which_lam.n_elem, arma::fill::zeros);
  for (arma::uword i = 0; i < which_lam.n_elem; i++)
  {
    snlam(i) = (i);
  }

  for (arma::uword i = 0; i < k; i++)
  {
    arma::uvec which = arma::find(foldid == i);
    arma::mat x_sub = X.rows(which);
    arma::mat pred = outlist[i].predict(x_sub, s);
    //outlist[i].alm.print();
    //std::cout << "----\n";
    predmat.submat(which, snlam) = pred;
  }

  //foldid.print();

  arma::vec cvm(which_lam.n_elem, arma::fill::zeros);
  arma::vec cvsd(which_lam.n_elem, arma::fill::zeros);

  for (arma::uword i = 0; i < which_lam.n_elem; i++)
  {
    arma::vec this_col = predmat.col(i);
    arma::vec mse_col = Y - this_col;
    mse_col.transform([](double x) {return x * x;});
    double mean = arma::mean(mse_col);
    //this_col.transform([&](double x){ return pow(x - mean, 2); });
    //double sd = arma::mean(this_col) / (N - 1);
    //sd = sqrt(sd);
    cvm(i) = mean;
    //cvsd(i) = sd;
  }

  double cvmin = arma::min(cvm);
  arma::uvec idmin = arma::find(cvm <= cvmin);
  double lambda_min = arma::max(alm(idmin));

  /*
  arma::vec semin = (cvm + cvsd);
  semin = semin(idmin);
  idmin = arma::find(cvm <= semin(0));
  double lambda_sd = arma::max(alm(idmin));*/

  glm_cv_t ret;
  ret.cvm = cvm;
  ret.cvsd = cvsd;
  ret.lambda_min = lambda_min;
  ret.lambda_sd = 0;

  return ret;
}


/*
int main(int argc, char**argv)
{

  arma::mat gm;
  gm.load("../../example/100_100_expr.tsv");
  arma::mat X = gm.submat(0,1,99,99);
  arma::vec Y = gm.col(0);
  Y = Y.subvec(0,99);
  std::cout << "Get GLM\n";

  glm ret = glm(X, Y, 6, 0.3, 0.3);
  std::cout << "Get beta\n";
  ret.calculate_beta();
  //ret.beta.print();
  //ret.alm.print();
  std::cout << "Predict1\n";
  arma::mat newx = gm.submat(0,1,9,99);
  std::cout << "Predict2\n";
  arma::mat p = ret.predict(newx);
  //p.print();

  glm_cv_t cv = ret.k_fold_cv(10);

  cv.cvm.print();

  arma::vec l{cv.lambda_min};

  l.print();

  glm best(X, Y, 6, 1.0, 0.3, l);
  best.calculate_beta();

  best.beta.print();

  return 0;
}


*/

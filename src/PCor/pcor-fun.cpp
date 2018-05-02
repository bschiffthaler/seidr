//seidr
#include <pcor-fun.h>
#include <BSlogger.h>
//external
#include <armadillo>
#include <iostream>
#include <cmath>
#include <fstream>

void fast_svd(arma::mat& U, arma::vec& s, arma::mat& V, arma::mat& m)
{

  if(m.n_rows * 2 < m.n_cols)
    {
      arma::mat B = m * m.t();
      
      arma::svd(U, s, V, B);
      double mrows = B.n_rows;
      double tol = mrows * arma::max(s) * arma::datum::eps;
      arma::uvec positive = arma::find(s > tol);
      arma::vec d = arma::sqrt(s(positive));
      
      U = U.cols(positive);
      
      V = m.t() * U * arma::diagmat(1 / d);
      s = d;
    }
  else if(m.n_rows > 2 * m.n_cols)
    {
      arma::mat B = m.t() * m;

      arma::svd(U, s, V, B);
      double mrows = B.n_rows;
      double tol = mrows * arma::max(s) * arma::datum::eps;
      arma::uvec positive = arma::find(s > tol);
      arma::vec d = arma::sqrt(s(positive));

      V = V.cols(positive);

      U = m * V * arma::diagmat(1 / d);
      s = d;
    }
  else
    {
      arma::svd(U, s, V, m);
      double mr = m.n_rows > m.n_cols ? m.n_rows : m.n_cols;
      double tol = mr * arma::max(s) * arma::datum::eps;
      arma::uvec positive = arma::find(s > tol);
      s = s(positive);
      U = U.cols(positive);
      V = V.cols(positive);
    }
}

arma::mat pcor(arma::mat& X)
{
  LOG_INIT_CLOG();
  arma::mat ret;
  log(LOG_INFO) << "Estimating lambda...\n";
  double lambda = estimate_lambda(X);
  log(LOG_INFO) << lambda << '\n';
  double nrows = X.n_rows;
  double h1 = nrows / (nrows - 1);
  double w = 1 / nrows;

  arma::mat U;
  arma::vec s;
  arma::mat V;

  log(LOG_INFO) << "Calculating partial correlations...\n";
  fast_svd(U, s, V, X);

  arma::uword m = s.n_elem;
  
  arma::mat UTWU = U.t() * (U * w);
  
  for(arma::uword i = 0; i < UTWU.n_rows; i++)
    {
      for(arma::uword j = 0; j < UTWU.n_cols; j++)
        {
          UTWU(i,j) *= s(j);
          UTWU(i,j) *= s(i);
        }
    }
    
  arma::mat C = (1 - lambda) * h1 * UTWU;

  C = (C + C.t()) / 2;
  
  if(lambda == 0)
    {
      //svdxs$v %*% tcrossprod( mpower(C, alpha), svdxs$v)
      ret = V * (arma::inv(C) * V.t());
    }
  else
    {
      //F = diag(m) - mpower(C/lambda + diag(m), alpha)
      //(diag(p) - svdxs$v %*% tcrossprod(F, svdxs$v))*(lambda)^alpha
      arma::mat F = arma::eye(m,m) -
        arma::inv(C / lambda + arma::eye(m,m));
      ret = (arma::eye(X.n_cols, X.n_cols) - (V * (F * V.t()))) *
        pow(lambda, -1);
    }

  ret *= -1;
  ret.diag() *= -1;

  arma::vec Is = arma::sqrt(1 / ret.diag());

  for(arma::uword i = 0; i < ret.n_rows; i++)
    for(arma::uword j = 0; j < ret.n_cols; j++)
      ret(i, j) *= Is(i);
  
  for(arma::uword i = 0; i < ret.n_rows; i++)
    for(arma::uword j = 0; j < ret.n_cols; j++)
      ret(i, j) *= Is(j);

  log(LOG_INFO) << "Writing data\n";
  
  return ret;
  
}

double estimate_lambda(arma::mat& X)
{
  double lambda = 0;
  double nrows = X.n_rows;
  
  double w2 = 1 / nrows;
  double h1w2 = w2/(1-w2);
  double sw = sqrt(w2);
  
  // Equivalent to xsw = sweep(xs, MARGIN=1, STATS=sw, FUN="*")
  // when weights are 1/nrow
  arma::mat xsw = X * sw;

  // SVD of positive values
  arma::mat U;
  arma::vec s;
  arma::mat V;

  fast_svd(U, s, V, xsw);
  
  // sE2R = sum(xsw*(sweep(xswsvd$u,2,xswsvd$d^3,'*')%*%t(xswsvd$v))) -
  // sum(colSums(xsw^2)^2)
  arma::vec d3 = arma::pow(s, 3);
  for(arma::uword i = 0; i < U.n_rows; i++)
    {
      for(arma::uword j = 0; j < U.n_cols; j++)
        {
          U(i,j) *= d3(j);
        }
    }

  arma::mat uvt = U * V.t();
  for(arma::uword i = 0; i < uvt.n_rows; i++)
    {
      for(arma::uword j = 0; j < uvt.n_cols; j++)
        {
          uvt(i,j) *= xsw(i,j);
        }
    }
  double sE2R = arma::accu(uvt);
  sE2R -= arma::accu(arma::pow(arma::sum(arma::pow(xsw, 2)), 2));


  //2*sum(xs2w[,(p-1):1] * t(apply(xs2w[,p:2, drop=FALSE],1,cumsum)))
  arma::mat xs2w = arma::pow(X, 2) * sw;
  arma::uvec xs2rhind(xs2w.n_cols - 1);
  arma::uvec xs2lhind(xs2w.n_cols - 1);
  
  arma::uword j = 0;
  for(arma::uword i = xs2w.n_cols - 1; i > 0; i--)
    xs2rhind(j++) = i;

  xs2lhind = xs2rhind - 1;
      
  arma::mat rhs = arma::cumsum(xs2w.cols(xs2rhind), 1);
  
  arma::mat lhs = xs2w.cols(xs2lhind);
  
  for(arma::uword i = 0; i < lhs.n_rows; i++)
    {
      for(arma::uword j = 0; j < lhs.n_cols; j++)
        {
          lhs(i,j) *= rhs(i,j);
        }
    }
  double sER2 = 2 * arma::accu(lhs);
  
  if(sE2R == 0)
    {
      lambda = 1;
    }
  else
    {
      double x = ((sER2 - sE2R) / sE2R) * h1w2;
      lambda = x > 1 ? 1 : x < 0 ? 0 : x;
    }
  
  return lambda;
}

void write_lm(arma::mat gm, std::string outfile, bool abs)
{
  std::ofstream ofs(outfile, std::ios::out);
  if(abs)
    gm = arma::abs(gm);
  if(! ofs)
    throw std::runtime_error("Could not write to file " + outfile);

  for(size_t i = 1; i < gm.n_cols; i++)
    {
      for(size_t j = 0; j < i; j++)
        {
          ofs << gm(i,j);
          ofs << (j == (i - 1) ? '\n' : '\t');
        }
    }

}

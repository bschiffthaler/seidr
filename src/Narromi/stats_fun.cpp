#include "stats_fun.h"
#include <armadillo>
#include <cmath>
#include "common.h"

/* MATLAB returns a 2x2 covariance table if the function is
invoked on two vectors, whereas armadillo's behaviour is
different. We calculate:

V = \left[\begin{matrix}
\frac 1n \sum_{i=1}^{n}x_i^2 - \bar x^2 & \frac 1n 
\sum_{i=1}^{n}x_iy_i - \bar x \bar y \\
\frac 1n \sum_{i=1}^{n}x_iy_i - \bar x \bar y & 
\frac 1n \sum_{i=1}^{n}y_i^2 - \bar y^2 \\
\end{matrix} \right]

as proposed here:
http://math.stackexchange.com/questions/524403/
covariance-matrix-for-2-vectors-with-elements-in-the-plane
*/
arma::mat mcov(arma::colvec x, arma::colvec y){
  arma::mat res(2,2);
  res.zeros();
  double n  = x.n_elem;
  double mx = arma::mean(x);
  double my = arma::mean(y);
  for(size_t i = 0; i < n; i++){
    res(0,0) += pow(x(i), 2) - pow(mx, 2);
    res(1,1) += pow(y(i), 2) - pow(my, 2);
    res(0,1) += x(i) * y(i) - mx * my;
    res(1,0) +=x(i) * y(i) - mx * my;
  }
  res.transform( [&n](double val) { return (val / (n - 1)); } );
  return res;
}

/* A simple calculator for mutual information using the
matrix determinants of covariance matrices. Straight from
the MATLAB implementation of narromi.
 */
double cmi(arma::colvec x, arma::colvec y){

  double c1   = arma::det(arma::cov(x));
  double c2   = arma::det(arma::cov(y));
  double c3   = arma::det(mcov(y, x));
  double cmiv = fabs(0.5 * log(c1 * c2 / c3));

  if(std::isnan(cmiv))
    cmiv = 1.0 ;

  return cmiv;
}

/* Helper function to calculate MI for a whole set of 
predictors and return a vector of the MIs.
 */
arma::vec fullMI(arma::mat X, arma::colvec Y){
  arma::vec res(X.n_cols);
  for(size_t f = 0; f < X.n_cols; f++){
    res(f) = cmi(X.col(f),Y);
  }
  return res;
}

/* Narromi estimates significance using a cumulative
distribution function. Here's a cpp implementation to
calculate phi as described here:

http://www.johndcook.com/blog/cpp_phi/

Precision is not optimal and two changes were made to 
enable the client to specify the mean mu and the 
standard deviation sigma.
 */
double phi(double x, double mu, double sigma)
{
  // constants
  double a1 =  0.254829592;
  double a2 = -0.284496736;
  double a3 =  1.421413741;
  double a4 = -1.453152027;
  double a5 =  1.061405429;
  double p  =  0.3275911;

  // Save the sign of x
  int sign = 1;
  if (x < mu)
    sign = -1;
  x = fabs(x-mu)/sqrt(2.0 * sigma * sigma);

  // A&S formula 7.1.26
  double t = 1.0/(1.0 + p*x);
  double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

  return 0.5*(1.0 + sign*y);
}
/* Wrapper function to calculate phi for an entire
vector of data.
 */
arma::vec normcdf(arma::vec x, double mu, double sigma){
  arma::vec res(x.n_elem);
  for(size_t i = 0; i < x.n_elem; i++){
    res(i) = phi(x(i), mu, sigma);
  }

  return res;
}


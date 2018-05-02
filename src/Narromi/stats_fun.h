#pragma once
#include <armadillo>
// covariance between two vectors
arma::mat mcov(arma::colvec x, arma::colvec y);
// calculate mutual information
double cmi(arma::colvec x, arma::colvec y);
// calculate MI for a whole set of predictors
arma::vec fullMI(arma::mat X, arma::colvec Y);
// cumulative distribution function(x)
double phi(double x, double mu, double sigma);
arma::vec normcdf(arma::vec x, double mu, double sigma);


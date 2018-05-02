#pragma once
#include <armadillo>
#include <string>

void fast_svd(arma::mat& U, arma::vec& s, arma::mat& V, arma::mat& m);
arma::mat pcor(arma::mat& X);
double estimate_lambda(arma::mat& X);
void write_lm(arma::mat gm, std::string outfile, bool abs);

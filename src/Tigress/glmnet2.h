#pragma once

#include <armadillo>

struct glm {
  arma::mat beta;
  arma::uvec sort_order;
  bool success;
};

glm glmnet(arma::mat, arma::vec, int nsteps, double fmin);


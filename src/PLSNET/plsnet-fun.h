#pragma once

#include <armadillo>
#include <vector>
#include <string>

struct simpls_t {
  arma::mat x;
  arma::mat y;
  arma::mat w;
};

struct plsreg_t {
  arma::mat p;
  arma::mat w;
};

simpls_t simpls(arma::mat& X, arma::vec& Y, arma::uword ncomp);
plsreg_t plsreg(arma::mat& X, arma::vec& Y, arma::uword ncomp);
arma::vec vip(arma::mat& X, arma::vec& Y, arma::uword ncomp);

void plsnet(arma::mat geneMatrix, std::vector<std::string> genes,
            std::vector<arma::uword> uvec, std::string outfile,
            int thread_id,
            arma::uword predictor_sample_size,
            arma::uword ensemble_size, arma::uword ncomp);

void plsnet_partial(arma::mat GM, std::vector<std::string> genes, size_t bs,
                    std::vector<std::string> targets, std::string outfile,
                    arma::uword predictor_sample_size,
                    arma::uword ensemble_size, arma::uword ncomp);

void plsnet_full(arma::mat GM, std::vector<std::string> genes, size_t bs,
                 std::string outfile,
                 arma::uword predictor_sample_size,
                 arma::uword ensemble_size, arma::uword ncomp);

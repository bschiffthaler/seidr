#pragma once

#include <vector>
#include <armadillo>
#include <string>

void el_ensemble(arma::mat geneMatrix, std::vector<std::string> genes,
                 std::vector<arma::uword> uvec, std::string outfile,
                 int thread_id, arma::uword min_sample_size,
                 arma::uword max_sample_size,
                 arma::uword predictor_sample_size_min,
                 arma::uword predictor_sample_size_max,
                 arma::uword ensemble_size, double alpha, double flmin,
                 arma::uword nlam);

void el_full(arma::mat GM, std::vector<std::string> genes, size_t bs,
             std::string outfile, arma::uword min_sample_size,
             arma::uword max_sample_size,
             arma::uword predictor_sample_size_min,
             arma::uword predictor_sample_size_max,
             arma::uword ensemble_size, double alpha, double flmin,
             arma::uword nlam);

void el_partial(arma::mat GM, std::vector<std::string> genes, size_t bs,
                std::vector<std::string> targets, std::string outfile,
                arma::uword min_sample_size,
                arma::uword max_sample_size,
                arma::uword predictor_sample_size_min,
                arma::uword predictor_sample_size_max,
                arma::uword ensemble_size,
                double alpha, double flmin,
                arma::uword nlam);

#pragma once

#include <vector>
#include <armadillo>
#include <string>
#include <svm.h>

void svm(arma::mat geneMatrix, std::vector<std::string> genes,
         std::vector<arma::uword> pred, std::string outfile,
         int thread_id, svm_parameter param,
         arma::uword min_sample_size, arma::uword max_sample_size,
         arma::uword predictor_sample_size_min,
         arma::uword predictor_sample_size_max,
         arma::uword ensemble_size);
void svm_full(arma::mat GM, std::vector<std::string> genes, size_t bs,
              std::string outfile, svm_parameter param,
              arma::uword min_sample_size, arma::uword max_sample_size,
              arma::uword predictor_sample_size_min,
              arma::uword predictor_sample_size_max,
              arma::uword ensemble_size);
void svm_partial(arma::mat GM, std::vector<std::string> genes, size_t bs,
                 std::vector<std::string> targets, std::string outfile,
                 svm_parameter param, arma::uword min_sample_size,
                 arma::uword max_sample_size,
                 arma::uword predictor_sample_size_min,
                 arma::uword predictor_sample_size_max,
                 arma::uword ensemble_size);

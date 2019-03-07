/* * Seidr - Create and operate on gene crowd networks
 * Copyright (C) 2016-2019 Bastian Schiffthaler <b.schiffthaler@gmail.com>
 *
 * This file is part of Seidr.
 *
 * Seidr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Seidr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Seidr.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <vector>
#include <armadillo>
#include <string>
#include <linear.h>

#define SVM_FULL 0
#define SVM_PARTIAL 1

#define LOG_NAME "llr-ensemble"

class seidr_mpi_svm;

struct seidr_llr_param_t
{
  std::string infile;
  std::string gene_file;
  std::string targets_file;
  bool do_scale = false;
  bool force = false;
  char row_delim = '\n';
  char field_delim = '\t';
  uint64_t bs;
  size_t mode = SVM_FULL;
  std::string outfile;
  std::string tempdir;
  arma::uword min_sample_size;
  arma::uword max_sample_size;
  arma::uword predictor_sample_size_min;
  arma::uword predictor_sample_size_max;
  arma::uword ensemble_size;
  unsigned int verbosity;
  parameter svparam;
  std::string solver;
  int nthreads;
};

void svm(const arma::mat& geneMatrix,
         const std::vector<std::string>& genes,
         const std::vector<arma::uword>& uvec,
         const std::string& tmpdir,
         parameter& param,
         const arma::uword& min_sample_size,
         const arma::uword& max_sample_size,
         const arma::uword& predictor_sample_size_min,
         const arma::uword& predictor_sample_size_max,
         const arma::uword& ensemble_size,
         seidr_mpi_svm * self);
void svm_full(const arma::mat& GM,
              const std::vector<std::string>& genes, 
              seidr_llr_param_t& param);
void svm_partial(const arma::mat& GM,
                 const std::vector<std::string>& genes,
                 const std::vector<std::string>& targets,
                 seidr_llr_param_t& param);

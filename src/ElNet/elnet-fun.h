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

#define EL_FULL 0
#define EL_PARTIAL 1

#define LOG_NAME "el-ensemble"

class seidr_mpi_elnet;

struct seidr_elnet_param_t
{
  std::string infile;
  std::string gene_file;
  std::string targets_file;
  bool do_scale = false;
  bool force = false;
  char row_delim = '\n';
  char field_delim = '\t';
  uint64_t bs;
  size_t mode = EL_FULL;
  std::string outfile;
  arma::uword min_sample_size;
  arma::uword max_sample_size;
  arma::uword predictor_sample_size_min;
  arma::uword predictor_sample_size_max;
  arma::uword ensemble_size;
  double alpha;
  double flmin;
  arma::uword nlam;
  unsigned verbosity;
  std::string tempdir;
  int nthreads;
};

void el_ensemble(const arma::mat& geneMatrix, 
                 const std::vector<std::string>& genes,
                 const std::vector<arma::uword>& uvec, 
                 const std::string& tmpdir,
                 const arma::uword min_sample_size,
                 const arma::uword max_sample_size,
                 const arma::uword predictor_sample_size_min,
                 const arma::uword predictor_sample_size_max,
                 const arma::uword ensemble_size, 
                 const double alpha, 
                 const double flmin,
                 const arma::uword nlam,
                 seidr_mpi_elnet * self);

void el_full(const arma::mat& GM, 
             const std::vector<std::string>& genes, 
             const seidr_elnet_param_t& param);

void el_partial(const arma::mat& GM, 
                const std::vector<std::string>& genes, 
                const std::vector<std::string>& targets, 
                const seidr_elnet_param_t& param);

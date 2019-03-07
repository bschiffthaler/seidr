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

#define TIGLM_FULL 0
#define TIGLM_PARTIAL 1

#define LOG_NAME "tigress"

class seidr_mpi_tigress;

struct seidr_tigress_param_t
{
  std::string infile;
  std::string gene_file;
  std::string targets_file;
  bool do_scale = false;
  bool force = false;
  char row_delim = '\n';
  char field_delim = '\t';
  uint64_t bs;
  size_t mode = TIGLM_FULL;
  std::string outfile;
  seidr_uword_t boot;
  double fmin;
  seidr_uword_t nsteps;
  std::string tempdir;
  unsigned verbosity;
  int nthreads;
};

std::vector<arma::uvec> shuffle(unsigned int);
void tiglm(const arma::mat& geneMatrix,
           const std::vector<std::string>& genes,
           const std::vector<arma::uword>& pred,
           const std::string& tmpdir,
           const seidr_uword_t& boot,
           const double& fmin,
           const seidr_uword_t& nsteps,
           seidr_mpi_tigress * self);
void tiglm_full(const arma::mat& gene_matrix,
                const std::vector<std::string>& genes,
                const seidr_tigress_param_t& param);
void tiglm_partial(const arma::mat& gene_matrix,
                   const std::vector<std::string>& genes,
                   const std::vector<std::string>& targets,
                   const seidr_tigress_param_t& param);

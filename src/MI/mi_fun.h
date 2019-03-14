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

#include <armadillo>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <cmath>
#include <string>

#define LOG_NAME "mi"

class aranode {
public:
  aranode(arma::uword a, double b) : i(a), v(b) {}
  arma::uword i;
  double v;
};

struct seidr_mi_param_t
{
  std::string infile;
  uint64_t bs;
  std::string outfile;
  std::string mi_file;
  size_t spline_order;
  size_t num_bins;
  std::string mode;
  bool force = false;
  std::string targets_file;
  std::string gene_file;
  bool use_existing = false;
  unsigned verbosity;
  char m = 0;
  std::string tempdir;
  int nthreads;
};

arma::vec knot_vector(size_t spline_order, size_t num_bins);
double percentile(arma::vec& data, size_t percentile);
double iqr(arma::vec& data);
double bin_width(arma::vec& data);
arma::uvec bin_count(const arma::mat& gm, size_t multiplier);
arma::vec to_z(arma::vec& x);
double basis_function(size_t i, size_t p, double t,
                      arma::vec& knot_vector, size_t num_bins);
void find_weights(const arma::mat& gm, arma::vec& knots, arma::mat& wm,
                  size_t spline_order, size_t num_bins, size_t i);
arma::vec hist1d(arma::vec& x, arma::vec& knots, arma::vec& weights,
                 size_t spline_order, size_t num_bins);
double log2d(double x);
double entropy1d(const arma::mat& gm, arma::vec& knots, arma::mat& wm,
                 size_t spline_order, size_t num_bins, size_t i);
void hist2d(arma::vec& x, arma::vec& y, arma::vec& knots,
            arma::vec& wx, arma::vec& wy, arma::mat& hist,
            size_t spline_order, size_t num_bins);
double entropy2d(const arma::mat& gm, arma::vec& knots,
                 arma::mat& wm, size_t spline_order,
                 size_t num_bins, size_t xi, size_t yi);
void mi_sub_matrix(const arma::mat& gm, size_t num_bins, size_t spline_order,
                   std::vector<arma::uword>& targets, 
                   const std::string& tmpfile);

void mi_full(const arma::mat& gm, 
             const std::vector<std::string>& genes,
             std::vector<std::string>& targets, 
             const seidr_mi_param_t& param);

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
#include <string>

struct seidr_pcor_param_t
{
  std::string infile;
  std::string outfile;
  std::string gene_file;
  std::string targets_file;
  bool do_scale = true;
  bool abs;
  bool force = false;
  unsigned verbosity;
};

void fast_svd(arma::mat& U, arma::vec& s, arma::mat& V, arma::mat& m);
arma::mat pcor(arma::mat& X);
double estimate_lambda(arma::mat& X);
void write_lm(arma::mat gm, std::string outfile, bool abs);

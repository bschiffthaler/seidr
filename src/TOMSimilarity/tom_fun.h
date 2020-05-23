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
#include <string>

constexpr unsigned TOM_DEF_VERBOSITY = 3;
constexpr uint64_t TOM_DEF_MAX_POW = 30;
constexpr double TOM_DEF_FIT = 0.8;

struct edge
{
  size_t i = 0;
  double w = 0;
  double r = 0;
};

enum tom_cor_t
{
  UNSIGNED,
  SIGNED,
  SIGNED_HYBRID
};

struct seidr_tom_param_t
{
  std::string infile;
  std::string outfile;
  std::string method;
  std::string gene_file;
  std::string targets_file;
  bool do_scale;
  bool abs;
  bool force = false;
  unsigned verbosity;
  std::string tom_type;
  uint64_t max_power;
  double sft_cutoff;
  uint64_t sft;
};

double
scale_free_fit(const arma::mat& gm,
               const uint64_t& p,
               const uint64_t& nbreaks,
               const tom_cor_t& cor);
arma::mat
tom_similarity(arma::mat& gm, const uint64_t& p, const tom_cor_t& cor);
bool
ascending(edge a, edge b);
void
rank_vector(std::vector<edge>& ev);
void
to_rank(arma::mat& GM);
void
write_lm(arma::mat& gm, const std::string& outfile, bool abs);
double
mad(const arma::vec& v);
arma::mat
bicor(const arma::mat& gm);

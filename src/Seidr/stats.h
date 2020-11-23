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
#include <string>
#include <vector>

constexpr double SEIDR_STATS_DEF_SAMPLES = 0.1;
constexpr double SEIDR_STATS_DEF_EV_TOL = 1e-8;

struct seidr_stats_param_t
{
  std::string infile;
  uint32_t tpos;
  bool trank;
  bool exact = false;
  uint32_t nsamples;
  std::string metrics;
  std::string tempdir;
  std::string tempfile;
  bool w_is_dist = false;
  bool directed = false;
  double ev_tol = SEIDR_STATS_DEF_EV_TOL;
};

int
stats(const std::vector<std::string>& args);
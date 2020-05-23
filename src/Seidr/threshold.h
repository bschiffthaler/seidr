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

constexpr uint16_t SEIDR_THRESHOLD_DEF_PREC = 8;
constexpr uint32_t SEIDR_THRESHOLD_DEF_NSTEP = 100;

struct seidr_threshold_param_t
{
  double min;
  double max;
  uint32_t nsteps;
  uint32_t tpos;
  std::string infile;
  bool trank = false;
  uint16_t precision;
  bool force;
  std::string outfile;
  int nthreads;
};

int
threshold(const std::vector<std::string>& args);
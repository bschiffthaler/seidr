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

constexpr uint16_t ASP_DEF_PRECISION = 8;
constexpr double ASP_DEF_MIN_WEIGHT = 1e-8;

struct seidr_asp_param_t
{
  uint32_t tpos;
  std::string infile;
  bool trank;
  uint16_t precision;
  bool invert;
  std::string outfile;
  bool force;
  bool absolute;
};

int
asp(const std::vector<std::string>& args);
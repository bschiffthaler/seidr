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

constexpr uint64_t SEIDR_NODE2VEC_DEF_WALK_LENGTH = 80;
constexpr uint64_t SEIDR_NODE2VEC_DEF_WALK_COUNT = 10;
constexpr double SEIDR_NODE2VEC_DEF_WALK_RETURN = 1;
constexpr double SEIDR_NODE2VEC_DEF_INOUT = 1;
constexpr uint32_t SEIDR_NODE2VEC_DEF_DIMENSION = 128;

struct seidr_node2vec_param_t
{
  std::string infile;
  uint32_t tpos;
  bool force;
  std::string outfile;
  bool directed;
  double walk_return;
  double walk_inout;
  uint64_t walk_length;
  uint64_t walk_count;
  uint32_t dimension;
};

int
node2vec(const std::vector<std::string>& args);

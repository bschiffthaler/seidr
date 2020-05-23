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

#include <string>
#include <vector>

constexpr uint32_t SEIDR_ROC_DEF_DATAP = 1000;

#pragma once
struct seidr_roc_param_t
{
  std::string gold;
  std::string gold_neg;
  std::string infile;
  uint32_t tpos;
  uint32_t nedg;
  std::string tfs;
  uint32_t datap;
  double fedg;
  bool force;
  std::string outfile;
  int nthreads;
  bool all;
};

int
roc(const std::vector<std::string>& args);

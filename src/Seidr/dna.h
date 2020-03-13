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

#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.hpp>

#include <vector>
#include <set>

struct seidr_dna_param_t
{
  std::string net1;
  std::string net2;
  std::string gene_to_pw;
  std::string target = "";
  std::string file_out;
  uint32_t tindex1;
  uint32_t tindex2;
  bool force = false;
};
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

#include <vector>
#include <string>

int view(int argc, char * argv[]);

struct seidr_param_view_t {
  std::string el_file;
  std::string ot_file;
  bool header;
  bool titles;
  seidr_score_t threshold;
  uint32_t tpos;
  uint32_t nlines;
  bool trank;
  bool binout;
  bool centrality;
  bool full;
  bool tags;
  bool print_supp;
  bool no_names;
  bool simple;
  std::string sc_delim;
  std::string delim;
  std::string filter;
  std::string nodes;
  std::vector<std::string> nodelist;
  bool force;
  std::string tempdir;
  std::string tempfile;
  std::string indexfile;
};

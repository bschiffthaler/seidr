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
#include <common.h>
#include <fstream>
#include <map>
#include <string>
#include <vector>

constexpr uint8_t SEIDR_CONVERT_DEF_PREC = 8;

typedef arma::mat mat_t;

struct seidr_conv_param_t
{
  std::string infile;
  std::string gene_file;
  std::string fill;
  std::string in_format;
  std::string in_sep;
  std::string out_format;
  std::string out_sep;
  uint16_t prec;
  bool force;
  std::string outfile;
};

void
parse_el(mat_t& m,
         std::map<std::string, size_t>& gm,
         std::istream& ifs,
         char sep);
void
parse_sm(mat_t& m, std::istream& ifs, char sep);
void
parse_ltri(mat_t& m, std::istream& ifs, bool diag, char sep);
void
parse_utri(mat_t& m, std::istream& ifs, bool diag, char sep);
void
parse_aracne(mat_t& m,
             std::map<std::string, size_t>& mp,
             std::istream& ifs,
             char sep);
void
write_el(mat_t& m, std::vector<std::string>& g, char sep);
void
write_sm(mat_t& m, char sep);
void
write_ltri(mat_t& m, bool diag, char sep);
void
write_utri(mat_t& m, bool diag, char sep);
void
write_aracne(mat_t& m, std::vector<std::string>& g, char sep);
int
convert(const std::vector<std::string>& args);

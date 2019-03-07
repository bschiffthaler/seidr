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
#include <BSlogger.h>

#include <vector>
#include <set>
#include <map>

struct seidr_compare_param_t
{
  std::string net1;
  std::string net2;
  std::string file_out;
  uint32_t tindex1;
  uint32_t tindex2;
  std::string trans;
  bool node_mode;
  bool force = false;
  std::string tempdir;
  std::string tempfile;
};

typedef std::pair<uint32_t, uint32_t> uint32_p;
typedef std::pair<std::string, std::string> string_p;

int compare(int argc, char * argv[]);

// Comparative index for comparison of networks
// Will create a set intersect of nodes and compute
// edges that **could** be in both datasets.
class cindex {
public:
  // ctor 
  cindex() = default;
  cindex(SeidrFileHeader& a, SeidrFileHeader& b);
  cindex(SeidrFileHeader& a, SeidrFileHeader& b, std::string& orth_dict);
  // Calculate the remapped indices of a new edge in the 
  // resulting network
  uint32_p reindex(SeidrFileHeader& h, SeidrFileEdge& e) { 
    return uint32_p(map_c.at(h.nodes[e.index.i]),
                    map_c.at(h.nodes[e.index.j])); 
  }
  // Get a vector of new node names after calculating the
  // set union of nodes of both networks
  std::vector<std::string>&  node_union(){ return node_uni; }
  // Check if the index needs to respect orthology mappings
  bool has_orth(){ return orth_mode; }
public:
  std::set<std::string> intersect;
  std::set<std::string> uni;
  std::set<std::string> diff_a;
  std::set<std::string> diff_b;
  std::map<std::string, uint32_t> map_a;
  std::map<std::string, uint32_t> map_b;
  std::map<std::string, uint32_t> map_c;
  std::map<std::string, std::string> dict_orth;
  std::map<size_t, string_p> dict_index;
  std::vector<std::string> node_uni;
  bool orth_mode = false;
};
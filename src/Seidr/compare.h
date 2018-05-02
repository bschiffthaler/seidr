#pragma once

#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>

#include <vector>
#include <set>
#include <map>

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
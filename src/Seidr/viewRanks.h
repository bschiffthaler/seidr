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
};
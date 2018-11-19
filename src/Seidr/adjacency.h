#pragma once
#include <string>

struct seidr_param_adjacency_t {
  std::string el_file;
  std::string out_file;
  std::string outfmt;
  std::string fill;
  std::vector<std::string> genes;
  bool force = false;
  uint32_t tpos;
  uint16_t prec;
  bool diag;
};

int adjacency(int argc, char * argv[]);

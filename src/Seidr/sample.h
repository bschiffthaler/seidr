#pragma once
#include <string>

struct seidr_param_random_t {
  std::string el_file;
  uint64_t nedges;
  double prop;
  uint32_t tpos;
  uint16_t prec;
  bool replacement;
};

int random(int argc, char * argv[]);

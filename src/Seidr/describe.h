#pragma once
#include <string>

struct seidr_param_describe_t {
  std::string el_file;
  uint32_t bins;
  bool parseable;
  std::string quantiles;
};

int describe(int argc, char * argv[]);

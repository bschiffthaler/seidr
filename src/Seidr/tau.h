#pragma once
#include <string>

struct seidr_param_tau_t {
  std::string network_a;
  std::string network_b;
  std::string out_file;
  std::string dict;
  bool force = false;
  uint32_t tpos_a;
  uint32_t tpos_b;
};

int tau(int argc, char * argv[]);

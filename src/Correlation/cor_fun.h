#pragma once

#include <armadillo>
#include <vector>

struct edge {
  size_t i = 0;
  double w = 0;
  double r = 0;
};

bool ascending(edge a, edge b);
void rank_vector(std::vector<edge>& ev);
void to_rank(arma::mat& GM);
void write_lm(arma::mat gm, std::string outfile, bool abs);

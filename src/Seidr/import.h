#pragma once

#include <vector>
#include <BSlogger.h>

#define PRINTING_MOD 100000000

class edge {
public:
  uint32_t i;
  uint32_t j;
  mutable seidr_score_t w = 0;
  mutable seidr_score_t r = 0;
  mutable char d = 0;
  friend bool operator==(const edge& lhs, const edge& rhs){
    return lhs.i == rhs.i && lhs.j == rhs.j;}
  friend bool operator!=(const edge& lhs, const edge& rhs){
    return ! (lhs == rhs);}
  friend bool operator<(const edge& lhs, const edge& rhs){
    if (lhs.i < rhs.i)
      return true;
    if (lhs.i == rhs.i)
      return lhs.j < rhs.j;
    return false;
  } 
};

struct reduced_edge {
  seidr_score_t w = 0;
  char d = 0;
};

struct shadow_t {
  shadow_t (void) : found(0), offset(0) {}
  uint8_t found = 0;
  uint64_t offset = 0;
};

class comb_hash {
public:
  size_t operator()(const edge& e) const
  {
    // Computes index in vector storing a lower triangular matrix
    // A perfect hash function if i != j and i > j, which is guaranteed fo our 
    // data
    size_t row = e.i + 1;
    size_t col = e.j + 1;
    return (row * (row - 1) / 2) - (row - col);
  }
};

class lt_map {
public:
  void insert(std::pair<uint32_t, uint32_t>&, reduced_edge&, bool&, bool&);
  std::vector<edge>& to_vec();
  void reserve(uint64_t n);
  comb_hash _hasher;
  std::vector<edge> _ev;
  std::vector<shadow_t> _shadow;
};

class read_logger {
public:
  read_logger(logger& l, size_t n = 100000) : _i(0), _n(n), _l(l) {};
  read_logger operator++(int);
  read_logger& operator++();
  size_t _i;
  size_t _n;
  logger& _l;
};

bool ascending(const edge a, const edge b);
bool descending(const edge a, const edge b);
bool abs_ascending(const edge a, const edge b);
bool abs_descending(const edge a, const edge b);
void rank_vector(std::vector<edge>&, bool, bool);

int import(int argc, char * argv[]);

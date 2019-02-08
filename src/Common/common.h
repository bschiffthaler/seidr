#pragma once

#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include <fstream>
#include <boost/tokenizer.hpp>

#define RANK_TYPE 0
#define AGGR_TYPE 1
#define FLT_EPSILON 0.000000001
#define ORDER_WEIGHT 0

#define SEIDR_MALLOC(type,n) (type *)malloc((n)*sizeof(type))

#define _XSTR(s) _STR(s)
#define _STR(s) #s

#if defined(SEIDR_PARALLEL_SORT) && defined(__ICC)
  #include <pstl/execution>
  #define SORT(start, end) std::sort(pstl::execution::par_unseq, start, end)
  #define SORTWCOMP(start, end, comp) std::sort(pstl::execution::par_unseq, start, end, comp)
#elif defined(SEIDR_PARALLEL_SORT) && defined(__GNUG__)
  #define SORT(start, end) __gnu_parallel::sort(start, end)
  #define SORTWCOMP(start, end, comp) __gnu_parallel::sort(start, end, comp)
#else
  #define SORT(start, end) std::sort(start, end)
  #define SORTWCOMP(start, end, comp) std::sort(start, end, comp)
#endif


#ifndef SEIDR_SCORE_DOUBLE
typedef float seidr_score_t;
#else
typedef double seidr_score_t;
#endif

#ifndef SEIDR_SCORE_DOUBLE
typedef arma::mat seidr_mat_t;
#else
typedef arma::fmat seidr_mat_t;
#endif

typedef arma::uword seidr_uword_t;
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

// Set difference for a single number
arma::uvec get_i(arma::uword ind, size_t s);

// Read gene names
std::vector<std::string> read_genes(std::string inputs, 
				    char row_delim = '\n', char field_delim = '\t');

bool is_gzip(std::string input);

bool almost_equal(seidr_score_t, seidr_score_t);

seidr_score_t unity_stand(seidr_score_t xmin, seidr_score_t xmax, seidr_score_t xi);

void scale(arma::mat& x);

static const std::string version = _XSTR(SEIDR_VERSION);

template <typename T>
inline void SWAP(T &x, T &y)
{
  T tmp = x;
  x = y; y = tmp;
}

std::vector<std::string> tokenize_delim(std::string nodes, std::string delim);

void merge_files(std::string outfile, std::string outfilebase,
                 std::string tempdir, bool targeted, int id,
                 std::vector<std::string>& genes);

bool any_const_expr(arma::mat& inp);

void verify_matrix(arma::mat& inp);

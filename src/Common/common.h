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

#include <Serialize.h>

#include <algorithm>
#include <armadillo>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

namespace po = boost::program_options;

#define RANK_TYPE 0
#define AGGR_TYPE 1
#define SFLT_EPSILON 0.000000001
#define ORDER_WEIGHT 0

#define SEIDR_MALLOC(type, n) (type*)malloc((n) * sizeof(type))

#define _XSTR(s) _STR(s)
#define _STR(s) #s
#define _AT __FILE__ ":" _XSTR(__LINE__)
#define _BUG(x)                                                                \
  throw std::runtime_error(                                                    \
    "Looks like you have encountered "                                         \
    "a bug.\nPlease report this message and steps to reproduce is to:\n"       \
    "https://github.com/bschiffthaler/seidr/issues\n"                          \
    "Error message: " x)

#if defined(SEIDR_PSTL)
#define SORT(start, end) std::sort(pstl::execution::par, start, end)
#define SORTWCOMP(start, end, comp)                                            \
  std::sort(pstl::execution::par, start, end, comp)
#define SET_NUM_PSTL_THREADS(x)                                                \
  tbb::global_control(tbb::global_control::max_allowed_parallelism, x);
#define GET_MAX_PSTL_THREADS() tbb::task_scheduler_init::default_num_threads()
#else
#define SORT(start, end) std::sort(start, end)
#define SORTWCOMP(start, end, comp) std::sort(start, end, comp)
#define SET_NUM_PSTL_THREADS(x)
#define GET_MAX_PSTL_THREADS() 1
#define INIT_TBB_CONTROL() int tbb_control = -1;
#endif

#ifdef SEIDR_WITH_MPI
#define SEIDR_MPI_BARRIER() MPI_Barrier(MPI_COMM_WORLD);
#define SEIDR_MPI_FINALIZE() MPI_Finalize();
#define SEIDR_MPI_INIT()                                                       \
  MPI_Init(&argc, &argv);                                                      \
  int rank = 0;                                                                \
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
#define SEIDR_MPI_BARRIER() nullptr;
#define SEIDR_MPI_FINALIZE() nullptr;
#define SEIDR_MPI_INIT() int rank = 0;
#endif

using seidr_score_t = double;
using seidr_mat_t = arma::mat;
using seidr_uword_t = arma::uword;
using tokenizer = boost::tokenizer<boost::char_separator<char>>;

#define SEIDR_DEFAULT_SSEQ                                                     \
  {                                                                            \
    3, 1, 4, 1, 5, 9, 2, 6, 5                                                  \
  }
#define SEIDR_DEFAULT_SEED 314159265

// Set difference for a single number
arma::uvec
get_i(arma::uword ind, size_t s);

// Read gene names
std::vector<std::string>
read_genes(const std::string& inputs,
           char row_delim = '\n',
           char field_delim = '\t');

// Check if value is in sorted range
template<typename T>
bool
in_sorted_range(const T& query, const std::vector<T>& subject)
{
  return std::binary_search(subject.cbegin(), subject.cend(), query);
}

bool
is_gzip(const std::string& input);

bool almost_equal(seidr_score_t, seidr_score_t);

seidr_score_t
unity_stand(seidr_score_t xmin, seidr_score_t xmax, seidr_score_t xi);

void
scale(arma::mat& x);

static const std::string version = _XSTR(SEIDR_VERSION);

template<typename T>
inline void
SWAP(T& x, T& y)
{
  T tmp = x;
  x = y;
  y = tmp;
}

std::vector<std::string>
tokenize_delim(const std::string& nodes, const std::string& delim);

void
merge_files(const std::string& outfile,
            std::string tempdir,
            bool targeted,
            int id,
            const std::vector<std::string>& genes);

bool
any_const_expr(const arma::mat& inp);

void
verify_matrix(const arma::mat& inp);
void
verify_matrix(const arma::mat& inp, const std::vector<std::string>& genes);

template<typename T>
void
assert_in_range(const T& num,
                const T& min,
                const T& max,
                const std::string& nam = "")
{
  if (num < min || num > max) {
    std::stringstream err;
    err << "Value";
    if (!nam.empty()) {
      err << " of " << nam << " ";
    }
    err << " (=" + std::to_string(num) << ") not in allowed range of ["
        << std::to_string(min) << "," << std::to_string(max) << "]";
    throw std::runtime_error(err.str());
  }
}

void
assert_mutually_exclusive(const po::variables_map& vm,
                          const std::vector<std::string>& targets);

std::string
str_join(const std::vector<std::string>& source, const std::string& delim);

void
make_tpos(uint32_t& tpos, const SeidrFileHeader& h);

template<typename T>
void
assert_arg_constraint(std::vector<T> allowed, T arg)
{
  auto hit = std::find(allowed.begin(), allowed.end(), arg);
  if (hit == allowed.end()) {
    std::stringstream ss;
    ss << "Argument '" << arg << "' is not allowed. Allowed values are ";
    ss << '[';
    for (uint64_t i = 0; i < allowed.size(); i++) {
      ss << "'" << allowed[i] << "'";
      if (i < (allowed.size() - 1)) {
        ss << ',';
      }
    }
    ss << ']';
    throw std::runtime_error(ss.str());
  }
}

uint64_t
get_mpi_nthread();

uint64_t
guess_batch_size(uint64_t const& set_size, uint64_t const& task_n);

// Prepare args to be passed to boost::po subprograms
inline std::vector<std::string>
shift_args(const int& argc, char* argv[])
{
  std::vector<std::string> args;
  if (argc < 2) {
    return args;
  }
  for (int i = 2; i < argc; i++) {
    args.push_back(argv[i]);
  }
  return args;
}

// Use for shared_ptr in case the object cannot be deleted
// e.g. std::cout
inline void
no_delete(void*)
{}
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

#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>

#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#define RANK_TYPE 0
#define AGGR_TYPE 1
#define SFLT_EPSILON 0.000000001
#define ORDER_WEIGHT 0

#define SEIDR_MALLOC(type,n) (type *)malloc((n)*sizeof(type))

#define _XSTR(s) _STR(s)
#define _STR(s) #s

#if defined(SEIDR_PSTL)
  #define SORT(start, end) \
    std::sort(pstl::execution::par, start, end)
  #define SORTWCOMP(start, end, comp) \
    std::sort(pstl::execution::par, start, end, comp)
  #define SET_NUM_PSTL_THREADS(x) \
    tbb::global_control(tbb::global_control::max_allowed_parallelism, x);
  #define GET_MAX_PSTL_THREADS() \
    tbb::task_scheduler_init::default_num_threads()
  #define INIT_TBB_CONTROL() tbb::global_control \
    tbb_control(tbb::global_control::max_allowed_parallelism, \
                tbb::task_scheduler_init::default_num_threads());
#else
  #define SORT(start, end) std::sort(start, end)
  #define SORTWCOMP(start, end, comp) std::sort(start, end, comp)
  #define SET_NUM_PSTL_THREADS(x)
  #define GET_MAX_PSTL_THREADS() 1
  #define INIT_TBB_CONTROL() int tbb_control = nullptr;
#endif

#ifdef SEIDR_WITH_MPI
    #define SEIDR_MPI_BARRIER() MPI_Barrier(MPI_COMM_WORLD);
    #define SEIDR_MPI_FINALIZE() MPI_Finalize();
    #define SEIDR_MPI_INIT() MPI_Init(&argc, &argv); int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    #define SEIDR_MPI_BARRIER() nullptr;
    #define SEIDR_MPI_FINALIZE() nullptr;
    #define SEIDR_MPI_INIT() int rank = 0;
#endif

typedef double seidr_score_t;

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
std::vector<std::string> read_genes(const std::string& inputs,
                                    char row_delim = '\n', char field_delim = '\t');

bool is_gzip(const std::string& input);

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

std::vector<std::string> tokenize_delim(const std::string& nodes,
                                        const std::string& delim);

void merge_files(const std::string& outfile,
                 std::string tempdir, bool targeted, int id,
                 const std::vector<std::string>& genes);

bool any_const_expr(arma::mat& inp);

void verify_matrix(arma::mat& inp);
void verify_matrix(arma::mat& inp, std::vector<std::string>& genes);

template<typename T>
void assert_in_range(const T& num, const T& min, const T& max,
                     const std::string& nam = "")
{
  if (num < min || num > max)
  {
    if (! nam.empty())
    {
      throw std::runtime_error("Value of " + nam + " (=" + std::to_string(num) +
                               ") not in allowed range of "
                               "[" + std::to_string(min) + "," +
                               std::to_string(max) + "]");
    }
    else
    {
      throw std::runtime_error("Value (=" + std::to_string(num) +
                               ") not in allowed range of "
                               "[" + std::to_string(min) + "," +
                               std::to_string(max) + "]");
    }
  }
}

void assert_mutually_exclusive(const po::variables_map& vm,
                               const std::vector<std::string> targets);

std::string str_join(const std::vector<std::string>& source,
                     const std::string& delim);

void make_tpos(uint32_t& tpos, const SeidrFileHeader& h);

template<typename T>
void assert_arg_constraint(std::vector<T> allowed, T arg)
{
  auto hit = std::find(allowed.begin(), allowed.end(), arg);
  if (hit == allowed.end())
  {
    std::stringstream ss;
    ss << "Argument '" << arg << "' is not allowed. Allowed values are ";
    ss << '[';
    for (uint64_t i = 0; i < allowed.size(); i++)
    {
      ss << "'" << allowed[i] << "'";
      if (i < (allowed.size() - 1))
      {
        ss << ',';
      }
    }
    ss << ']';
    throw std::runtime_error(ss.str());
  }
}

uint64_t get_mpi_nthread();
uint64_t guess_batch_size(uint64_t const & set_size,
                          uint64_t const & task_n);
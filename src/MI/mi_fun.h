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

#include <algorithm>
#include <armadillo>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <string>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#define LOG_NAME "mi"

constexpr unsigned MI_DEF_VERBOSITY = 3;
constexpr uint64_t MI_DEF_SPLINE = 3;
constexpr uint64_t MI_DEF_NUM_BINS = 0;
constexpr uint64_t MI_DEF_BS = 0;
constexpr uint64_t MI_MIN_BINS = 2;
constexpr uint64_t MI_WARN_BINS = 15;

class aranode
{
public:
  aranode(arma::uword a, double b)
    : i(a)
    , v(b)
  {}
  arma::uword i;
  double v;
};

class seidr_mpi_mi;

struct seidr_mi_param_t
{
  friend class boost::serialization::access;
  template<typename Archive>
  void serialize(Archive& ar, const unsigned int version) // NOLINT(clang-diagnostic-unused-parameter)
  {
    ar& BOOST_SERIALIZATION_NVP(infile);
    ar& BOOST_SERIALIZATION_NVP(gene_file);
    ar& BOOST_SERIALIZATION_NVP(targets_file);
    ar& BOOST_SERIALIZATION_NVP(mi_file);
    ar& BOOST_SERIALIZATION_NVP(force);
    ar& BOOST_SERIALIZATION_NVP(bs);
    ar& BOOST_SERIALIZATION_NVP(mode);
    ar& BOOST_SERIALIZATION_NVP(outfile);
    ar& BOOST_SERIALIZATION_NVP(tempdir);
    ar& BOOST_SERIALIZATION_NVP(verbosity);
    ar& BOOST_SERIALIZATION_NVP(nthreads);
    ar& BOOST_SERIALIZATION_NVP(use_existing);
    ar& BOOST_SERIALIZATION_NVP(spline_order);
    ar& BOOST_SERIALIZATION_NVP(num_bins);
    ar& BOOST_SERIALIZATION_NVP(m);
  }
  std::string infile;
  uint64_t bs;
  std::string outfile;
  std::string mi_file;
  size_t spline_order;
  size_t num_bins;
  std::string mode;
  bool force = false;
  std::string targets_file;
  std::string gene_file;
  bool use_existing = false;
  unsigned verbosity;
  char m = 0;
  std::string tempdir;
  bool resuming;
  int nthreads;
  std::vector<uint64_t> good_idx;
  std::string cmd_file;
};

arma::vec
knot_vector(size_t spline_order, size_t num_bins);
double
percentile(arma::vec& data, size_t percentile);
double
iqr(arma::vec& data);
double
bin_width(arma::vec& data);
arma::uvec
bin_count(const arma::mat& gm, size_t multiplier);
arma::vec
to_z(arma::vec& x);
double
basis_function(size_t i,
               size_t p,
               double t,
               arma::vec& knot_vector,
               size_t num_bins);
void
find_weights(const arma::mat& gm,
             arma::vec& knots,
             arma::mat& wm,
             size_t spline_order,
             size_t num_bins,
             size_t i);
arma::vec
hist1d(arma::vec& x,
       arma::vec& weights,
       size_t num_bins);
double
log2d(double x);
double
entropy1d(const arma::mat& gm,
          arma::mat& wm,
          size_t num_bins,
          size_t i);
void
hist2d(arma::vec& x,
       arma::vec& wx,
       arma::vec& wy,
       arma::mat& hist,
       size_t num_bins);
double
entropy2d(const arma::mat& gm,
          arma::mat& wm,
          size_t num_bins,
          size_t xi,
          size_t yi);
void
mi_sub_matrix(const arma::mat& gm,
              size_t num_bins,
              size_t spline_order,
              std::vector<arma::uword>& targets,
              const std::string& tmpdir,
              seidr_mpi_mi* self);

void
mi_full(const arma::mat& gm,
        const std::vector<std::string>& genes,
        std::vector<std::string>& targets,
        const seidr_mi_param_t& param);

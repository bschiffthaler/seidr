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

#include <armadillo>
#include <common.h>
#include <string>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#define EL_FULL 0
#define EL_PARTIAL 1

#define LOG_NAME "el-ensemble"

constexpr unsigned ELNET_DEF_VERBOSITY = 3;
constexpr seidr_uword_t ELNET_DEF_NLAM = 10;
constexpr double ELNET_DEF_FLMIN = 0.3;
constexpr double ELNET_DEF_ALPHA = 0.3;
constexpr uint64_t ELNET_DEF_BS = 0;
constexpr seidr_uword_t ELNET_DEF_ENSEMBLE = 1000;
constexpr seidr_uword_t ELNET_DEF_SAMPLE_SIZE_MIN = 0;
constexpr seidr_uword_t ELNET_DEF_SAMPLE_SIZE_MAX = 0;
constexpr seidr_uword_t ELNET_DEF_PREDICTOR_SIZE_MIN = 0;
constexpr seidr_uword_t ELNET_DEF_PREDICTOR_SIZE_MAX = 0;

constexpr double ELNET_DEF_MIN_RATIO = 0.2;
constexpr double ELNET_DEF_MAX_RATIO = 0.8;

constexpr arma::uword ELNET_DEF_MAX_K_CV = 10;
constexpr arma::uword ELNET_DEF_TOP_N_CONSIDERED = 20;

class seidr_mpi_elnet;

struct seidr_elnet_param_t
{
  friend class boost::serialization::access;
  template<typename Archive>
  // NOLINTNEXTLINE(clang-diagnostic-unused-parameter): version currently unused
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& BOOST_SERIALIZATION_NVP(infile);
    ar& BOOST_SERIALIZATION_NVP(gene_file);
    ar& BOOST_SERIALIZATION_NVP(targets_file);
    ar& BOOST_SERIALIZATION_NVP(do_scale);
    ar& BOOST_SERIALIZATION_NVP(force);
    ar& BOOST_SERIALIZATION_NVP(row_delim);
    ar& BOOST_SERIALIZATION_NVP(field_delim);
    ar& BOOST_SERIALIZATION_NVP(bs);
    ar& BOOST_SERIALIZATION_NVP(mode);
    ar& BOOST_SERIALIZATION_NVP(outfile);
    ar& BOOST_SERIALIZATION_NVP(tempdir);
    ar& BOOST_SERIALIZATION_NVP(min_sample_size);
    ar& BOOST_SERIALIZATION_NVP(max_sample_size);
    ar& BOOST_SERIALIZATION_NVP(predictor_sample_size_min);
    ar& BOOST_SERIALIZATION_NVP(predictor_sample_size_max);
    ar& BOOST_SERIALIZATION_NVP(ensemble_size);
    ar& BOOST_SERIALIZATION_NVP(verbosity);
    ar& BOOST_SERIALIZATION_NVP(alpha);
    ar& BOOST_SERIALIZATION_NVP(flmin);
    ar& BOOST_SERIALIZATION_NVP(nlam);
    ar& BOOST_SERIALIZATION_NVP(k);
    ar& BOOST_SERIALIZATION_NVP(nthreads);
  }
  std::string infile;
  std::string gene_file;
  std::string targets_file;
  bool do_scale = false;
  bool force = false;
  char row_delim = '\n';
  char field_delim = '\t';
  uint64_t bs;
  size_t mode = EL_FULL;
  std::string outfile;
  arma::uword min_sample_size;
  arma::uword max_sample_size;
  arma::uword predictor_sample_size_min;
  arma::uword predictor_sample_size_max;
  arma::uword ensemble_size;
  double alpha;
  double flmin;
  arma::uword nlam;
  uint16_t k;
  unsigned verbosity;
  std::string tempdir;
  bool resuming;
  int nthreads;
  std::vector<uint64_t> good_idx;
  std::string cmd_file;
};

void
el_ensemble(const arma::mat& geneMatrix,
            const std::vector<std::string>& genes,
            const std::vector<arma::uword>& uvec,
            const std::string& tmpdir,
            const arma::uword& min_sample_size,
            const arma::uword& max_sample_size,
            const arma::uword& predictor_sample_size_min,
            const arma::uword& predictor_sample_size_max,
            const arma::uword& ensemble_size,
            const double& alpha,
            const double& flmin,
            const arma::uword& nlam,
            const uint16_t& k,
            seidr_mpi_elnet* self);

void
el_full(const arma::mat& GM,
        const std::vector<std::string>& genes,
        const seidr_elnet_param_t& param);

void
el_partial(const arma::mat& GM,
           const std::vector<std::string>& genes,
           const std::vector<std::string>& targets,
           const seidr_elnet_param_t& param);

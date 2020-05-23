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
#include <string>
#include <svm.h>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#define SVM_FULL 0
#define SVM_PARTIAL 1

constexpr unsigned SVM_DEF_VERBOSITY = 3;
constexpr arma::uword SVM_DEF_ENSEMBLE = 1000;
constexpr int SVM_DEF_POLY_DEGREE = 3;
constexpr double SVM_DEF_EPS = 0.001;
constexpr double SVM_DEF_GAMMA = 0.01;
constexpr double SVM_DEF_COEF0 = 0.01;
constexpr double SVM_DEF_C = 1;
constexpr double SVM_DEF_P = 0.1;
constexpr double SVM_DEF_NU = 0.5;

#define LOG_NAME "svm-ensemble"

class seidr_mpi_svm;

struct seidr_svm_param_t
{
  friend class boost::serialization::access;
  template<typename Archive>
  void serialize(Archive& ar, const unsigned int version) // NOLINT: version
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
    ar& BOOST_SERIALIZATION_NVP(svm_type);
    ar& BOOST_SERIALIZATION_NVP(kernel);
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
  size_t mode = SVM_FULL;
  std::string outfile;
  arma::uword min_sample_size;
  arma::uword max_sample_size;
  arma::uword predictor_sample_size_min;
  arma::uword predictor_sample_size_max;
  arma::uword ensemble_size;
  unsigned int verbosity;
  std::string tempdir;
  svm_parameter svparam;
  std::string svm_type;
  std::string kernel;
  std::string cmd_file;
  bool resuming;
  int nthreads;
  std::vector<uint64_t> good_idx;
};

void
svm(const arma::mat& geneMatrix,
    const std::vector<arma::uword>& uvec,
    const std::string& tmpdir,
    svm_parameter& param,
    const arma::uword& min_sample_size,
    const arma::uword& max_sample_size,
    const arma::uword& predictor_sample_size_min,
    const arma::uword& predictor_sample_size_max,
    const arma::uword& ensemble_size,
    seidr_mpi_svm* self);
void
svm_full(const arma::mat& GM,
         const std::vector<std::string>& genes,
         seidr_svm_param_t& param);
void
svm_partial(const arma::mat& GM,
            const std::vector<std::string>& genes,
            const std::vector<std::string>& targets,
            seidr_svm_param_t& param);

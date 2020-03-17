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
#include <vector>
#include <string>


#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

#define PLSNET_FULL 0
#define PLSNET_PARTIAL 1

#define LOG_NAME "plsnet"

class seidr_mpi_plsnet;

struct seidr_plsnet_param_t
{
  friend class boost::serialization::access;
  template<typename Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & BOOST_SERIALIZATION_NVP(infile);
    ar & BOOST_SERIALIZATION_NVP(gene_file);
    ar & BOOST_SERIALIZATION_NVP(targets_file);
    ar & BOOST_SERIALIZATION_NVP(do_scale);
    ar & BOOST_SERIALIZATION_NVP(force);
    ar & BOOST_SERIALIZATION_NVP(row_delim);
    ar & BOOST_SERIALIZATION_NVP(field_delim);
    ar & BOOST_SERIALIZATION_NVP(bs);
    ar & BOOST_SERIALIZATION_NVP(mode);
    ar & BOOST_SERIALIZATION_NVP(outfile);
    ar & BOOST_SERIALIZATION_NVP(tempdir);
    ar & BOOST_SERIALIZATION_NVP(verbosity);
    ar & BOOST_SERIALIZATION_NVP(nthreads);
    ar & BOOST_SERIALIZATION_NVP(ensemble_size);
    ar & BOOST_SERIALIZATION_NVP(ncomp);
  }
  std::string infile;
  std::string gene_file;
  std::string targets_file;
  bool do_scale = false;
  bool force = false;
  char row_delim = '\n';
  char field_delim = '\t';
  uint64_t bs;
  size_t mode = PLSNET_FULL;
  std::string outfile;
  arma::uword predictor_sample_size;
  arma::uword ensemble_size;
  arma::uword ncomp;
  std::string tempdir;
  unsigned verbosity;
  bool resuming;
  int nthreads;
  std::vector<uint64_t> good_idx;
  std::string cmd_file;
};

struct simpls_t {
  arma::mat x;
  arma::mat y;
  arma::mat w;
};

struct plsreg_t {
  arma::mat p;
  arma::mat w;
};

simpls_t simpls(arma::mat& X, arma::vec& Y, arma::uword ncomp);
plsreg_t plsreg(arma::mat& X, arma::vec& Y, arma::uword ncomp);
arma::vec vip(arma::mat& X, arma::vec& Y, arma::uword ncomp);

void plsnet(const arma::mat& geneMatrix, 
            const std::vector<std::string>& genes,
            const std::vector<arma::uword>& uvec, 
            const std::string& tmpdir,
            const arma::uword& predictor_sample_size,
            const arma::uword& ensemble_size, 
            const arma::uword& ncomp,
            seidr_mpi_plsnet * self);

void plsnet_partial(const arma::mat& GM,
                    const std::vector<std::string>& genes,
                    const std::vector<std::string>& targets,
                    const seidr_plsnet_param_t& param);

void plsnet_full(const arma::mat& GM,
                 const std::vector<std::string>& genes,
                 const seidr_plsnet_param_t& param);

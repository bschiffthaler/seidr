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
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#define NARROMI_FULL 0
#define NARROMI_PARTIAL 1

#define LOG_NAME "narromi"

class seidr_mpi_narromi;

struct seidr_narromi_param_t
{
  friend class boost::serialization::access;
  template<typename Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& BOOST_SERIALIZATION_NVP(infile);
    ar& BOOST_SERIALIZATION_NVP(gene_file);
    ar& BOOST_SERIALIZATION_NVP(targets_file);
    ar& BOOST_SERIALIZATION_NVP(alpha);
    ar& BOOST_SERIALIZATION_NVP(force);
    ar& BOOST_SERIALIZATION_NVP(row_delim);
    ar& BOOST_SERIALIZATION_NVP(field_delim);
    ar& BOOST_SERIALIZATION_NVP(bs);
    ar& BOOST_SERIALIZATION_NVP(mode);
    ar& BOOST_SERIALIZATION_NVP(outfile);
    ar& BOOST_SERIALIZATION_NVP(tempdir);
    ar& BOOST_SERIALIZATION_NVP(verbosity);
    ar& BOOST_SERIALIZATION_NVP(nthreads);
    ar& BOOST_SERIALIZATION_NVP(t);
    ar& BOOST_SERIALIZATION_NVP(al);
  }
  std::string al;
  double alpha = 0.05;
  double t = 0.6;
  bool force = false;
  std::string infile;
  std::string gene_file;
  std::string targets_file;
  size_t mode = NARROMI_FULL;
  char row_delim = '\n';
  char field_delim = '\t';
  uint64_t bs;
  std::string outfile;
  std::string tempdir;
  unsigned verbosity;
  bool resuming;
  int nthreads;
  std::vector<uint64_t> good_idx;
  std::string cmd_file;
};

class NaResult
{
  arma::vec s;
  arma::vec n;
  arma::uvec r;
  arma::uword I;

public:
  NaResult(arma::vec, arma::vec, arma::uvec, arma::uword);
  arma::vec sig();
  arma::vec nv();
  arma::uvec re();
  arma::uword i();
};

std::pair<arma::vec, arma::vec>
reoptim(const arma::mat& X,
        const arma::rowvec& Y,
        const std::string& al,
        const double& alpha);

NaResult
narromi(const arma::mat& GM,
        const std::string& al,
        const double& alpha,
        const double& t,
        const arma::uword& i);

arma::vec
get_sig(arma::vec z);

void
full_narromi(const arma::mat& gene_matrix,
             const std::vector<std::string>& genes,
             const seidr_narromi_param_t& param);

void
partial_narromi(const arma::mat& gene_matrix,
                const std::vector<std::string>& genes,
                const std::vector<std::string>& targets,
                const seidr_narromi_param_t& param);

void
narromi_thread(const arma::mat& gene_matrix,
               const std::string& al,
               const double& alpha,
               const double& t,
               const std::vector<arma::uword>& ind,
               const std::vector<std::string>& genes,
               const std::string tmpdir,
               seidr_mpi_narromi* self);
